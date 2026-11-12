"""DuMuX/Python comparison on the heterogeneous 5 x 5 network with five top pores.

The liquid is fully saturated and the bottom DGF row is a fixed-pressure, fixed-composition
reservoir.  The top row evaporates through the same 1.6 mm stagnant boundary layer used by the
one-dimensional DuMuX reference case.  Equilibrium in-pore precipitation is opt-in but enabled
by the command-line case used for the comparison.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from ..materials import brine
from ..materials.salt import RHO_CR
from ..pnm import (DGF_DUMUX_5X5_VARTHROAT, PoreClogged, PoreLayer,
                   PorePrecipitation, boundary_layer_evaporation, network_step)
from .dumux1d import DELTA_BL, RH_ENV, T_DUMUX

X_RES = 0.08
C_RES = float(brine.mass_fraction(X_RES))
X_SOLUBILITY = float(brine.mole_fraction(brine.C_SAT))
P_ABS_REFERENCE = 1.0e5
TOP_FILE_INDICES = np.array([25, 26, 27, 28, 29], dtype=int)


def load_flux_replay(path):
    """Read a wide CSV with ``time_s`` and ``pore_25`` ... ``pore_29`` water fluxes [kg/s]."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = data.dtype.names or ()
    required = ("time_s",) + tuple(f"pore_{idx}" for idx in TOP_FILE_INDICES)
    missing = [name for name in required if name not in names]
    if missing:
        raise ValueError(f"{path}: missing replay columns {missing}")
    times = np.atleast_1d(data["time_s"]).astype(float)
    flux = np.column_stack([np.atleast_1d(data[f"pore_{idx}"]) for idx in TOP_FILE_INDICES])
    if np.any(np.diff(times) <= 0.0):
        raise ValueError(f"{path}: replay times must be strictly increasing")
    return times, flux.astype(float)


def _replay_at(replay, t):
    times, flux = replay
    return np.array([np.interp(t, times, flux[:, j]) for j in range(flux.shape[1])])


def _full_pore_state(layer, dynamic_values, reservoir_value):
    out = np.full(len(layer.dgf_vertices), reservoir_value, dtype=float)
    out[layer.file_vertex_index] = np.asarray(dynamic_values, float)
    return out


def _full_throat_state(layer, Q_lat, Q_res, p_full):
    n_edges = len(layer.dgf_edges)
    Q = np.zeros(n_edges)
    radius = np.asarray(layer.dgf_edges[:, 2], float).copy()
    for local_idx, file_edge_idx in enumerate(layer.lat_file_edge_index):
        Q[file_edge_idx] = Q_lat[local_idx]
        radius[file_edge_idx] = layer.lat[local_idx][2]
    for dynamic_idx, file_edge_idx in enumerate(layer.res_file_edge_index):
        if file_edge_idx < 0:
            continue
        a, b = layer.dgf_edges[file_edge_idx, :2].astype(int)
        a_is_res = a in set(layer.reservoir_file_indices.tolist())
        Q[file_edge_idx] = Q_res[dynamic_idx] if a_is_res else -Q_res[dynamic_idx]
        radius[file_edge_idx] = layer.rt_res[dynamic_idx]
    ends = layer.dgf_edges[:, :2].astype(int)
    p_mean = 0.5 * (p_full[ends[:, 0]] + p_full[ends[:, 1]])
    return Q, radius, p_mean


def run(dt=5.0, t_end=40000.0, output_interval=100.0, dgf=DGF_DUMUX_5X5_VARTHROAT,
        precipitation=True, replay_flux=None, max_solid_fraction=0.99):
    """Run the comparison case and return fixed-size checkpoint arrays.

    ``replay_flux`` is ``(times, per_top_pore_water_flux)`` and bypasses the Python vapour
    closure.  Positive flux means water leaves a top pore.  With ``None``, each top pore uses
    its own current water activity.
    """
    dt = float(dt)
    t_end = float(t_end)
    output_interval = float(output_interval)
    if dt <= 0.0 or t_end <= 0.0 or output_interval <= 0.0:
        raise ValueError("dt, t_end and output_interval must be positive")
    if abs(round(output_interval / dt) * dt - output_interval) > 1e-12 * output_interval:
        raise ValueError("output_interval must be an integer multiple of dt")

    layer = PoreLayer.from_dgf_dirichlet_boundary(dgf, boundary_axis=1, boundary_side="min")
    if not np.array_equal(layer.file_vertex_index[layer.surface], TOP_FILE_INDICES):
        raise ValueError("the selected DGF does not have the expected top pore indices 25..29")

    m = brine.salt_mass(C_RES, layer.V, T=T_DUMUX)
    initial_salt = float(m.sum())
    initial_liquid = float(np.sum(brine.density(C_RES, T=T_DUMUX) * layer.V))
    initial_water = initial_liquid - initial_salt
    precip = (PorePrecipitation(layer, brine.C_SAT, RHO_CR,
                                max_solid_fraction=max_solid_fraction, T=T_DUMUX)
              if precipitation else None)
    g_lat, g_res = layer.conductances()
    d_lat, d_res = layer.diffusive()
    n_output = int(round(output_interval / dt))

    records = {name: [] for name in (
        "time_s", "pore_x", "pore_w", "pore_pressure_gauge_pa", "pore_solid_fraction",
        "pore_volume_m3", "top_evaporation_kg_s", "throat_flow_m3_s",
        "throat_radius_m", "throat_pressure_gauge_pa", "salt_reservoir_in_kg",
        "salt_balance_residual_kg", "liquid_reservoir_in_kg",
        "water_reservoir_in_kg", "water_balance_residual_kg",
        "cumulative_evaporation_kg"
    )}
    salt_in = 0.0
    liquid_in = 0.0
    water_in = 0.0
    water_out = 0.0
    first_saturation_time = None
    stop_reason = "t_end"
    last_flow = (np.zeros(layer.N), np.zeros(len(layer.lat)), np.zeros(layer.N))
    last_top_flux = np.zeros(layer.n_surface)

    def save(t):
        w = brine.mass_fraction_from_salt_mass(m, layer.V, T=T_DUMUX)
        x = brine.mole_fraction(w)
        p, Q_lat, Q_res = last_flow
        p_full = _full_pore_state(layer, p, 0.0)
        x_full = _full_pore_state(layer, x, X_RES)
        w_full = _full_pore_state(layer, w, C_RES)
        phi = np.zeros(layer.N) if precip is None else precip.solid_volume_fraction
        phi_full = _full_pore_state(layer, phi, 0.0)
        V_full = np.asarray(layer.dgf_vertices[:, 3], float).copy()
        V_full[layer.file_vertex_index] = layer.V
        Q_full, r_full, p_throat = _full_throat_state(layer, Q_lat, Q_res, p_full)
        solid = 0.0 if precip is None else float(precip.m_solid.sum())
        residual = float(m.sum() + solid - initial_salt - salt_in)
        liquid_inventory = float(np.sum(brine.density(w, T=T_DUMUX) * layer.V))
        water_inventory = liquid_inventory - float(m.sum())
        water_residual = water_inventory - initial_water - water_in + water_out
        values = {
            "time_s": float(t), "pore_x": x_full, "pore_w": w_full,
            "pore_pressure_gauge_pa": p_full, "pore_solid_fraction": phi_full,
            "pore_volume_m3": V_full, "top_evaporation_kg_s": last_top_flux.copy(),
            "throat_flow_m3_s": Q_full, "throat_radius_m": r_full,
            "throat_pressure_gauge_pa": p_throat, "salt_reservoir_in_kg": salt_in,
            "salt_balance_residual_kg": residual, "liquid_reservoir_in_kg": liquid_in,
            "water_reservoir_in_kg": water_in, "water_balance_residual_kg": water_residual,
            "cumulative_evaporation_kg": water_out,
        }
        for key, value in values.items():
            records[key].append(value)

    save(0.0)
    n_steps = int(np.ceil(t_end / dt))
    for step_idx in range(1, n_steps + 1):
        t0 = (step_idx - 1) * dt
        step_dt = min(dt, t_end - t0)
        if step_dt <= 0.0:
            break
        w_before = brine.mass_fraction_from_salt_mass(m, layer.V, T=T_DUMUX)
        x_before = brine.mole_fraction(w_before)
        rho_before = brine.density(w_before, T=T_DUMUX)
        if replay_flux is None:
            top_flux = np.asarray(boundary_layer_evaporation(
                w_before[layer.surface], layer.A_open[layer.surface], DELTA_BL, RH_ENV), float)
        else:
            top_flux = _replay_at(replay_flux, t0)
        if np.any(top_flux < 0.0):
            raise ValueError("this evaporation comparison does not admit a negative replay flux")
        q_mass = np.zeros(layer.N)
        q_mass[layer.surface] = top_flux
        if precip is not None:
            q_mass -= precip.last_volume_rate * rho_before

        out = network_step(layer, m, C_RES, q_mass, g_lat, g_res, d_lat, d_res,
                           step_dt, T=T_DUMUX)
        m_raw = out.m
        w_raw = brine.mass_fraction_from_salt_mass(m_raw, layer.V, T=T_DUMUX)
        x_raw = brine.mole_fraction(w_raw)
        if first_saturation_time is None and float(x_raw.max()) >= X_SOLUBILITY:
            x0_max = float(x_before.max())
            x1_max = float(x_raw.max())
            fraction = ((X_SOLUBILITY - x0_max) / (x1_max - x0_max)
                        if x1_max > x0_max else 1.0)
            first_saturation_time = t0 + step_dt * float(np.clip(fraction, 0.0, 1.0))

        m = m_raw
        salt_in += out.salt_in
        rho_res = float(brine.density(C_RES, T=T_DUMUX))
        q_res = np.asarray(out.Q_res, float)
        liquid_in_step = step_dt * float(np.sum(
            np.maximum(q_res, 0.0) * rho_res
            - np.maximum(-q_res, 0.0) * rho_before))
        liquid_in += liquid_in_step
        # The salt ledger includes both advection and Fickian diffusion.  In the mass-average
        # component balance the latter has an equal-and-opposite solvent contribution, hence
        # reservoir water is total liquid advection minus the *net* salt transfer, not merely
        # minus the advective part of the salt transfer.
        water_in += liquid_in_step - out.salt_in
        water_out += float(top_flux.sum()) * step_dt
        last_flow = (out.p, out.Q_lat, out.Q_res)
        last_top_flux = top_flux
        if precip is not None:
            try:
                m = precip.step(m, step_dt)
            except PoreClogged as error:
                stop_reason = str(error)
                save(t0 + step_dt)
                break
            if np.any(precip.last_volume_change > 0.0):
                g_lat, g_res = layer.conductances()
                d_lat, d_res = layer.diffusive()
        if step_idx % n_output == 0 or t0 + step_dt >= t_end:
            save(t0 + step_dt)

    arrays = {key: np.asarray(value) for key, value in records.items()}
    arrays.update({
        "dgf_vertices": layer.dgf_vertices,
        "dgf_edges": layer.dgf_edges,
        "file_vertex_index": layer.file_vertex_index,
        "reservoir_file_indices": layer.reservoir_file_indices,
        "lat_file_edge_index": layer.lat_file_edge_index,
        "res_file_edge_index": layer.res_file_edge_index,
    })
    final_solid = 0.0 if precip is None else float(precip.m_solid.sum())
    final_residual = float(arrays["salt_balance_residual_kg"][-1])
    final_water_residual = float(arrays["water_balance_residual_kg"][-1])
    scale = max(abs(initial_salt + salt_in), 1e-300)
    water_scale = max(abs(initial_water + water_in), 1e-300)
    summary = {
        "dt_s": dt, "requested_t_end_s": t_end, "actual_t_end_s": float(arrays["time_s"][-1]),
        "output_interval_s": output_interval, "x_reservoir": X_RES,
        "w_reservoir": C_RES, "temperature_K": T_DUMUX,
        "boundary_layer_m": DELTA_BL, "relative_humidity": RH_ENV,
        "first_saturation_time_s": first_saturation_time,
        "final_solid_mass_kg": final_solid, "initial_dissolved_salt_kg": initial_salt,
        "reservoir_salt_in_kg": salt_in, "salt_balance_residual_kg": final_residual,
        "salt_balance_relative_residual": final_residual / scale,
        "reservoir_liquid_in_kg": liquid_in, "reservoir_water_in_kg": water_in,
        "water_balance_residual_kg": final_water_residual,
        "water_balance_relative_residual": final_water_residual / water_scale,
        "cumulative_evaporation_kg": water_out, "stop_reason": stop_reason,
        "boundary_mode": "dumux_flux_replay" if replay_flux is not None else "python_self_driven",
        "precipitation": bool(precipitation),
    }
    return arrays, summary


def write_outputs(output_dir, arrays, summary):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "python_results.npz", **arrays)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    with (output_dir / "top_evaporation.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", *(f"pore_{idx}" for idx in TOP_FILE_INDICES), "total_kg_s"])
        for t, values in zip(arrays["time_s"], arrays["top_evaporation_kg_s"]):
            writer.writerow([f"{t:.17g}", *(f"{v:.17g}" for v in values), f"{values.sum():.17g}"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dt", type=float, default=5.0)
    parser.add_argument("--t-end", type=float, default=40000.0)
    parser.add_argument("--output-interval", type=float, default=100.0)
    parser.add_argument("--dgf", default=DGF_DUMUX_5X5_VARTHROAT)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--flux-replay")
    parser.add_argument("--no-precipitation", action="store_true")
    args = parser.parse_args(argv)
    replay = load_flux_replay(args.flux_replay) if args.flux_replay else None
    arrays, summary = run(dt=args.dt, t_end=args.t_end, output_interval=args.output_interval,
                          dgf=args.dgf, precipitation=not args.no_precipitation,
                          replay_flux=replay)
    write_outputs(args.output_dir, arrays, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

