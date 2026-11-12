"""Normalize DuMuX VTP output, compare it with Python, and draw pore/throat pressure maps."""

from __future__ import annotations

import argparse
import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Circle

M_H2O = 0.01801518
M_NACL = 0.05844
RHO_CR = 2170.0
X_RES = 0.08
X_SOL = 0.09956503729854478
RESERVOIR = np.array([0, 5, 10, 15, 20], dtype=int)
TOP = np.array([25, 26, 27, 28, 29], dtype=int)


def _array_map(vtp_path):
    piece = ET.parse(vtp_path).getroot().find(".//Piece")
    arrays = {}
    for node in piece.findall(".//DataArray"):
        name = node.attrib.get("Name")
        if name:
            arrays[name] = np.fromstring(node.text or "", sep=" ")
    return piece, arrays


def _pvd_entries(path):
    root = ET.parse(path).getroot()
    return [(float(node.attrib["timestep"]), node.attrib["file"])
            for node in root.findall(".//DataSet")]


def reconstruct_gauge_pressure(coordinates, edges, volume_flux, transmissibility):
    """Recover sub-Float32 pressure differences from DuMuX Q/T diagnostics.

    VTP stores absolute liquid pressure as Float32, which rounds the O(1e-2 Pa) network drop out
    of a 1e5 Pa baseline.  ``volumeFluxW/transmissibilityW`` retains each edge's pressure-drop
    magnitude.  Vertical flow is reservoir-to-surface; horizontal orientations are inferred from
    the vertical-column potentials and iterated through a weighted graph least-squares solve.
    """
    coordinates = np.asarray(coordinates, float)
    edges = np.asarray(edges, int)
    q_abs = np.abs(np.asarray(volume_flux, float))
    trans = np.asarray(transmissibility, float)
    dp_abs = np.divide(q_abs, trans, out=np.zeros_like(q_abs), where=trans > 0.0)
    bottom = np.flatnonzero(coordinates[:, 1] == coordinates[:, 1].min())
    is_bottom = np.zeros(len(coordinates), bool); is_bottom[bottom] = True
    dynamic = np.flatnonzero(~is_bottom)
    to_unknown = np.full(len(coordinates), -1, int); to_unknown[dynamic] = np.arange(dynamic.size)

    p = np.zeros(len(coordinates))
    for x in np.unique(coordinates[:, 0]):
        column = np.flatnonzero(coordinates[:, 0] == x)
        column = column[np.argsort(coordinates[column, 1])]
        for low, high in zip(column[:-1], column[1:]):
            hit = np.flatnonzero(((edges[:, 0] == low) & (edges[:, 1] == high)) |
                                 ((edges[:, 0] == high) & (edges[:, 1] == low)))
            if hit.size == 1:
                p[high] = p[low] - dp_abs[hit[0]]

    for _ in range(5):
        rows, rhs, weights = [], [], []
        for edge_idx, (a, b) in enumerate(edges):
            if is_bottom[a] and is_bottom[b]:
                continue
            delta_y = coordinates[b, 1] - coordinates[a, 1]
            if delta_y != 0.0:
                sign = 1.0 if delta_y > 0.0 else -1.0
            else:
                sign = 1.0 if p[a] >= p[b] else -1.0
            row = np.zeros(dynamic.size)
            if not is_bottom[a]: row[to_unknown[a]] += 1.0
            if not is_bottom[b]: row[to_unknown[b]] -= 1.0
            rows.append(row)
            rhs.append(sign * dp_abs[edge_idx])
            weights.append(max(trans[edge_idx], np.max(trans) * 1e-12))
        A = np.asarray(rows); b = np.asarray(rhs); w = np.sqrt(np.asarray(weights) / np.max(weights))
        sol = np.linalg.lstsq(A * w[:, None], b * w, rcond=None)[0]
        p[dynamic] = sol; p[bottom] = 0.0
    predicted = trans * (p[edges[:, 0]] - p[edges[:, 1]])
    relative_magnitude_residual = float(
        np.linalg.norm(np.abs(predicted) - q_abs) / max(np.linalg.norm(q_abs), 1e-300)
    )
    return p, predicted, relative_magnitude_residual


def read_evaporation(path):
    rows = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            values = [float(value) for value in line.split(",")]
            if len(values) < 24:
                continue
            times = values[0]
            per_pore_molar = np.asarray(values[15::2][:5], float)
            rows.append((times, *(-per_pore_molar * M_H2O)))
    data = np.asarray(rows)
    order = np.argsort(data[:, 0]); data = data[order]
    keep = np.r_[True, np.diff(data[:, 0]) > 0.0]
    return data[keep, 0], data[keep, 1:]


def normalize_dumux(raw_dir, output_npz):
    raw_dir = Path(raw_dir)
    pvd = next(raw_dir.glob("*.pvd"))
    evap_t, evap = read_evaporation(raw_dir / "evaporation.csv")
    records = {key: [] for key in (
        "time_s", "pore_x", "pore_liquid_saturation", "pore_pressure_gauge_pa", "pore_solid_fraction",
        "pore_volume_m3", "throat_flow_m3_s", "throat_radius_m",
        "throat_pressure_gauge_pa", "pressure_reference_pa",
        "pressure_reconstruction_residual"
    )}
    coords0 = edges0 = initial_volume = None
    for time_s, filename in _pvd_entries(pvd):
        _, arrays = _array_map(raw_dir / filename)
        coords = arrays["Coordinates"].reshape(-1, 3)[:, :2]
        edges = arrays["connectivity"].astype(int).reshape(-1, 2)
        pressure, flow, residual = reconstruct_gauge_pressure(
            coords, edges, arrays["volumeFluxW"], arrays["transmissibilityW"])
        bottom = np.flatnonzero(coords[:, 1] == coords[:, 1].min())
        pressure_reference = float(np.mean(arrays["p_liq"][bottom]))
        if coords0 is None:
            coords0, edges0 = coords, edges
            initial_volume = arrays["poreVolume"].copy()
        records["time_s"].append(time_s)
        records["pore_x"].append(arrays["x^NaCl_liq"])
        records["pore_liquid_saturation"].append(arrays["S_liq"])
        records["pore_pressure_gauge_pa"].append(pressure)
        records["pore_solid_fraction"].append(arrays["precipitateVolumeFraction^NaCl"])
        records["pore_volume_m3"].append(arrays["poreVolume"])
        records["throat_flow_m3_s"].append(flow)
        records["throat_radius_m"].append(arrays["throatInscribedRadius"])
        records["throat_pressure_gauge_pa"].append(
            0.5 * (pressure[edges[:, 0]] + pressure[edges[:, 1]]))
        records["pressure_reference_pa"].append(pressure_reference)
        records["pressure_reconstruction_residual"].append(residual)
    result = {key: np.asarray(value) for key, value in records.items()}
    result["dgf_vertices"] = np.column_stack([coords0, np.zeros((len(coords0), 3))])
    result["dgf_edges"] = np.column_stack([edges0, np.zeros((len(edges0), 3))])
    result["initial_pore_volume_m3"] = initial_volume
    result["top_evaporation_kg_s"] = np.column_stack([
        np.interp(result["time_s"], evap_t, evap[:, j]) for j in range(5)])
    total = evap.sum(axis=1)
    cumulative = np.zeros_like(evap_t)
    cumulative[1:] = np.cumsum(0.5 * (total[1:] + total[:-1]) * np.diff(evap_t))
    result["cumulative_evaporation_kg"] = np.interp(result["time_s"], evap_t, cumulative)
    solid_mass = (result["pore_solid_fraction"] * initial_volume[None, :] * RHO_CR).sum(axis=1)
    result["solid_mass_kg"] = solid_mass
    np.savez_compressed(output_npz, **result)
    replay_path = Path(output_npz).with_name(Path(output_npz).stem + "_flux_replay.csv")
    with replay_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", *(f"pore_{idx}" for idx in TOP)])
        for t, values in zip(evap_t, evap):
            writer.writerow([f"{t:.17g}", *(f"{v:.17g}" for v in values)])
    return result, replay_path


def load_double_diagnostics(pore_csv, throat_csv, geometry):
    """Load the validation-only Float64 state and signed flux CSVs emitted by DuMuX."""
    pore = np.genfromtxt(pore_csv, delimiter=",", names=True)
    throat = np.genfromtxt(throat_csv, delimiter=",", names=True)
    times = np.unique(pore["time_s"])
    result = {key: [] for key in (
        "pore_x", "pore_liquid_saturation", "pore_pressure_gauge_pa",
        "pore_pressure_absolute_pa", "pore_solid_fraction", "pore_volume_m3",
        "pressure_reference_pa", "throat_flow_m3_s", "throat_radius_m",
        "throat_pressure_gauge_pa", "pressure_reconstruction_residual"
    )}
    edges = np.asarray(geometry["dgf_edges"][:, :2], int)
    for time_s in times:
        ps = np.sort(pore[pore["time_s"] == time_s], order="pore_index")
        ts = np.sort(throat[throat["time_s"] == time_s], order="throat_index")
        if len(ps) != 30 or len(ts) != 45:
            raise ValueError(f"incomplete double diagnostic at {time_s:g} s")
        absolute = np.asarray(ps["p_liq_pa"], float)
        reference = float(np.mean(absolute[RESERVOIR]))
        gauge = absolute - reference
        csv_edges = np.column_stack([ts["pore_a"], ts["pore_b"]]).astype(int)
        if not np.array_equal(csv_edges, edges):
            raise ValueError(f"throat ordering mismatch at {time_s:g} s")
        result["pore_x"].append(ps["x_nacl_liq"])
        result["pore_liquid_saturation"].append(ps["s_liq"])
        result["pore_pressure_gauge_pa"].append(gauge)
        result["pore_pressure_absolute_pa"].append(absolute)
        result["pore_solid_fraction"].append(ps["solid_volume_fraction"])
        result["pore_volume_m3"].append(ps["pore_volume_m3"])
        result["pressure_reference_pa"].append(reference)
        result["throat_flow_m3_s"].append(ts["q_liq_a_to_b_m3_s"])
        result["throat_radius_m"].append(ts["throat_radius_m"])
        result["throat_pressure_gauge_pa"].append(
            0.5 * (gauge[edges[:, 0]] + gauge[edges[:, 1]]))
        result["pressure_reconstruction_residual"].append(0.0)
    result = {key: np.asarray(value) for key, value in result.items()}
    result["time_s"] = times
    for key in ("dgf_vertices", "dgf_edges", "initial_pore_volume_m3"):
        if key in geometry:
            result[key] = geometry[key]
    return result


def reconstruct_dumux_ledger(pore_csv, throat_csv, evaporation_csv, output_csv):
    """Reconstruct the dynamic-pore water/NaCl mole books from 1 s Float64 diagnostics."""
    pore = np.genfromtxt(pore_csv, delimiter=",", names=True)
    throat = np.genfromtxt(throat_csv, delimiter=",", names=True)
    times = np.unique(pore["time_s"])
    reservoir = set(RESERVOIR.tolist())
    dynamic = np.array([idx for idx in range(30) if idx not in reservoir], int)
    water_inventory, salt_inventory = [], []
    reservoir_water_flux, reservoir_salt_flux = [], []
    for time_s in times:
        ps = np.sort(pore[pore["time_s"] == time_s], order="pore_index")
        ts = np.sort(throat[throat["time_s"] == time_s], order="throat_index")
        water_inventory.append(float(np.sum(ps["water_inventory_mol"][dynamic])))
        salt_inventory.append(float(np.sum(
            ps["dissolved_nacl_inventory_mol"][dynamic]
            + ps["solid_nacl_inventory_mol"][dynamic])))
        water_flux = salt_flux = 0.0
        for row in ts:
            a, b = int(row["pore_a"]), int(row["pore_b"])
            orientation = (1.0 if a in reservoir and b not in reservoir
                           else -1.0 if b in reservoir and a not in reservoir else 0.0)
            water_flux += orientation * float(row["total_water_all_phases_mol_s"])
            salt_flux += orientation * float(row["total_nacl_all_phases_mol_s"])
        reservoir_water_flux.append(water_flux)
        reservoir_salt_flux.append(salt_flux)
    water_inventory = np.asarray(water_inventory)
    salt_inventory = np.asarray(salt_inventory)
    reservoir_water_flux = np.asarray(reservoir_water_flux)
    reservoir_salt_flux = np.asarray(reservoir_salt_flux)
    dt = np.diff(times)
    # The 1 s diagnostic run is forced to land on every checkpoint.  DuMuX uses backward Euler,
    # so the end-of-step (right rectangle) flux is the discrete flux book to compare with storage.
    reservoir_water = np.r_[0.0, np.cumsum(reservoir_water_flux[1:] * dt)]
    reservoir_salt = np.r_[0.0, np.cumsum(reservoir_salt_flux[1:] * dt)]
    evap_t, evap_kg_s = read_evaporation(evaporation_csv)
    evap_mol_s = evap_kg_s.sum(axis=1) / M_H2O
    cumulative_evap = np.r_[0.0, np.cumsum(evap_mol_s[1:] * np.diff(evap_t))]
    cumulative_evap = np.interp(times, evap_t, cumulative_evap)
    water_residual = water_inventory - water_inventory[0] - reservoir_water + cumulative_evap
    salt_residual = salt_inventory - salt_inventory[0] - reservoir_salt
    water_scale = np.maximum(np.abs(water_inventory[0] + reservoir_water), 1e-300)
    salt_scale = np.maximum(np.abs(salt_inventory[0] + reservoir_salt), 1e-300)
    water_relative = water_residual / water_scale
    salt_relative = salt_residual / salt_scale
    with Path(output_csv).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "water_inventory_mol", "salt_inventory_mol",
                         "reservoir_water_flux_mol_s", "reservoir_salt_flux_mol_s",
                         "cumulative_reservoir_water_mol", "cumulative_reservoir_salt_mol",
                         "cumulative_evaporation_water_mol", "water_residual_mol",
                         "salt_residual_mol", "water_relative_residual",
                         "salt_relative_residual"])
        writer.writerows(zip(times, water_inventory, salt_inventory,
                             reservoir_water_flux, reservoir_salt_flux,
                             reservoir_water, reservoir_salt, cumulative_evap,
                             water_residual, salt_residual, water_relative, salt_relative))
    return {
        "dumux_ledger_end_time_s": float(times[-1]),
        "dumux_water_ledger_max_abs_relative_residual_0_1000s": float(
            np.max(np.abs(water_relative))),
        "dumux_salt_ledger_max_abs_relative_residual_0_1000s": float(
            np.max(np.abs(salt_relative))),
        "dumux_water_ledger_relative_residual_940s": float(
            water_relative[_at({"time_s": times}, 940.0)]),
        "dumux_salt_ledger_relative_residual_940s": float(
            salt_relative[_at({"time_s": times}, 940.0)]),
    }


def load_npz(path):
    with np.load(path) as data:
        return {key: data[key] for key in data.files}


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _at(data, time_s):
    return int(np.argmin(np.abs(data["time_s"] - time_s)))


def _interp(data, key, time_s):
    """Linearly interpolate one scalar/vector checkpoint field at a physical time."""
    times = np.asarray(data["time_s"], float)
    values = np.asarray(data[key], float)
    if values.ndim == 1:
        return float(np.interp(time_s, times, values))
    shape = values.shape[1:]
    flat = values.reshape(len(times), -1)
    out = np.array([np.interp(time_s, times, flat[:, j]) for j in range(flat.shape[1])])
    return out.reshape(shape)


def merge_time_series(primary, supplement):
    """Add high-frequency diagnostic checkpoints without changing static geometry fields."""
    result = dict(primary)
    for key in primary.keys() & supplement.keys():
        a, b = np.asarray(primary[key]), np.asarray(supplement[key])
        if a.ndim and b.ndim and a.shape[0] == len(primary["time_s"]) \
                and b.shape[0] == len(supplement["time_s"]) and a.shape[1:] == b.shape[1:]:
            times = np.concatenate([primary["time_s"], supplement["time_s"]])
            values = np.concatenate([a, b], axis=0)
            order = np.argsort(times, kind="stable")
            times, values = times[order], values[order]
            keep = np.r_[np.diff(times) != 0.0, True]  # supplement wins at duplicate times
            result[key] = values[keep]
            if key == "time_s":
                result[key] = times[keep]
    return result


def first_solid_time(data):
    solid = np.max(data["pore_solid_fraction"], axis=1)
    hit = np.flatnonzero(solid > 0.0)
    return None if hit.size == 0 else float(data["time_s"][hit[0]])


def first_solubility_crossing(data):
    maximum = np.max(data["pore_x"], axis=1)
    hit = np.flatnonzero(maximum >= X_SOL)
    if hit.size == 0:
        return None
    i = int(hit[0])
    if i == 0 or maximum[i] <= maximum[i - 1]:
        return float(data["time_s"][i])
    fraction = (X_SOL - maximum[i - 1]) / (maximum[i] - maximum[i - 1])
    return float(data["time_s"][i - 1] + fraction
                 * (data["time_s"][i] - data["time_s"][i - 1]))


def _network_panel(ax, coordinates, edges, pore_radius, throat_radius, pore_value, throat_value,
                   norm, cmap, title):
    xy = coordinates * 1e3
    segments = np.stack([xy[edges[:, 0]], xy[edges[:, 1]]], axis=1)
    widths = 0.7 + 4.0 * throat_radius / max(float(np.max(throat_radius)), 1e-300)
    lines = LineCollection(segments, cmap=cmap, norm=norm, linewidths=widths, zorder=1)
    lines.set_array(throat_value); ax.add_collection(lines)
    patches = [Circle((x, y), r * 1e3) for (x, y), r in zip(xy, pore_radius)]
    pores = PatchCollection(patches, cmap=cmap, norm=norm, edgecolor="black", linewidth=0.5, zorder=2)
    pores.set_array(pore_value); ax.add_collection(pores)
    for pore in RESERVOIR:
        ax.add_patch(Circle(tuple(xy[pore]), pore_radius[pore] * 1e3, facecolor="none",
                            edgecolor="#2166ac", linewidth=1.5, zorder=3))
    for pore in TOP:
        ax.add_patch(Circle(tuple(xy[pore]), pore_radius[pore] * 1e3, facecolor="none",
                            edgecolor="#b2182b", linewidth=1.5, zorder=3))
    ax.set_aspect("equal"); ax.set_xlim(-0.45, 3.25); ax.set_ylim(-0.45, 3.95)
    ax.set_title(title, fontsize=9); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    return pores


def pressure_figure(dumux, python, times, output_base):
    coordinates = dumux["dgf_vertices"][:, :2]
    edges = dumux["dgf_edges"][:, :2].astype(int)
    pore_radius = python["dgf_vertices"][:, 2]
    fig, axes = plt.subplots(len(times), 3, figsize=(10.8, 3.7 * len(times)), constrained_layout=True)
    for row, time_s in enumerate(times):
        dp = _interp(dumux, "pore_pressure_gauge_pa", time_s)
        pp = _interp(python, "pore_pressure_gauge_pa", time_s)
        dtp = _interp(dumux, "throat_pressure_gauge_pa", time_s)
        ptp = _interp(python, "throat_pressure_gauge_pa", time_s)
        low, high = min(float(dp.min()), float(pp.min())), max(float(dp.max()), float(pp.max()))
        if high <= low: high = low + 1e-12
        norm = Normalize(low, high)
        m0 = _network_panel(axes[row, 0], coordinates, edges, pore_radius,
                            _interp(dumux, "throat_radius_m", time_s), dp, dtp, norm, "viridis",
                            f"DuMuX, t={time_s:g} s")
        _network_panel(axes[row, 1], coordinates, edges, pore_radius,
                       _interp(python, "throat_radius_m", time_s), pp, ptp, norm, "viridis",
                       f"Python, t={time_s:g} s")
        diff_p, diff_t = pp - dp, ptp - dtp
        bound = max(float(np.max(np.abs(diff_p))), float(np.max(np.abs(diff_t))), 1e-12)
        diff_norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
        m2 = _network_panel(axes[row, 2], coordinates, edges, pore_radius,
                            _interp(python, "throat_radius_m", time_s), diff_p, diff_t, diff_norm, "coolwarm",
                            f"Python - DuMuX, t={time_s:g} s")
        fig.colorbar(m0, ax=axes[row, :2], label="gauge pressure (Pa)", shrink=0.82)
        fig.colorbar(m2, ax=axes[row, 2], label="pressure difference (Pa)", shrink=0.82)
    fig.suptitle("5 x 5 pore-network pressure distribution (pores and throats)\n"
                 r"colors: $\Delta p=p-p_{res}$; model-specific $p_{res}$ and absolute pressure are in the CSV",
                 fontsize=12)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(f"{output_base}.{suffix}", dpi=600 if suffix == "png" else None)
    plt.close(fig)


def pressure_range_figure(dumux, python, output_base, dumux_direct=None):
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 6.5), sharex=True, constrained_layout=True)
    series = [(python, "Python", "-"), (dumux, "DuMuX Q/T reconstructed", "--")]
    if dumux_direct is not None:
        series.insert(0, (dumux_direct, "DuMuX direct Float64", "-"))
    for data, label, linestyle in series:
        span = np.ptp(data["pore_pressure_gauge_pa"], axis=1)
        axes[0].plot(data["time_s"], span, linestyle=linestyle, label=label)
        positive = span > 0.0
        axes[1].semilogy(data["time_s"][positive], span[positive],
                         linestyle=linestyle, label=label)
    axes[0].set_ylabel("pressure range (Pa)")
    axes[0].set_title("Full linear scale")
    axes[1].set_xlabel("time (s)"); axes[1].set_ylabel("pressure range (Pa)")
    axes[1].set_title("Log scale (resolves the pre-failure and Python ranges)")
    for ax in axes:
        ax.grid(alpha=0.25); ax.legend()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(f"{output_base}.{suffix}", dpi=600 if suffix == "png" else None)
    plt.close(fig)


def pressure_single_figure(data, time_s, output_base, title):
    coordinates = data["dgf_vertices"][:, :2]
    edges = data["dgf_edges"][:, :2].astype(int)
    pore_radius = data["dgf_vertices"][:, 2]
    pore_pressure = _interp(data, "pore_pressure_gauge_pa", time_s)
    throat_pressure = _interp(data, "throat_pressure_gauge_pa", time_s)
    throat_radius = _interp(data, "throat_radius_m", time_s)
    low, high = float(pore_pressure.min()), float(pore_pressure.max())
    if high <= low:
        high = low + 1e-12
    fig, ax = plt.subplots(figsize=(5.2, 5.0), constrained_layout=True)
    mappable = _network_panel(ax, coordinates, edges, pore_radius, throat_radius,
                              pore_pressure, throat_pressure, Normalize(low, high), "viridis", title)
    fig.colorbar(mappable, ax=ax, label="gauge pressure (Pa)")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(f"{output_base}.{suffix}", dpi=600 if suffix == "png" else None)
    plt.close(fig)


def write_snapshot_csv(dumux, python, times, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    coordinates = dumux["dgf_vertices"][:, :2]
    edges = dumux["dgf_edges"][:, :2].astype(int)
    pore_radius = python["dgf_vertices"][:, 2]
    for time_s in times:
        for data, model in ((dumux, "dumux"), (python, "python")):
            pore_pressure = _interp(data, "pore_pressure_gauge_pa", time_s)
            pore_x = _interp(data, "pore_x", time_s)
            pore_solid = _interp(data, "pore_solid_fraction", time_s)
            pore_volume = _interp(data, "pore_volume_m3", time_s)
            throat_radius = _interp(data, "throat_radius_m", time_s)
            throat_pressure = _interp(data, "throat_pressure_gauge_pa", time_s)
            throat_flow = _interp(data, "throat_flow_m3_s", time_s)
            pressure_reference = (_interp(data, "pressure_reference_pa", time_s)
                                  if "pressure_reference_pa" in data else 1.0e5)
            with (out_dir / f"pores_{model}_{time_s:g}s.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["file_pore", "x_m", "y_m", "radius_m", "pressure_reference_pa",
                                 "pressure_gauge_pa", "pressure_absolute_pa",
                                 "x_nacl", "solid_volume_fraction", "pore_volume_m3"])
                for pore in range(len(coordinates)):
                    writer.writerow([pore, *coordinates[pore], pore_radius[pore],
                                     pressure_reference, pore_pressure[pore],
                                     pressure_reference + pore_pressure[pore],
                                     pore_x[pore], pore_solid[pore],
                                     pore_volume[pore]])
            with (out_dir / f"throats_{model}_{time_s:g}s.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["file_throat", "pore_a", "pore_b", "radius_m",
                                 "pressure_reference_pa", "pressure_mean_gauge_pa",
                                 "pressure_mean_absolute_pa", "pressure_a_gauge_pa",
                                 "pressure_b_gauge_pa", "pressure_a_absolute_pa",
                                 "pressure_b_absolute_pa", "flow_a_to_b_m3_s"])
                for edge_idx, (a, b) in enumerate(edges):
                    writer.writerow([edge_idx, a, b, throat_radius[edge_idx],
                                     pressure_reference, throat_pressure[edge_idx],
                                     pressure_reference + throat_pressure[edge_idx],
                                     pore_pressure[a], pore_pressure[b],
                                     pressure_reference + pore_pressure[a],
                                     pressure_reference + pore_pressure[b],
                                     throat_flow[edge_idx]])


def _relative(a, b):
    return float(abs(float(a) - float(b)) / max(abs(float(b)), 1e-300))


def _solid_mass_python(data, index):
    return float((data["pore_solid_fraction"][index]
                  * data["dgf_vertices"][:, 3] * RHO_CR).sum())


def _edge_alignment(source_edges, target_edges):
    """Map target DGF edges to source edges and retain orientation for signed flow."""
    lookup = {}
    for idx, (a, b) in enumerate(np.asarray(source_edges)[:, :2].astype(int)):
        lookup[tuple(sorted((int(a), int(b))))] = (idx, int(a), int(b))
    indices, signs = [], []
    for a, b in np.asarray(target_edges)[:, :2].astype(int):
        idx, sa, sb = lookup[tuple(sorted((int(a), int(b))))]
        indices.append(idx)
        signs.append(1.0 if (sa, sb) == (int(a), int(b)) else -1.0)
    return np.asarray(indices, int), np.asarray(signs, float)


def comparison_metrics(d10, d5, p10, p5, p2p5=None, replay=None, d_pre=None,
                       dumux_first_saturation_time=None,
                       python_first_saturation_time=None,
                       python_summaries=None):
    common_end = min(d5["time_s"][-1], p5["time_s"][-1])
    final_d5, final_p5 = _at(d5, common_end), _at(p5, common_end)
    # DuMuX first reaches saturation at about 936 s, so 900 s is the latest common saved
    # checkpoint that is unambiguously single-phase in both models.
    pre = min(900.0, common_end)
    d_pre = d5 if d_pre is None else d_pre
    pre_d, pre_d_coarse, pre_p = _at(d_pre, pre), _at(d5, pre), _at(p5, pre)
    metrics = {
        "common_end_s": float(common_end),
        "dumux_first_saved_solid_time_s": first_solid_time(d5),
        "python_first_solid_time_s": first_solid_time(p5),
        "pre_saturation_comparison_time_s": pre,
        "self_driven_pre_saturation_x_max_abs_difference": float(np.max(np.abs(
            d_pre["pore_x"][pre_d] - p5["pore_x"][pre_p]))),
        "final_pressure_weighted_relative_error": float(np.linalg.norm(
            p5["pore_pressure_gauge_pa"][final_p5] - d5["pore_pressure_gauge_pa"][final_d5]) /
            max(np.linalg.norm(d5["pore_pressure_gauge_pa"][final_d5]), 1e-300)),
        "final_solid_mass_dumux_kg": float(d5["solid_mass_kg"][final_d5]),
        "final_solid_mass_python_kg": _solid_mass_python(p5, final_p5),
        "final_min_throat_radius_dumux_m": float(d5["throat_radius_m"][final_d5].min()),
        "final_min_throat_radius_python_m": float(p5["throat_radius_m"][final_p5].min()),
        "dumux_pressure_reconstruction_max_relative_residual": float(
            np.max(d5["pressure_reconstruction_residual"])),
        "dumux_pressure_reconstruction_pre_saturation_relative_residual": float(
            d_pre["pressure_reconstruction_residual"][pre_d]),
        "dumux_min_liquid_saturation_100s": float(
            np.min(d5["pore_liquid_saturation"][_at(d5, 100.0)])),
        "dumux_min_liquid_saturation_pre_comparison": float(
            np.min(d_pre["pore_liquid_saturation"][pre_d])),
        "dumux_min_liquid_saturation_last_saved": float(
            np.min(d5["pore_liquid_saturation"][final_d5])),
    }
    metrics["self_driven_pre_saturation_cumulative_evap_relative_difference"] = _relative(
        p5["cumulative_evaporation_kg"][pre_p],
        d5["cumulative_evaporation_kg"][pre_d_coarse])
    metrics["final_solid_mass_relative_difference"] = _relative(
        metrics["final_solid_mass_python_kg"], metrics["final_solid_mass_dumux_kg"])
    metrics["final_pore_solid_fraction_max_abs_difference"] = float(np.max(np.abs(
        p5["pore_solid_fraction"][final_p5] - d5["pore_solid_fraction"][final_d5])))
    metrics["final_min_throat_radius_relative_difference"] = _relative(
        metrics["final_min_throat_radius_python_m"], metrics["final_min_throat_radius_dumux_m"])

    edge_idx, edge_sign = _edge_alignment(d_pre["dgf_edges"], p5["dgf_edges"])
    dq = d_pre["throat_flow_m3_s"][pre_d, edge_idx] * edge_sign
    pq = p5["throat_flow_m3_s"][pre_p]
    effective = np.abs(dq) > max(float(np.max(np.abs(dq))) * 1e-6, 1e-30)
    metrics["pre_saturation_effective_throat_count"] = int(np.count_nonzero(effective))
    metrics["pre_saturation_effective_flow_direction_agreement_fraction"] = float(
        np.mean(np.sign(dq[effective]) == np.sign(pq[effective]))) if np.any(effective) else 1.0
    metrics["pre_saturation_flow_weighted_relative_error"] = float(
        np.sum(np.abs(pq[effective] - dq[effective]))
        / max(np.sum(np.abs(dq[effective])), 1e-300))
    metrics["pre_saturation_near_zero_flow_max_abs_error_m3_s"] = float(
        np.max(np.abs(pq[~effective] - dq[~effective]))) if np.any(~effective) else 0.0
    metrics["pre_saturation_pressure_weighted_relative_error"] = float(np.linalg.norm(
        p5["pore_pressure_gauge_pa"][pre_p] - d_pre["pore_pressure_gauge_pa"][pre_d])
        / max(np.linalg.norm(d_pre["pore_pressure_gauge_pa"][pre_d]), 1e-300))

    if dumux_first_saturation_time is not None and python_first_saturation_time is not None:
        metrics["dumux_first_saturation_time_s"] = float(dumux_first_saturation_time)
        metrics["python_interpolated_first_saturation_time_s"] = float(python_first_saturation_time)
        metrics["first_saturation_time_relative_difference"] = _relative(
            python_first_saturation_time, dumux_first_saturation_time)

    if replay is not None:
        replay_end = min(pre, replay["time_s"][-1])
        ri, di = _at(replay, replay_end), _at(d_pre, replay_end)
        metrics["flux_replay_pre_saturation_x_max_abs_difference"] = float(np.max(np.abs(
            replay["pore_x"][ri] - d_pre["pore_x"][di])))

    if python_summaries:
        for name, summary in python_summaries.items():
            metrics[name + "_salt_balance_relative_residual"] = float(
                summary["salt_balance_relative_residual"])
            metrics[name + "_water_balance_relative_residual"] = float(
                summary["water_balance_relative_residual"])
    for prefix, a, b in (("dumux_dt10_vs_dt5", d10, d5), ("python_dt10_vs_dt5", p10, p5)):
        end = min(a["time_s"][-1], b["time_s"][-1]); ia, ib = _at(a, end), _at(b, end)
        metrics[prefix + "_x_max_abs"] = float(np.max(np.abs(a["pore_x"][ia] - b["pore_x"][ib])))
        metrics[prefix + "_cumulative_evap_relative"] = float(abs(
            a["cumulative_evaporation_kg"][ia] - b["cumulative_evaporation_kg"][ib]) /
            max(abs(b["cumulative_evaporation_kg"][ib]), 1e-300))
        if prefix.startswith("dumux"):
            sa, sb = a["solid_mass_kg"][ia], b["solid_mass_kg"][ib]
        else:
            sa = _solid_mass_python(a, ia)
            sb = _solid_mass_python(b, ib)
        metrics[prefix + "_solid_mass_relative"] = float(abs(sa - sb) / max(abs(sb), 1e-300))
    if p2p5 is not None:
        end = min(p5["time_s"][-1], p2p5["time_s"][-1])
        i5, i25 = _at(p5, end), _at(p2p5, end)
        metrics["python_dt5_vs_dt2p5_x_max_abs"] = float(np.max(np.abs(
            p5["pore_x"][i5] - p2p5["pore_x"][i25])))
        metrics["python_dt5_vs_dt2p5_cumulative_evap_relative"] = _relative(
            p5["cumulative_evaporation_kg"][i5], p2p5["cumulative_evaporation_kg"][i25])
        metrics["python_dt5_vs_dt2p5_solid_mass_relative"] = _relative(
            _solid_mass_python(p5, i5), _solid_mass_python(p2p5, i25))
    return metrics


def validation_criteria(metrics):
    """Machine-readable evaluation; missing evidence is never silently treated as a pass."""
    checks = {}
    def add(name, value, limit, relation="lt"):
        passed = value < limit if relation == "lt" else value == limit
        checks[name] = {"value": value, "limit": limit, "relation": relation, "pass": bool(passed)}
    for model in ("dumux", "python"):
        prefix = model + "_dt10_vs_dt5"
        add(prefix + "_x", metrics[prefix + "_x_max_abs"], 1e-4)
        add(prefix + "_evap", metrics[prefix + "_cumulative_evap_relative"], 0.005)
        add(prefix + "_solid", metrics[prefix + "_solid_mass_relative"], 0.01)
    checks["dumux_completed_40000s_and_newton_converged"] = {
        "value": bool(metrics["dumux_completed_40000s_and_newton_converged"]),
        "required": True,
        "pass": bool(metrics["dumux_completed_40000s_and_newton_converged"]),
    }
    checks["dumux_reconstructed_mass_ledger"] = {
        "value": max(metrics["dumux_water_ledger_max_abs_relative_residual_0_1000s"],
                     metrics["dumux_salt_ledger_max_abs_relative_residual_0_1000s"]),
        "limit": 1e-6,
        "status": "failed on 0-1000 s diagnostic; full failed runs predate signed diagnostics",
        "pass": False,
    }
    if "python_dt5_vs_dt2p5_x_max_abs" in metrics:
        add("python_dt5_vs_dt2p5_x", metrics["python_dt5_vs_dt2p5_x_max_abs"], 1e-4)
        add("python_dt5_vs_dt2p5_evap",
            metrics["python_dt5_vs_dt2p5_cumulative_evap_relative"], 0.005)
        add("python_dt5_vs_dt2p5_solid",
            metrics["python_dt5_vs_dt2p5_solid_mass_relative"], 0.01)
    add("self_driven_pre_saturation_x",
        metrics["self_driven_pre_saturation_x_max_abs_difference"], 5e-4)
    add("self_driven_pre_saturation_evap",
        metrics["self_driven_pre_saturation_cumulative_evap_relative_difference"], 0.02)
    if "flux_replay_pre_saturation_x_max_abs_difference" in metrics:
        add("flux_replay_pre_saturation_x",
            metrics["flux_replay_pre_saturation_x_max_abs_difference"], 1e-4)
    if "first_saturation_time_relative_difference" in metrics:
        add("first_saturation_time", metrics["first_saturation_time_relative_difference"], 0.03)
    add("flow_direction", 1.0 - metrics[
        "pre_saturation_effective_flow_direction_agreement_fraction"], 0.0, relation="eq")
    add("flow_weighted_relative", metrics["pre_saturation_flow_weighted_relative_error"], 0.02)
    add("pressure_weighted_relative",
        metrics["pre_saturation_pressure_weighted_relative_error"], 0.02)
    add("final_solid_mass", metrics["final_solid_mass_relative_difference"], 0.05)
    add("final_pore_solid_fraction",
        metrics["final_pore_solid_fraction_max_abs_difference"], 0.02)
    add("final_min_throat_radius",
        metrics["final_min_throat_radius_relative_difference"], 0.02)
    if "python_dt5_salt_balance_relative_residual" in metrics:
        add("python_salt_ledger", abs(metrics["python_dt5_salt_balance_relative_residual"]), 1e-8)
        add("python_water_ledger", abs(metrics["python_dt5_water_balance_relative_residual"]), 1e-8)
    checks["overall_go"] = bool(all(item["pass"] for item in checks.values()))
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args(argv)
    root = args.run_root
    normalized = root / "normalized"; normalized.mkdir(exist_ok=True)
    figures = root / "figures"; figures.mkdir(exist_ok=True)
    tables = root / "tables"; tables.mkdir(exist_ok=True)
    d10, _ = normalize_dumux(root / "dumux_dt10" / "raw", normalized / "dumux_dt10.npz")
    d5, replay = normalize_dumux(root / "dumux_dt5" / "raw", normalized / "dumux_dt5.npz")
    d_double = load_double_diagnostics(
        root / "dumux_diagnostics_1000s_1s" / "validation_pores.csv",
        root / "dumux_diagnostics_1000s_1s" / "validation_throats.csv", d5)
    np.savez_compressed(normalized / "dumux_double_diagnostics_0_1000s.npz", **d_double)
    d_plot = merge_time_series(d5, d_double)
    p10 = load_npz(root / "python_dt10" / "python_results.npz")
    p5 = load_npz(root / "python_dt5" / "python_results.npz")
    p2p5_path = root / "python_dt2p5" / "python_results.npz"
    p2p5 = load_npz(p2p5_path) if p2p5_path.exists() else None
    replay_result_path = root / "python_replay_dt5" / "python_results.npz"
    replay_result = load_npz(replay_result_path) if replay_result_path.exists() else None
    summaries = {
        "python_dt10": load_json(root / "python_dt10" / "summary.json"),
        "python_dt5": load_json(root / "python_dt5" / "summary.json"),
        "python_dt2p5": load_json(root / "python_dt2p5" / "summary.json"),
    }
    common_last = min(d5["time_s"][-1], p5["time_s"][-1])
    times = [1000.0, 940.0, float(common_last)]
    pressure_figure(d_plot, p5, times, figures / "pressure_network_dumux_python")
    pressure_range_figure(d5, p5, figures / "pressure_range_vs_time", dumux_direct=d_double)
    write_snapshot_csv(d_plot, p5, times, tables)
    pressure_single_figure(p5, 40000.0, figures / "pressure_network_python_40000s",
                           "Python only, t=40000 s (DuMuX unavailable)")
    ledger_metrics = reconstruct_dumux_ledger(
        root / "dumux_diagnostics_1000s_1s" / "validation_pores.csv",
        root / "dumux_diagnostics_1000s_1s" / "validation_throats.csv",
        root / "dumux_diagnostics_1000s_1s" / "evaporation.csv",
        tables / "dumux_reconstructed_ledger_0_1000s.csv")
    metrics = comparison_metrics(
        d10, d5, p10, p5, p2p5=p2p5, replay=replay_result, d_pre=d_double,
        dumux_first_saturation_time=first_solubility_crossing(d_double),
        python_first_saturation_time=summaries["python_dt5"]["first_saturation_time_s"],
        python_summaries=summaries)
    metrics.update({
        "requested_end_s": 40000.0,
        "dumux_dt5_last_converged_time_s": 27340.0,
        "dumux_dt10_last_forward_time_s": 27374.0,
        "dumux_last_saved_checkpoint_s": float(d5["time_s"][-1]),
        "dumux_completed_40000s_and_newton_converged": False,
        "geometry_original_pores": int(len(p5["dgf_vertices"])),
        "geometry_dynamic_pores": int(len(p5["file_vertex_index"])),
        "geometry_reservoir_pores": int(len(p5["reservoir_file_indices"])),
        "geometry_internal_throats": int(len(p5["lat_file_edge_index"])),
        "geometry_reservoir_throats": int(np.count_nonzero(p5["res_file_edge_index"] >= 0)),
        "geometry_evaporation_pores": int(len(TOP)),
        "dumux_reservoir_x_max_abs_drift_0_1000s": float(np.max(np.abs(
            d_double["pore_x"][:, RESERVOIR] - X_RES))),
        "dumux_reservoir_sliq_min_0_1000s": float(np.min(
            d_double["pore_liquid_saturation"][:, RESERVOIR])),
    })
    metrics.update(ledger_metrics)
    metrics["pressure_figure_times_s"] = times
    metrics["flux_replay_csv"] = str(replay)
    (root / "comparison_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    criteria = validation_criteria(metrics)
    (root / "validation_criteria.json").write_text(
        json.dumps(criteria, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"metrics": metrics, "criteria": criteria}, indent=2))


if __name__ == "__main__":
    main()

