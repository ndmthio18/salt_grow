"""M4 fixed-material histories and travelling-profile diagnostics.

Run from prat_dynamic_growth with ``python -m prat_growth.m4_analysis``.
"""
from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from .campaign import load_result, save_result
from .dynamic import DynamicALEModel, snapshot_diagnostics
from .material_tracking import (
    WindowResult, evaluate_marker, exact_region_average, marker_history,
    reconstruct, restart_window, state_from_result,
)
from .parameters import NumericalParameters
from .provenance import cache_matches, identity, source_hash


BIRTH_TIMES = (300.0, 3600.0, 10800.0, 18000.0)
PROFILE_TIMES = (300.0, 3600.0, 10800.0, 18000.0, 21600.0)


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def exact_index(result, time: float) -> int:
    index = int(np.argmin(abs(result.time-time)))
    if abs(float(result.time[index])-time) > 1e-7:
        raise ValueError(f"checkpoint {time} not present")
    return index


def l95_and_wake(model, eps: np.ndarray, concentration: np.ndarray, height: float) -> dict:
    p = model.p
    state = model.flux_state(eps, concentration, height)
    reaction = state["reaction"]
    volumes = height*model.widths
    depth = height*(1-model.centers)
    order = np.argsort(depth)
    weights = np.maximum(reaction[order], 0.0)*volumes[order]
    total = float(np.sum(weights))
    if total > 0:
        cumulative = np.r_[0.0, np.cumsum(weights)]
        depth_nodes = np.r_[0.0, depth[order]]
        l95 = float(np.interp(.95*total, cumulative, depth_nodes))
    else:
        l95 = np.nan
    delta = np.sqrt(p.epsilon_ref*p.D_eff_ref/(p.a_v_ref*p.k_reaction*(1-p.C_eq)))
    wake_low = max(3*l95, 5*delta) if np.isfinite(l95) else 5*delta
    wake_high = wake_low+10*delta
    lower, upper = max(0.0, height-wake_high), max(0.0, height-wake_low)
    wake, wake_length = exact_region_average(eps, height*model.faces, lower, upper)
    mask = (depth >= wake_low) & (depth <= wake_high)
    eps_range = float(np.ptp(eps[mask])) if np.count_nonzero(mask) >= 2 else np.nan
    reaction_beyond = float(np.sum(np.maximum(reaction[depth >= wake_low], 0)*volumes[depth >= wake_low]))
    reaction_fraction = reaction_beyond/total if total > 0 else np.nan
    detected = bool(wake_length > .8*(wake_high-wake_low)
                    and np.isfinite(eps_range) and eps_range < .002
                    and reaction_fraction < .01)
    g_internal = float(np.sum(reaction*volumes))
    g_front = p.rho_crystal*(1-p.epsilon_birth)*float(state["velocity"])
    chi = g_internal/(g_internal+g_front)
    predicted = (p.epsilon_birth-chi)/(1-chi)
    return {
        "L95_m": l95, "wake_low_m": wake_low, "wake_high_m": wake_high,
        "epsilon_wake": wake, "wake_epsilon_range": eps_range,
        "positive_reaction_beyond_wake_fraction": reaction_fraction,
        "wake_detected": detected, "G_internal_kg_m2_s": g_internal,
        "G_front_kg_m2_s": g_front, "chi_internal_instantaneous": chi,
        "epsilon_wake_predicted": predicted,
        "wake_relation_error": wake-predicted if np.isfinite(wake) else np.nan,
        "signed_dissolution_flux_kg_m2_s": float(np.sum(np.minimum(reaction, 0)*volumes)),
    }


def spatial_statistics(result, closure: str, n_cells: int) -> list[dict]:
    model = DynamicALEModel(result.physical, result.numerical, closure)
    rows = []
    for i, time in enumerate(result.time):
        eps, c, h = result.epsilon[:, i], result.concentration[:, i], float(result.height[i])
        physical_faces = h*result.faces
        initial_mean, initial_length = exact_region_average(eps, physical_faces, 0, min(result.physical.h0, h))
        new_mean, new_length = exact_region_average(eps, physical_faces, result.physical.h0, h)
        minimum = int(np.argmin(eps))
        state = model.flux_state(eps, c, h)
        rows.append({
            "closure": closure, "n_cells": n_cells, "time_s": float(time),
            "height_m": h, "velocity_m_s": float(state["velocity"]),
            "epsilon_min": float(eps[minimum]),
            "z_min_m": float(h*result.centers[minimum]),
            "d_min_m": float(h*(1-result.centers[minimum])),
            "reaction_at_min_kg_m3_s": float(state["reaction"][minimum]),
            "epsilon_initial_region_mean": initial_mean,
            "initial_region_length_m": initial_length,
            "epsilon_new_region_mean": new_mean,
            "new_region_length_m": new_length,
            "epsilon_whole_mean": float(np.sum(eps*h*result.widths)/h),
            **l95_and_wake(model, eps, c, h),
        })
    return rows


def sample_window(window: WindowResult, stride: int) -> WindowResult:
    return WindowResult(window.time[::stride], window.states[:, ::stride],
                        window.accepted_time, window.nfev, window.njev, window.nlu,
                        window.warnings)


def completion_time(rows: list[dict], fraction: float = .95) -> float:
    start, final = rows[0]["epsilon"], rows[-1]["epsilon"]
    drop = start-final
    if drop <= 0:
        return np.nan
    threshold = fraction*drop
    for row in rows:
        if start-row["epsilon"] >= threshold:
            return float(row["material_age_s"])
    return np.nan


def run_long_b1_480(out: Path, template) -> object:
    p = template.physical
    n = template.numerical.with_updates(n_cells=480)
    name = "B1_480_long"
    meta, data = out/"runs"/(name+".json"), out/"runs"/(name+".npz")
    expected = identity(p, n, "B1")
    if meta.exists() and data.exists():
        record = json.loads(meta.read_text(encoding="utf-8"))
        # Analysis/plotting modules change the package-wide hash without changing
        # the solver. Reuse only after configuration and NPZ byte verification.
        if (record.get("identity", {}).get("config_sha256") == expected["config_sha256"]
                and record.get("data_sha256") == file_sha(data)):
            return load_result(out, name)
    result = DynamicALEModel(p, n, "B1").run()
    if result.status != 0 or result.time[-1] < n.t_end:
        raise RuntimeError(result.message)
    save_result(out, name, result)
    return result


def analyse_windows(runs: dict, out: Path, delta_ref: float, initial_sources: dict):
    histories, balances, summaries, window_records, restart_checks = [], [], [], [], []
    for (closure, n_cells), result in runs.items():
        model = DynamicALEModel(result.physical, result.numerical, closure)
        marker_specs = [(time, f"born_{int(time)}s", None) for time in BIRTH_TIMES]
        # Initial-structure markers are secondary controls. Use the M3 matched
        # quasi-steady concentration initializations to avoid re-solving the
        # millisecond concentration startup; M3 showed the structural paths
        # converge after that startup. Newborn cohorts use both grids.
        if n_cells == 640:
            marker_specs.append((0.0, "initial_group", "initial"))
        for start, label, marker_type in marker_specs:
            source_result = initial_sources[closure] if marker_type == "initial" else result
            index = exact_index(source_result, start)
            y0 = state_from_result(source_result, index)
            eps0, c0, h0, _ = model.unpack(y0)
            state0 = model.flux_state(eps0, c0, h0)
            tau_delta = delta_ref/float(state0["velocity"])
            coarse_target = min(1.0, tau_delta/50)
            duration = max(300.0, 10*tau_delta)
            coarse_count = int(np.ceil(duration/coarse_target))
            fine_dt = duration/(2*coarse_count)
            window = restart_window(model, start, y0, duration, fine_dt, max_step=5.0,
                                    rtol=1e-8 if start == 0 else 1e-9)
            print(f"window {closure} n={n_cells} t={start:g}s: {window.accepted_time.size-1} accepted steps", flush=True)
            refined_window = None
            if closure == "B2" and n_cells == 640 and start == 3600.0:
                refined_window = restart_window(model, start, y0, duration, fine_dt, max_step=1.0)
            window_records.append({
                "closure": closure, "n_cells": n_cells, "start_time_s": start,
                "duration_s": duration, "fine_output_interval_s": fine_dt,
                "coarse_output_interval_s": 2*fine_dt,
                "max_step_s": 5.0,
                "accepted_steps": window.accepted_time.size-1, "nfev": window.nfev,
                "njev": window.njev, "nlu": window.nlu, "warnings": window.warnings,
                "accepted_step_min_s": float(np.min(np.diff(window.accepted_time))),
                "accepted_step_max_s": float(np.max(np.diff(window.accepted_time))),
            })
            if refined_window is not None:
                window_records.append({
                    "closure": closure, "n_cells": n_cells, "start_time_s": start,
                    "duration_s": duration, "fine_output_interval_s": fine_dt,
                    "coarse_output_interval_s": 2*fine_dt, "max_step_s": 1.0,
                    "accepted_steps": refined_window.accepted_time.size-1,
                    "nfev": refined_window.nfev, "njev": refined_window.njev,
                    "nlu": refined_window.nlu, "warnings": refined_window.warnings,
                    "accepted_step_min_s": float(np.min(np.diff(refined_window.accepted_time))),
                    "accepted_step_max_s": float(np.max(np.diff(refined_window.accepted_time))),
                })
            end_index = exact_index(source_result, start+300.0)
            reference = state_from_result(source_result, end_index)
            restart_checks.append({
                "closure": closure, "n_cells": n_cells, "start_time_s": start,
                "height_difference_m": float(window.states[2*n_cells, -1]-reference[2*n_cells]),
                "epsilon_linf": float(np.max(abs(window.states[:n_cells, -1]-reference[:n_cells]))),
                "concentration_excess_linf": float(np.max(abs(window.states[n_cells:2*n_cells, -1]-reference[n_cells:2*n_cells]))),
            })
            if marker_type == "initial":
                positions = [
                    (result.physical.h0-delta_ref, "initial_h0_minus_delta"),
                    (result.physical.h0-3*delta_ref, "initial_h0_minus_3delta"),
                    (result.physical.h0/2, "initial_midheight"),
                ]
            else:
                positions = [(h0, label)]
            for z0, marker in positions:
                initial_eps = evaluate_marker(model, y0, z0, "linear")["epsilon"]
                fine_linear = marker_history(model, window, z0, None if marker_type == "initial" else start,
                                             "linear", initial_eps)
                fine_pc = marker_history(model, window, z0, None if marker_type == "initial" else start,
                                         "piecewise_constant", initial_eps)
                coarse = sample_window(window, 2)
                coarse_linear = marker_history(model, coarse, z0, None if marker_type == "initial" else start,
                                               "linear", initial_eps)
                refined_linear = (marker_history(model, refined_window, z0, start,
                                                  "linear", initial_eps)
                                  if refined_window is not None else None)
                for row, pc in zip(fine_linear, fine_pc):
                    row.update({"closure": closure, "n_cells": n_cells,
                                "marker": marker, "marker_type": marker_type or "newborn",
                                "initial_concentration_profile": "quasi_steady" if marker_type == "initial" else "checkpoint",
                                "epsilon_piecewise_constant": pc["epsilon"],
                                "piecewise_material_balance_error": pc["material_balance_error"]})
                    histories.append(row)
                    balances.append({key: row[key] for key in (
                        "closure", "n_cells", "marker", "marker_type", "time_s",
                        "birth_time_s", "material_age_s", "z0_m", "epsilon",
                        "epsilon_reaction", "material_balance_error",
                        "epsilon_piecewise_constant", "piecewise_material_balance_error",
                        "cumulative_signed_reaction_kg_m3",
                        "cumulative_positive_precipitation_kg_m3",
                        "cumulative_dissolution_kg_m3")})
                drop = fine_linear[0]["epsilon"]-fine_linear[-1]["epsilon"]
                tail = fine_linear[int(.8*len(fine_linear)):]
                reaction_scale = max(abs(row["reaction_kg_m3_s"]) for row in fine_linear)
                platform = (np.ptp([row["epsilon"] for row in tail]) < 5e-4
                            and abs(fine_linear[-1]["reaction_kg_m3_s"])
                            < 1e-3*max(reaction_scale, 1e-30))
                coarse_t95 = completion_time(coarse_linear)
                fine_t95 = completion_time(fine_linear)
                summaries.append({
                    "closure": closure, "n_cells": n_cells, "marker": marker,
                    "marker_type": marker_type or "newborn",
                    "initial_concentration_profile": "quasi_steady" if marker_type == "initial" else "checkpoint",
                    "birth_time_s": "" if marker_type == "initial" else start,
                    "z0_m": z0, "final_observed_age_s": fine_linear[-1]["material_age_s"],
                    "final_epsilon": fine_linear[-1]["epsilon"],
                    "final_epsilon_piecewise": fine_pc[-1]["epsilon"],
                    "platform_detected": bool(platform),
                    "major_densification_t95_s": fine_t95,
                    "cumulative_signed_reaction_kg_m3": fine_linear[-1]["cumulative_signed_reaction_kg_m3"],
                    "cumulative_positive_precipitation_kg_m3": fine_linear[-1]["cumulative_positive_precipitation_kg_m3"],
                    "cumulative_dissolution_kg_m3": fine_linear[-1]["cumulative_dissolution_kg_m3"],
                    "final_material_balance_error": fine_linear[-1]["material_balance_error"],
                    "final_balance_error_relative_to_drop": abs(fine_linear[-1]["material_balance_error"])/max(abs(drop), 1e-30),
                    "max_abs_material_balance_error": max(abs(row["material_balance_error"]) for row in fine_linear),
                    "relative_balance_error_to_drop": max(abs(row["material_balance_error"]) for row in fine_linear)/max(abs(drop), 1e-30),
                    "material_balance_target_met_all_ages": max(abs(row["material_balance_error"]) for row in fine_linear)/max(abs(drop), 1e-30) <= .01,
                    "final_reconstruction_difference": fine_linear[-1]["epsilon"]-fine_pc[-1]["epsilon"],
                    "output_halving_reaction_difference_relative": abs(fine_linear[-1]["cumulative_signed_reaction_kg_m3"]-coarse_linear[-1]["cumulative_signed_reaction_kg_m3"])/max(abs(fine_linear[-1]["cumulative_signed_reaction_kg_m3"]), 1e-30),
                    "output_halving_t95_difference_relative": abs(fine_t95-coarse_t95)/max(abs(fine_t95), 1e-30),
                    "integration_step_final_epsilon_difference": "" if refined_linear is None else abs(fine_linear[-1]["epsilon"]-refined_linear[-1]["epsilon"]),
                    "integration_step_reaction_difference_relative": "" if refined_linear is None else abs(fine_linear[-1]["cumulative_signed_reaction_kg_m3"]-refined_linear[-1]["cumulative_signed_reaction_kg_m3"])/max(abs(refined_linear[-1]["cumulative_signed_reaction_kg_m3"]), 1e-30),
                    "integration_step_t95_difference_relative": "" if refined_linear is None else abs(fine_t95-completion_time(refined_linear))/max(abs(completion_time(refined_linear)), 1e-30),
                    "space_converged": "pending", "time_converged": "pending",
                })
    # Cross-grid status uses identical closure and material label/age.
    for row in summaries:
        partner = next((item for item in summaries if item["closure"] == row["closure"]
                        and item["marker"] == row["marker"] and item["n_cells"] != row["n_cells"]), None)
        if partner is None:
            row["cross_grid_final_epsilon_difference"] = ""
            row["space_converged"] = "not_assessed"
            row["time_converged"] = (row["output_halving_reaction_difference_relative"] <= .01
                                     and row["output_halving_t95_difference_relative"] <= .01)
            continue
        spatial_difference = abs(row["final_epsilon"]-partner["final_epsilon"])
        row["cross_grid_final_epsilon_difference"] = spatial_difference
        row["space_converged"] = spatial_difference <= .001
        row["time_converged"] = (row["output_halving_reaction_difference_relative"] <= .01
                                 and row["output_halving_t95_difference_relative"] <= .01)
    return histories, balances, summaries, window_records, restart_checks


def travelling_profiles(runs: dict, delta_ref: float):
    rows, metrics = [], []
    depth_grid = np.linspace(0, 15e-6, 301)
    profile_cache = {}
    for (closure, n_cells), result in runs.items():
        model = DynamicALEModel(result.physical, result.numerical, closure)
        for time in PROFILE_TIMES:
            i = exact_index(result, time)
            eps, c, h = result.epsilon[:, i], result.concentration[:, i], float(result.height[i])
            state = model.flux_state(eps, c, h)
            eps_d, c_d = [], []
            for depth in depth_grid:
                z = h-depth
                eps_d.append(reconstruct(z, h, model.faces, model.centers, eps,
                                         result.physical.epsilon_birth, "linear"))
                c_d.append(reconstruct(z, h, model.faces, model.centers, c,
                                       state["c_top"], "linear"))
            eps_d, c_d = np.array(eps_d), np.array(c_d)
            av = np.full_like(eps_d, result.physical.a_v_ref) if closure == "B1" else 36*(1-eps_d)/result.physical.d_sub
            reaction = result.physical.rho_solution*av*result.physical.k_reaction*(c_d-result.physical.C_eq)
            deps_dd = np.gradient(eps_d, depth_grid, edge_order=2)
            if time == result.time[-1]:
                i0, i1 = i-1, i
            else:
                i0, i1 = max(0, i-1), min(len(result.time)-1, i+1)
            eps_other = []
            for j in (i0, i1):
                hj = float(result.height[j])
                sj = model.flux_state(result.epsilon[:, j], result.concentration[:, j], hj)
                eps_other.append(np.array([reconstruct(hj-depth, hj, model.faces, model.centers,
                                                       result.epsilon[:, j], result.physical.epsilon_birth, "linear")
                                                   for depth in depth_grid]))
            deps_dt_d = (eps_other[1]-eps_other[0])/(result.time[i1]-result.time[i0])
            source = reaction/result.physical.rho_crystal
            steady_residual = float(state["velocity"])*deps_dd+source
            general_residual = deps_dt_d+steady_residual
            steady_integral = result.physical.epsilon_birth-cumulative_trapezoid(
                reaction/(result.physical.rho_crystal*float(state["velocity"])),
                depth_grid, initial=0.0)
            steady_integral_error = eps_d-steady_integral
            scale = max(float(np.max(abs(source))), 1e-30)
            profile_cache[(closure, n_cells, time)] = (eps_d, reaction)
            for k, depth in enumerate(depth_grid):
                rows.append({
                    "closure": closure, "n_cells": n_cells, "time_s": time,
                    "distance_from_front_m": depth,
                    "distance_over_delta_ref": depth/delta_ref,
                    "epsilon": eps_d[k], "supersaturation": c_d[k]-result.physical.C_eq,
                    "reaction_kg_m3_s": reaction[k],
                    "V_depsilon_dd_s_inv": float(state["velocity"])*deps_dd[k],
                    "partial_epsilon_partial_t_at_d_s_inv": deps_dt_d[k],
                    "steady_relation_residual_s_inv": steady_residual[k],
                    "general_relation_residual_s_inv": general_residual[k],
                    "steady_integral_epsilon": steady_integral[k],
                    "steady_integral_error": steady_integral_error[k],
                })
            mask = depth_grid <= 5e-6
            metrics.append({
                "closure": closure, "n_cells": n_cells, "time_s": time,
                "steady_residual_linf_relative": float(np.max(abs(steady_residual[mask]))/scale),
                "general_residual_linf_relative": float(np.max(abs(general_residual[mask]))/scale),
                "explicit_time_term_linf_relative": float(np.max(abs(deps_dt_d[mask]))/scale),
                "steady_integral_error_linf_relative_to_drop": float(
                    np.max(abs(steady_integral_error[mask]))
                    / max(result.physical.epsilon_birth-np.min(eps_d[mask]), 1e-30)),
            })
    for row in metrics:
        eps, reaction = profile_cache[(row["closure"], row["n_cells"], row["time_s"])]
        eps6, reaction6 = profile_cache[(row["closure"], row["n_cells"], 21600.0)]
        mask5, mask10 = depth_grid <= 5e-6, depth_grid <= 10e-6
        row["epsilon_linf_to_6h_5um"] = float(np.max(abs(eps[mask5]-eps6[mask5])))
        row["epsilon_linf_to_6h_10um"] = float(np.max(abs(eps[mask10]-eps6[mask10])))
        row["reaction_linf_to_6h_relative"] = float(np.max(abs(reaction[mask10]-reaction6[mask10]))/max(np.max(abs(reaction6[mask10])), 1e-30))
    return rows, metrics


def save_figure(fig, out: Path, name: str):
    fig.savefig(out/(name+".png"), dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(out/(name+".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_figures(out: Path, histories, balances, stats, profiles):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    colors = {"B1": "#0072B2", "B2": "#D55E00"}
    s640 = [r for r in stats if r["n_cells"] == 640]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for closure in ("B1", "B2"):
        data = [r for r in s640 if r["closure"] == closure]
        t = np.array([r["time_s"] for r in data])/3600
        axes[0].plot(t, [r["epsilon_min"] for r in data], color=colors[closure], label=f"{closure} min")
        axes[0].plot(t, [r["epsilon_wake"] for r in data], "--", color=colors[closure], label=f"{closure} wake")
        axes[1].plot(t, [r["epsilon_initial_region_mean"] for r in data], ":", color=colors[closure], label=f"{closure} initial")
        axes[1].plot(t, [r["epsilon_new_region_mean"] for r in data], "--", color=colors[closure], label=f"{closure} new")
        axes[1].plot(t, [r["epsilon_whole_mean"] for r in data], color=colors[closure], label=f"{closure} whole")
    for ax in axes:
        ax.set(xlabel="Time (h)", ylabel="Porosity")
        ax.legend(ncol=2)
    save_figure(fig, out, "spatial_statistics")

    new = [r for r in histories if r["n_cells"] == 640 and r["marker_type"] == "newborn"]
    for xkey, name, xlabel in (("time_s", "material_absolute_time", "Absolute time (h)"),
                               ("material_age_s", "material_age", "Material age (s)")):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True, sharey=True)
        for ax, closure in zip(axes, ("B1", "B2")):
            for marker in sorted(set(r["marker"] for r in new)):
                data = [r for r in new if r["closure"] == closure and r["marker"] == marker]
                x = np.array([r[xkey] for r in data])/(3600 if xkey == "time_s" else 1)
                ax.plot(x, [r["epsilon"] for r in data], label=marker.replace("born_", "tb="))
            ax.set(title=closure, xlabel=xlabel, ylabel="Fixed-material porosity")
            ax.legend()
        save_figure(fig, out, name)

    for ykey, name, ylabel in (("supersaturation", "material_supersaturation", "Local excess mass fraction"),
                               ("reaction_kg_m3_s", "material_reaction", r"Local reaction ($kg\,m^{-3}\,s^{-1}$)")):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        for ax, closure in zip(axes, ("B1", "B2")):
            for marker in sorted(set(r["marker"] for r in new)):
                data = [r for r in new if r["closure"] == closure and r["marker"] == marker]
                ax.plot([r["material_age_s"] for r in data], [r[ykey] for r in data], label=marker.replace("born_", "tb="))
            ax.set(title=closure, xlabel="Material age (s)", ylabel=ylabel)
            ax.legend()
        save_figure(fig, out, name)

    bnew = [r for r in balances if r["n_cells"] == 640 and r["marker_type"] == "newborn"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for closure in ("B1", "B2"):
        marker = "born_3600s"
        data = [r for r in bnew if r["closure"] == closure and r["marker"] == marker]
        axes[0].plot([r["material_age_s"] for r in data], [r["epsilon"] for r in data], color=colors[closure], label=f"{closure} field")
        axes[0].plot([r["material_age_s"] for r in data], [r["epsilon_reaction"] for r in data], "--", color=colors[closure], label=f"{closure} reaction integral")
        axes[1].plot([r["material_age_s"] for r in data], [r["material_balance_error"] for r in data], color=colors[closure], label=closure)
    axes[0].set(xlabel="Material age (s)", ylabel="Porosity", title="1 h birth cohort")
    axes[1].set(xlabel="Material age (s)", ylabel=r"$\epsilon_{field}-\epsilon_{reaction}$")
    axes[0].legend(); axes[1].legend()
    save_figure(fig, out, "material_balance")

    p640 = [r for r in profiles if r["n_cells"] == 640]
    for ykey, name, ylabel in (("epsilon", "epsilon_front_profiles", "Porosity"),
                               ("reaction_kg_m3_s", "reaction_front_profiles", r"Reaction ($kg\,m^{-3}\,s^{-1}$)")):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        for ax, closure in zip(axes, ("B1", "B2")):
            for time in PROFILE_TIMES:
                data = [r for r in p640 if r["closure"] == closure and r["time_s"] == time and r["distance_from_front_m"] <= 8e-6]
                ax.plot(np.array([r["distance_from_front_m"] for r in data])*1e6,
                        [r[ykey] for r in data], label=f"{time/3600:g} h" if time != 300 else "5 min")
            ax.set(title=closure, xlabel="Distance below front (micrometres)", ylabel=ylabel)
            ax.legend()
        save_figure(fig, out, name)


def run(output: Path = Path("outputs/m4")):
    output.mkdir(parents=True, exist_ok=True)
    (output/"runs").mkdir(exist_ok=True)
    m3 = Path("outputs/m3")
    b2_640 = load_result(m3, "grid_640")
    runs = {
        ("B2", 640): b2_640,
        ("B2", 480): load_result(m3, "grid_480"),
        ("B1", 640): load_result(m3, "T0"),
    }
    runs[("B1", 480)] = run_long_b1_480(output, b2_640)
    p = b2_640.physical
    delta_ref = np.sqrt(p.epsilon_ref*p.D_eff_ref/(p.a_v_ref*p.k_reaction*(1-p.C_eq)))
    stats = []
    for (closure, n_cells), result in runs.items():
        stats.extend(spatial_statistics(result, closure, n_cells))
    initial_sources = {"B1": load_result(m3, "T1"), "B2": load_result(m3, "T3")}
    histories, balances, cohorts, windows, restart_checks = analyse_windows(
        runs, output, delta_ref, initial_sources)
    profiles, profile_metrics = travelling_profiles(runs, delta_ref)
    write_csv(output/"material_histories.csv", histories)
    write_csv(output/"cohort_summary.csv", cohorts)
    write_csv(output/"spatial_statistics.csv", stats)
    write_csv(output/"material_balance.csv", balances)
    write_csv(output/"travelling_profile_checks.csv", profiles)
    write_csv(output/"travelling_profile_summary.csv", profile_metrics)
    make_figures(output, histories, balances, stats, profiles)
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    manifest = {
        "reference_commit": "a37759f90f0aa92b3fb37df7b4106cf711da2d02",
        "analysis_head_before_m4_commit": head,
        "code_sha256": source_hash(), "physical": asdict(p),
        "delta_ref_m": delta_ref,
        "tau_delta_at_6h_s": delta_ref/snapshot_diagnostics(b2_640, -1)["velocity"],
        "tau_95_at_6h_s": l95_and_wake(DynamicALEModel(p, b2_640.numerical, "B2"), b2_640.epsilon[:, -1], b2_640.concentration[:, -1], b2_640.height[-1])["L95_m"]/snapshot_diagnostics(b2_640, -1)["velocity"],
        "input_files": {path.as_posix(): file_sha(path)
                        for path in sorted((m3/"runs").glob("*"))
                        if path.is_file() and path.stem in
                        ("grid_480", "grid_640", "T0", "T1", "T3")},
        "windows": windows, "restart_checks": restart_checks,
        "diagnostic_thresholds": {
            "cross_grid_final_epsilon_absolute": .001,
            "material_balance_relative_to_drop": .01,
            "output_halving_reaction_and_t95_relative": .01,
        },
        "software": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "notes": [
            "Dense output interval is not the BDF accepted step size.",
            "Linear and piecewise-constant reconstructions are both retained.",
            "No physical parameter, birth rule, hydraulic feedback or reaction closure was changed.",
        ],
    }
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"M4 complete: {len(histories)} material samples, {len(profiles)} profile samples")


if __name__ == "__main__":
    run()
