from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib
import numpy as np
import scipy

from .analytic import solve_stage_a_finite_volume
from .dynamic import DynamicALEModel, DynamicResult, snapshot_diagnostics
from .parameters import NumericalParameters, PhysicalParameters
from .plotting import make_required_figures
from .provenance import identity, cache_matches
from hashlib import sha256


def load_config(path: Path) -> tuple[dict, PhysicalParameters, NumericalParameters]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    p = PhysicalParameters(**raw["physical"])
    n = NumericalParameters(**raw["numerical"])
    return raw, p, n


def write_environment(output: Path, raw: dict) -> None:
    record = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "paper": raw["paper"],
    }
    (output / "environment.json").write_text(json.dumps(record, indent=2), encoding="utf-8")


def stage_a(output: Path, p: PhysicalParameters) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for n_cells in (80, 160, 320, 640, 1280):
        result = solve_stage_a_finite_volume(p, p.h0, n_cells, grid_power=4.0)
        rows.append({
            "n_cells": n_cells,
            "grid_power": 4.0,
            "Pe": result.groups.Pe,
            "Da": result.groups.Da,
            "Da_top": result.groups.Da_top,
            "Da_star": result.groups.Da_star,
            "delta_C_m": result.groups.delta_C,
            "velocity_numerical_m_s": result.velocity,
            "velocity_general_m_s": result.velocity_exact,
            "velocity_large_Da_m_s": result.velocity_large_da,
            "relative_velocity_error": result.relative_velocity_error,
            "max_profile_error": result.max_profile_error,
        })
    with (output / "stage_a_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def save_result(output: Path, case_id: str, result: DynamicResult) -> None:
    run_dir = output / "runs"
    run_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        run_dir / f"{case_id}.npz",
        time=result.time,
        faces=result.faces,
        centers=result.centers,
        widths=result.widths,
        epsilon=result.epsilon,
        concentration=result.concentration,
        height=result.height,
        top_solid=result.top_solid,
        internal_solid=result.internal_solid,
        bottom_salt=result.bottom_salt,
        bottom_solution=result.bottom_solution,
        evaporated=result.evaporated,
    )
    metadata = {
        "case_id": case_id,
        "closure": result.closure,
        "physical": asdict(result.physical),
        "numerical": asdict(result.numerical),
        "solver_status": result.status,
        "solver_message": result.message,
        "initial_profile": result.initial_profile,
        "solver_warnings": result.solver_warnings,
        "identity": identity(result.physical, result.numerical, result.closure, result.initial_profile),
        "data_sha256": sha256((run_dir / f"{case_id}.npz").read_bytes()).hexdigest(),
    }
    (run_dir / f"{case_id}.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def load_result(output: Path, case_id: str) -> DynamicResult:
    run_dir = output / "runs"
    metadata = json.loads((run_dir / f"{case_id}.json").read_text(encoding="utf-8"))
    arrays = np.load(run_dir / f"{case_id}.npz")
    return DynamicResult(
        closure=metadata["closure"],
        physical=PhysicalParameters(**metadata["physical"]),
        numerical=NumericalParameters(**metadata["numerical"]),
        faces=arrays["faces"],
        centers=arrays["centers"],
        widths=arrays["widths"],
        time=arrays["time"],
        epsilon=arrays["epsilon"],
        concentration=arrays["concentration"],
        height=arrays["height"],
        top_solid=arrays["top_solid"],
        internal_solid=arrays["internal_solid"],
        bottom_salt=arrays["bottom_salt"],
        bottom_solution=arrays["bottom_solution"],
        evaporated=arrays["evaporated"],
        status=metadata["solver_status"],
        message=metadata["solver_message"],
        initial_profile=metadata.get("initial_profile", "equilibrium"),
        solver_warnings=metadata.get("solver_warnings", []),
    )


def result_rows(case_id: str, result: DynamicResult) -> list[dict[str, float | str | int]]:
    rows = []
    for i, time in enumerate(result.time):
        d = snapshot_diagnostics(result, i)
        rows.append({
            "case_id": case_id,
            "closure": result.closure,
            "solver_status": result.status,
            "solver_message": result.message,
            "n_cells": result.numerical.n_cells,
            "max_step_s": result.numerical.max_step,
            "time_s": time,
            "height_m": result.height[i],
            "velocity_m_s": d["velocity"],
            "epsilon_min": d["epsilon_min"],
            "epsilon_top_cell": d["epsilon_top_cell"],
            "concentration_min": d["concentration_min"],
            "concentration_max": d["concentration_max"],
            "top_solid_kg_m2": result.top_solid[i],
            "internal_solid_kg_m2": result.internal_solid[i],
            "dissolved_salt_kg_m2": d["dissolved_salt"],
            "bottom_salt_kg_m2": result.bottom_salt[i],
            "salt_residual_normalized": d["salt_residual_normalized"],
            "bulk_residual_normalized": d["bulk_residual_normalized"],
            "L95_m": d["L95"],
            "delta_C_m": d["delta_C"],
            "cells_across_delta": d["cells_across_delta"],
            "Pe": d["Pe"],
            "Da": d["Da"],
            "Da_top": d["Da_top"],
            "Da_star": d["Da_star"],
        })
    return rows


def run_case(output: Path, case_id: str, p: PhysicalParameters,
             n: NumericalParameters, closure: str,
             reuse_existing: bool = False,
             allow_stopped: bool = False) -> tuple[DynamicResult, list[dict]]:
    metadata_path = output / "runs" / f"{case_id}.json"
    data_path = output / "runs" / f"{case_id}.npz"
    if reuse_existing and cache_matches(metadata_path, data_path, identity(p, n, closure)):
        result = load_result(output, case_id)
    else:
        result = DynamicALEModel(p, n, closure).run()
    if result.time[-1] < n.t_end and not allow_stopped:
        raise RuntimeError(f"{case_id} stopped at {result.time[-1]} s: {result.message}")
    save_result(output, case_id, result)
    return result, result_rows(case_id, result)


def write_rows(output: Path, rows: list[dict]) -> None:
    if not rows:
        return
    time_series_path = output / "run_time_series.csv"
    if time_series_path.exists():
        with time_series_path.open("r", newline="", encoding="utf-8") as stream:
            existing = list(csv.DictReader(stream))
        replaced = {str(row["case_id"]) for row in rows}
        rows = [row for row in existing if row["case_id"] not in replaced] + rows
    fieldnames = list(dict.fromkeys(key for row in rows for key in row.keys()))
    with time_series_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    terminal = [row for i, row in enumerate(rows)
                if i == len(rows) - 1 or rows[i + 1]["case_id"] != row["case_id"]]
    with (output / "run_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(terminal)


def baseline(output: Path, p: PhysicalParameters, n: NumericalParameters):
    results: dict[str, DynamicResult] = {}
    rows: list[dict] = []
    for closure in ("B1", "B2"):
        result, case_rows = run_case(output, f"baseline_{closure}", p, n, closure)
        results[closure] = result
        rows.extend(case_rows)
    return results, rows


def refinement(output: Path, p: PhysicalParameters, n: NumericalParameters):
    rows: list[dict] = []
    for cells in (80, 120, 160):
        n_case = n.with_updates(n_cells=cells, grid_power=2.0)
        _, case_rows = run_case(output, f"grid_B2_N{cells}", p, n_case, "B2",
                                reuse_existing=True)
        rows.extend(case_rows)
    for power, label in ((2.5, "2p5"), (3.0, "3p0")):
        n_clustered = n.with_updates(grid_power=power)
        _, case_rows = run_case(output, f"grid_B2_N160_power{label}", p, n_clustered, "B2",
                                reuse_existing=True)
        rows.extend(case_rows)
    for max_step in (120.0, 30.0):
        n_case = n.with_updates(max_step=max_step, grid_power=3.0)
        _, case_rows = run_case(output, f"time_B2_power3_dt{int(max_step)}", p, n_case, "B2",
                                reuse_existing=True)
        rows.extend(case_rows)
    return rows


def sensitivity(output: Path, p: PhysicalParameters, n: NumericalParameters):
    rows: list[dict] = []
    cases = [
        ("S1_dsub10um", p.with_updates(d_sub=10e-6), n),
        ("S1_dsub100um", p.with_updates(d_sub=100e-6), n),
        ("S2_evap0p1cmday", p.with_updates(U_evap=0.1e-2 / 86400.0), n),
        # Lower k_r produces a thicker reaction layer; quadratic clustering is
        # sufficient and avoids unnecessary tiny-cell stiffness.
        ("S3_kr2p3e6_power2", p.with_updates(k_reaction=2.3e-6), n.with_updates(grid_power=2.0)),
        ("S3_kr2p3e5_power2", p.with_updates(k_reaction=2.3e-5), n.with_updates(grid_power=2.0)),
        ("S3_kr2p3e4_power2", p.with_updates(k_reaction=2.3e-4), n.with_updates(grid_power=2.0)),
        ("S6_h0p5mm_power2", p.with_updates(h0=0.5e-3), n.with_updates(grid_power=2.0)),
        ("S6_h2mm_power2p5", p.with_updates(h0=2e-3), n.with_updates(grid_power=2.5)),
    ]
    for case_id, p_case, n_case in cases:
        _, case_rows = run_case(output, case_id, p_case, n_case, "B2",
                                reuse_existing=True, allow_stopped=True)
        rows.extend(case_rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("stage-a", "baseline", "refinement", "sensitivity", "all"))
    parser.add_argument("--config", type=Path, default=Path("configs/baseline.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw, p, n = load_config(args.config)
    write_environment(args.output, raw)
    rows: list[dict] = []
    results = None
    if args.mode in ("stage-a", "all"):
        stage_a(args.output, p)
    if args.mode in ("baseline", "all"):
        results, baseline_rows = baseline(args.output, p, n)
        rows.extend(baseline_rows)
        make_required_figures(results, args.output / "figures")
    if args.mode in ("refinement", "all"):
        rows.extend(refinement(args.output, p, n))
    if args.mode in ("sensitivity", "all"):
        rows.extend(sensitivity(args.output, p, n))
    write_rows(args.output, rows)


if __name__ == "__main__":
    main()
