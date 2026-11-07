"""Run the focused Phase 2A.3b-4R front-error attribution audit."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
from time import perf_counter
import tracemalloc
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from run_phase2a3b4_validation import (  # noqa: E402
    FAMILIES,
    FRONT_NAMES,
    NORMALS,
    build_planar_system,
    central_plane_mask,
)
from salt_phase2.dynamic_interface_flux import DynamicInterfaceFluxSystem, FrontRegion  # noqa: E402
from salt_phase2.front_error_attribution import (  # noqa: E402
    controlled_variant_weights,
    fixed_cylinder_mask,
    fixed_region_inventory,
    front_height_statistics,
    interface_filter_audit,
    reservoir_rows,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def translated_plane_mask(cells: int, normal: np.ndarray, offset_m: float, dx_m: float) -> np.ndarray:
    coordinates = (-0.5 * cells + 0.5 + np.arange(cells, dtype=float)) * dx_m
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij", sparse=True)
    return np.asarray(normal[0] * x + normal[1] * y + normal[2] * z <= offset_m)


def single_front_system(
    mask: np.ndarray,
    normal: np.ndarray,
    dx_m: float,
    fitting_radius_m: float,
    radius_m: float,
    minimum_alignment: float,
    density: float,
) -> DynamicInterfaceFluxSystem:
    return DynamicInterfaceFluxSystem(
        [FrontRegion(
            name="front",
            solid_mask=mask,
            reference_normal=normal,
            centre_m=np.zeros(3),
            radius_m=radius_m,
            minimum_alignment=minimum_alignment,
        )],
        voxel_size_m=dx_m,
        fitting_radius_m=fitting_radius_m,
        solid_density_kg_m3=density,
        grain_to_lab=np.eye(3),
        alpha=0.0,
        rng_seed=0,
    )


def translated_plane_benchmarks(section: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for voxel_um in section["voxel_sizes_um"]:
        dx = float(voxel_um) * 1e-6
        cells = int(round(float(section["domain_length_um"]) / float(voxel_um)))
        for family, normal in zip(FAMILIES, NORMALS, strict=True):
            initial_mask = translated_plane_mask(cells, normal, 0.0, dx)
            region = fixed_cylinder_mask(
                initial_mask.shape,
                dx,
                normal,
                float(section["measurement_radius_um"]) * 1e-6,
                float(section["axial_half_width_um"]) * 1e-6,
            )
            baseline = single_front_system(
                initial_mask,
                normal,
                dx,
                float(section["fitting_radius_um"]) * 1e-6,
                float(section["growth_radius_um"]) * 1e-6,
                float(section["minimum_normal_alignment"]),
                float(section["solid_density_kg_m3"]),
            )
            base_snapshot = baseline.front_snapshot("front")
            base_height = front_height_statistics(
                baseline, "front", float(section["measurement_radius_um"]) * 1e-6
            )
            for offset_voxels in section["translation_in_voxels"]:
                offset_m = float(offset_voxels) * dx
                mask = translated_plane_mask(cells, normal, offset_m, dx)
                system = single_front_system(
                    mask,
                    normal,
                    dx,
                    float(section["fitting_radius_um"]) * 1e-6,
                    float(section["growth_radius_um"]) * 1e-6,
                    float(section["minimum_normal_alignment"]),
                    float(section["solid_density_kg_m3"]),
                )
                snapshot = system.front_snapshot("front")
                height = front_height_statistics(
                    system, "front", float(section["measurement_radius_um"]) * 1e-6
                )
                gained = int(np.count_nonzero(region & mask & ~initial_mask))
                lost = int(np.count_nonzero(region & initial_mask & ~mask))
                realized_volume = (gained - lost) * dx**3
                reference_area = math.pi * (float(section["measurement_radius_um"]) * 1e-6) ** 2
                current_advance = snapshot.area_weighted_position_m - base_snapshot.area_weighted_position_m
                volume_advance = realized_volume / reference_area
                mean_advance = float(height["mean_m"]) - float(base_height["mean_m"])
                median_advance = float(height["median_m"]) - float(base_height["median_m"])
                rows.append({
                    "voxel_size_um": voxel_um,
                    "orientation_family": family,
                    "translation_voxels": offset_voxels,
                    "true_translation_um": offset_m * 1e6,
                    "current_area_weighted_advance_um": current_advance * 1e6,
                    "fixed_region_realized_volume_advance_um": volume_advance * 1e6,
                    "interior_height_mean_advance_um": mean_advance * 1e6,
                    "interior_height_median_advance_um": median_advance * 1e6,
                    "interior_height_std_um": float(height["std_m"]) * 1e6,
                    "interior_height_iqr_um": float(height["iqr_m"]) * 1e6,
                    "interior_height_p05_um": float(height["p05_m"]) * 1e6,
                    "interior_height_p95_um": float(height["p95_m"]) * 1e6,
                    "current_measurement_bias_um": (current_advance - offset_m) * 1e6,
                    "volume_measurement_bias_um": (volume_advance - offset_m) * 1e6,
                    "height_mean_bias_um": (mean_advance - offset_m) * 1e6,
                    "height_median_bias_um": (median_advance - offset_m) * 1e6,
                    "reconstructed_area_m2": snapshot.reconstructed_area_m2,
                    "binary_realized_volume_m3": realized_volume,
                    "fractional_deposited_volume_m3": 0.0,
                })
    return rows


def growth_case_worker(case: dict[str, object]) -> dict[str, object]:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    started = perf_counter()
    tracemalloc.start()
    system = build_planar_system(
        domain_length_um=float(case["domain_length_um"]),
        front_radius_um=float(case["growth_radius_um"]),
        voxel_size_um=float(case["voxel_size_um"]),
        fitting_radius_um=float(case["fitting_radius_um"]),
        minimum_alignment=float(case["minimum_normal_alignment"]),
        density=float(case["solid_density_kg_m3"]),
        alpha=float(case["alpha"]),
        rng_seed=int(case["rng_seed"]),
    )
    initial_masks = {name: system.fronts[name].solid_mask.copy() for name in FRONT_NAMES}
    initial_snapshots = {name: system.front_snapshot(name) for name in FRONT_NAMES}
    initial_heights = {
        name: front_height_statistics(system, name, float(case["measurement_radius_um"]) * 1e-6)
        for name in FRONT_NAMES
    }
    initial_keys = {
        name: set(system._surface_cache[name].geometry)
        for name in FRONT_NAMES
    }
    regions = {
        name: fixed_cylinder_mask(
            system.fronts[name].solid_mask.shape,
            system.voxel_size_m,
            system.fronts[name].reference_normal,
            float(case["measurement_radius_um"]) * 1e-6,
            float(case["axial_half_width_um"]) * 1e-6,
        )
        for name in FRONT_NAMES
    }
    growth_radius_m = float(case["growth_radius_um"]) * 1e-6
    growth_area = math.pi * growth_radius_m**2
    total_volume = 3.0 * growth_area * float(case["reference_advance_um"]) * 1e-6
    total_mass = total_volume * float(case["solid_density_kg_m3"])
    prior_fraction = 0.0
    stage_rows: list[dict[str, object]] = []
    variant_rows: list[dict[str, object]] = []
    filter_rows: list[dict[str, object]] = []
    update_runtime_s = 0.0
    for stage_fraction in case["stage_fractions"]:
        increment = float(stage_fraction) - prior_fraction
        update = system.advance(total_mass * increment)
        update_runtime_s += update.runtime_s
        prior_fraction = float(stage_fraction)
        expected_advance_um = float(case["reference_advance_um"]) * float(stage_fraction)
        stage_input_volume = total_volume * float(stage_fraction)
        controls_by_alpha = {}
        for alpha in (float(case["alpha"]), float(case["anisotropic_probe_alpha"])):
            controls_by_alpha[alpha] = controlled_variant_weights(
                system, FRONT_NAMES, growth_area, alpha
            )
        for index, (name, family) in enumerate(zip(FRONT_NAMES, FAMILIES, strict=True)):
            snapshot = system.front_snapshot(name)
            height = front_height_statistics(
                system, name, float(case["measurement_radius_um"]) * 1e-6
            )
            inventory = fixed_region_inventory(
                system, name, initial_masks[name], regions[name]
            )
            measurement_area = math.pi * (float(case["measurement_radius_um"]) * 1e-6) ** 2
            centroid_advance = (
                snapshot.area_weighted_position_m
                - initial_snapshots[name].area_weighted_position_m
            ) * 1e6
            global_allocated_advance = snapshot.cumulative_allocated_volume_m3 / growth_area * 1e6
            global_realized_advance = (
                snapshot.realized_voxels * system.voxel_volume_m3 / growth_area * 1e6
            )
            interior_deposited_advance = float(inventory["deposited_volume_m3"]) / measurement_area * 1e6
            interior_realized_advance = float(inventory["realized_volume_m3"]) / measurement_area * 1e6
            stage_rows.append({
                "case": case["case"],
                "voxel_size_um": case["voxel_size_um"],
                "growth_radius_um": case["growth_radius_um"],
                "measurement_radius_um": case["measurement_radius_um"],
                "orientation_family": family,
                "stage_fraction": stage_fraction,
                "expected_advance_um": expected_advance_um,
                "global_fractional_deposited_advance_um": global_allocated_advance,
                "global_binary_realized_advance_um": global_realized_advance,
                "current_area_weighted_advance_um": centroid_advance,
                "interior_fractional_deposited_advance_um": interior_deposited_advance,
                "interior_binary_realized_advance_um": interior_realized_advance,
                "interior_height_mean_advance_um": (float(height["mean_m"]) - float(initial_heights[name]["mean_m"])) * 1e6,
                "interior_height_median_advance_um": (float(height["median_m"]) - float(initial_heights[name]["median_m"])) * 1e6,
                "interior_height_std_um": float(height["std_m"]) * 1e6,
                "interior_height_iqr_um": float(height["iqr_m"]) * 1e6,
                "interior_binary_lag_um": float(inventory["binary_lag_volume_m3"]) / measurement_area * 1e6,
                "reconstructed_active_area_m2": snapshot.reconstructed_area_m2,
                "initial_reconstructed_active_area_m2": initial_snapshots[name].reconstructed_area_m2,
                "realized_voxels": snapshot.realized_voxels,
                "active_fractional_volume_m3": snapshot.active_fractional_volume_m3,
                "inactive_fractional_volume_m3": snapshot.inactive_fractional_volume_m3,
                "cumulative_allocated_volume_m3": snapshot.cumulative_allocated_volume_m3,
                "global_closure_error_m3": system.closure_error_m3(),
                "cumulative_deferred_volume_m3": system.cumulative_deferred_volume_m3,
            })
            audit = interface_filter_audit(system, name, initial_keys[name])
            filter_rows.append({
                "case": case["case"],
                "voxel_size_um": case["voxel_size_um"],
                "growth_radius_um": case["growth_radius_um"],
                "orientation_family": family,
                "stage_fraction": stage_fraction,
                **audit,
            })
            for alpha, controls in controls_by_alpha.items():
                control = controls[index]
                variant_rows.append({
                    "case": case["case"],
                    "voxel_size_um": case["voxel_size_um"],
                    "growth_radius_um": case["growth_radius_um"],
                    "orientation_family": family,
                    "stage_fraction": stage_fraction,
                    "alpha": alpha,
                    **{key: value for key, value in control.items() if key != "front"},
                    "variant_a_frozen_advance_um": stage_input_volume * float(control["variant_a_share"]) / growth_area * 1e6,
                    "variant_b_frozen_advance_um": stage_input_volume * float(control["variant_b_share"]) / growth_area * 1e6,
                    "variant_c_frozen_advance_um": stage_input_volume * float(control["variant_c_share"]) / growth_area * 1e6,
                    "actual_history_integrated_fractional_advance_um": global_allocated_advance,
                    "variant_d_binary_realized_advance_um": global_realized_advance,
                    "measured_centroid_advance_um": centroid_advance,
                })
    final_reservoir_rows = reservoir_rows(system, FRONT_NAMES)
    for row in final_reservoir_rows:
        row.update({
            "case": case["case"],
            "voxel_size_um": case["voxel_size_um"],
            "growth_radius_um": case["growth_radius_um"],
        })
    _, peak_memory = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    runtime = {
        "case": case["case"],
        "voxel_size_um": case["voxel_size_um"],
        "growth_radius_um": case["growth_radius_um"],
        "wall_runtime_s": perf_counter() - started,
        "dynamic_update_runtime_s": update_runtime_s,
        "peak_tracemalloc_bytes": peak_memory,
        "event_batches": sum(update.event_batches for update in system.history),
        "local_patch_updates": system.local_patch_update_count,
        "local_recomputed_elements": system.local_recomputed_elements,
        "closure_error_m3": system.closure_error_m3(),
        "deferred_volume_m3": system.cumulative_deferred_volume_m3,
    }
    return {
        "stage_rows": stage_rows,
        "variant_rows": variant_rows,
        "filter_rows": filter_rows,
        "reservoir_rows": final_reservoir_rows,
        "runtime": runtime,
    }


def error_attribution(stage_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in stage_rows:
        if float(record["stage_fraction"]) != 1.0:
            continue
        ideal = float(record["expected_advance_um"])
        fractional = float(record["global_fractional_deposited_advance_um"])
        binary = float(record["global_binary_realized_advance_um"])
        measured = float(record["current_area_weighted_advance_um"])
        normal = 0.0
        area_allocation = fractional - ideal
        quantization = binary - fractional
        measurement = measured - binary
        total = measured - ideal
        rows.append({
            "case": record["case"],
            "voxel_size_um": record["voxel_size_um"],
            "growth_radius_um": record["growth_radius_um"],
            "orientation_family": record["orientation_family"],
            "ideal_advance_um": ideal,
            "normal_component_um": normal,
            "evolving_area_allocation_component_um": area_allocation,
            "binary_quantization_component_um": quantization,
            "measurement_definition_component_um": measurement,
            "total_error_um": total,
            "component_sum_um": normal + area_allocation + quantization + measurement,
            "component_closure_um": total - (normal + area_allocation + quantization + measurement),
            "area_allocation_relative_to_ideal": area_allocation / ideal,
            "binary_quantization_relative_to_ideal": quantization / ideal,
            "measurement_definition_relative_to_ideal": measurement / ideal,
            "total_relative_error_signed": total / ideal,
            "interior_fractional_deposited_advance_um": record["interior_fractional_deposited_advance_um"],
            "interior_binary_realized_advance_um": record["interior_binary_realized_advance_um"],
            "interior_height_mean_advance_um": record["interior_height_mean_advance_um"],
            "interior_height_median_advance_um": record["interior_height_median_advance_um"],
        })
    return rows


def finite_window_comparison(stage_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    final = [row for row in stage_rows if float(row["stage_fraction"]) == 1.0]
    lookup = {
        (float(row["voxel_size_um"]), str(row["orientation_family"]), float(row["growth_radius_um"])): row
        for row in final
    }
    radii = sorted({float(row["growth_radius_um"]) for row in final})
    if len(radii) != 2:
        raise ValueError("finite-window comparison requires exactly two growth radii")
    rows: list[dict[str, object]] = []
    for voxel in sorted({float(row["voxel_size_um"]) for row in final}, reverse=True):
        for family in FAMILIES:
            small = lookup[(voxel, family, radii[0])]
            large = lookup[(voxel, family, radii[1])]
            rows.append({
                "voxel_size_um": voxel,
                "orientation_family": family,
                "small_growth_radius_um": radii[0],
                "large_growth_radius_um": radii[1],
                "expected_advance_um": small["expected_advance_um"],
                "small_interior_fractional_advance_um": small["interior_fractional_deposited_advance_um"],
                "large_interior_fractional_advance_um": large["interior_fractional_deposited_advance_um"],
                "large_minus_small_interior_fractional_um": float(large["interior_fractional_deposited_advance_um"]) - float(small["interior_fractional_deposited_advance_um"]),
                "small_interior_binary_advance_um": small["interior_binary_realized_advance_um"],
                "large_interior_binary_advance_um": large["interior_binary_realized_advance_um"],
                "large_minus_small_interior_binary_um": float(large["interior_binary_realized_advance_um"]) - float(small["interior_binary_realized_advance_um"]),
                "small_interior_height_mean_um": small["interior_height_mean_advance_um"],
                "large_interior_height_mean_um": large["interior_height_mean_advance_um"],
                "large_minus_small_interior_height_mean_um": float(large["interior_height_mean_advance_um"]) - float(small["interior_height_mean_advance_um"]),
                "small_full_window_centroid_um": small["current_area_weighted_advance_um"],
                "large_full_window_centroid_um": large["current_area_weighted_advance_um"],
            })
    return rows


def locked_error_attribution(path: Path) -> list[dict[str, object]]:
    """Reuse the immutable alpha=0 locked rows in the same additive ledger."""

    with path.open(newline="", encoding="utf-8") as handle:
        source = list(csv.DictReader(handle))
    rows: list[dict[str, object]] = []
    for record in source:
        if float(record["alpha"]) != 0.0 or float(record["voxel_size_um"]) not in (10.0, 5.0):
            continue
        ideal = float(record["expected_continuum_advance_um"])
        fractional = float(record["partial_volume_equivalent_advance_um"])
        binary = float(record["realized_volume_equivalent_advance_um"])
        measured = float(record["actual_reconstructed_front_advance_um"])
        area_allocation = fractional - ideal
        quantization = binary - fractional
        measurement = measured - binary
        total = measured - ideal
        rows.append({
            "voxel_size_um": float(record["voxel_size_um"]),
            "orientation_family": record["orientation_family"],
            "ideal_advance_um": ideal,
            "normal_component_um": 0.0,
            "evolving_area_allocation_component_um": area_allocation,
            "binary_quantization_component_um": quantization,
            "measurement_definition_component_um": measurement,
            "total_error_um": total,
            "component_sum_um": area_allocation + quantization + measurement,
            "component_closure_um": total - (area_allocation + quantization + measurement),
            "area_allocation_relative_to_ideal": area_allocation / ideal,
            "binary_quantization_relative_to_ideal": quantization / ideal,
            "measurement_definition_relative_to_ideal": measurement / ideal,
            "total_relative_error_signed": total / ideal,
            "initial_reconstructed_area_m2": float(record["initial_reconstructed_area_m2"]),
            "final_reconstructed_area_m2": float(record["final_reconstructed_area_m2"]),
            "active_fractional_volume_m3": float(record["active_fractional_volume_m3"]),
            "inactive_fractional_volume_m3": float(record["inactive_fractional_volume_m3"]),
            "locked_runtime_s": float(record["runtime_s"]),
        })
    return rows


def write_checksums(output: Path) -> None:
    lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (output / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    benchmark_started = perf_counter()
    benchmark_rows = translated_plane_benchmarks(config["translated_planes"])
    benchmark_runtime = perf_counter() - benchmark_started
    section = config["targeted_growth"]
    cases = []
    for index, (voxel, radius) in enumerate(
        (pair for voxel in section["voxel_sizes_um"] for pair in [(voxel, radius) for radius in section["growth_radii_um"]])
    ):
        cases.append({
            **section,
            "case": f"dx{float(voxel):g}_r{float(radius):g}",
            "voxel_size_um": voxel,
            "growth_radius_um": radius,
            "rng_seed": int(section["rng_seed"]) + index,
        })
    workers = min(int(config["resource_limits"]["maximum_workers"]), len(cases))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(growth_case_worker, cases))
    stage_rows = [row for result in results for row in result["stage_rows"]]
    variant_rows = [row for result in results for row in result["variant_rows"]]
    filter_rows = [row for result in results for row in result["filter_rows"]]
    all_reservoir_rows = [row for result in results for row in result["reservoir_rows"]]
    runtime_rows = [result["runtime"] for result in results]
    attribution_rows = error_attribution(stage_rows)
    finite_rows = finite_window_comparison(stage_rows)
    locked_source = REPOSITORY_ROOT / config["locked_results_preserved"] / "evolving_planar_fronts.csv"
    locked_rows = locked_error_attribution(locked_source)
    write_csv(output / "translated_plane_benchmarks.csv", benchmark_rows)
    write_csv(output / "targeted_growth_stages.csv", stage_rows)
    write_csv(output / "controlled_variants.csv", variant_rows)
    write_csv(output / "interface_filtering.csv", filter_rows)
    write_csv(output / "reservoir_fractions.csv", all_reservoir_rows)
    write_csv(output / "error_attribution.csv", attribution_rows)
    write_csv(output / "finite_window_comparison.csv", finite_rows)
    write_csv(output / "locked_error_attribution.csv", locked_rows)
    write_json(output / "runtime_resources.json", {
        "resource_policy": config["resource_limits"],
        "visible_logical_cpu_count": os.cpu_count(),
        "gpu_used": False,
        "memory_capacity": "unknown; peak Python allocations reported per worker",
        "translated_plane_runtime_s": benchmark_runtime,
        "growth_cases": runtime_rows,
    })
    write_json(output / "mass_balance.json", {
        "maximum_absolute_closure_error_m3": max(abs(float(row["global_closure_error_m3"])) for row in stage_rows),
        "maximum_deferred_volume_m3": max(float(row["cumulative_deferred_volume_m3"]) for row in stage_rows),
        "cases": [
            {
                "case": row["case"],
                "closure_error_m3": row["closure_error_m3"],
                "deferred_volume_m3": row["deferred_volume_m3"],
            }
            for row in runtime_rows
        ],
    })
    write_json(output / "manifest.json", {
        "phase": config["phase"],
        "base_commit": config["base_commit"],
        "config": str(config_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
        "config_sha256": sha256(config_path),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "locked_results_preserved": config["locked_results_preserved"],
        "locked_source_sha256": sha256(locked_source),
        "diagnostic_status": config["diagnostic_status"],
        "scope_exclusions": config["scope_exclusions"],
    })
    write_checksums(output)
    print(json.dumps({
        "translated_plane_rows": len(benchmark_rows),
        "growth_stage_rows": len(stage_rows),
        "attribution_rows": len(attribution_rows),
        "maximum_absolute_closure_error_m3": max(abs(float(row["global_closure_error_m3"])) for row in stage_rows),
        "output": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
