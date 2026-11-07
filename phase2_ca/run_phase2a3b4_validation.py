"""Run Phase 2A.3b-4 dynamic conservative interface-flux verification."""

from __future__ import annotations

import argparse
import copy
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import tracemalloc
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from salt_phase2.dynamic_interface_flux import (  # noqa: E402
    DynamicInterfaceFluxSystem,
    FrontRegion,
    candidate_face_flux_audit,
    perpendicular_face_counterexample,
)
from salt_phase2.interface_flux import corrected_digital_interface  # noqa: E402
from salt_phase2.orientation_verification import ReferenceGeometry  # noqa: E402


FAMILIES = ("<100>", "<110>", "<111>")
FRONT_NAMES = ("f100", "f110", "f111")
NORMALS = np.asarray([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
NORMALS /= np.linalg.norm(NORMALS, axis=1)[:, None]


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
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def write_checksums(output: Path) -> None:
    lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (output / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def central_plane_mask(cells: int, normal: np.ndarray) -> np.ndarray:
    axis = np.arange(cells, dtype=float) - 0.5 * (cells - 1)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij", sparse=True)
    return np.asarray(normal[0] * x + normal[1] * y + normal[2] * z <= 0.0)


def geometry_masks(section: dict[str, object]) -> dict[str, np.ndarray]:
    length = float(section["domain_length_um"]) * 1e-6
    dx = float(section["voxel_size_um"]) * 1e-6
    cells = int(round(length / dx))
    coordinates = (-0.5 * cells + 0.5 + np.arange(cells)) * dx
    x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij", sparse=True)
    plane = ReferenceGeometry("plane", "plane", length, normal=(1.0, 1.0, 1.0))
    sphere = ReferenceGeometry(
        "sphere", "sphere", length,
        centre_m=(7e-6, -5e-6, 3e-6),
        radius_m=float(section["sphere_radius_um"]) * 1e-6,
    )
    cube = ReferenceGeometry(
        "rotated_cube", "cube", length,
        half_width_m=float(section["cube_half_width_um"]) * 1e-6,
        rotation_xyz_deg=tuple(section["cube_rotation_xyz_deg"]),
    )
    return {
        "plane": plane.voxelize(dx),
        "sphere": sphere.voxelize(dx),
        "convex_edge": np.asarray((x <= 0.0) & (y <= 0.0) & np.ones_like(z, dtype=bool)),
        "concave_corner": np.asarray(~((x > 0.0) & (y > 0.0) & (z > 0.0))),
        "rotated_cube": cube.voxelize(dx),
    }


def adapter_audit(config: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    section = config["audit"]
    dx = float(section["voxel_size_um"]) * 1e-6
    fitting = float(section["fitting_radius_um"]) * 1e-6
    rows = []
    counter = []
    for alpha in config["alpha_values"]:
        counter.append(perpendicular_face_counterexample(float(alpha)))
        for geometry, mask in geometry_masks(section).items():
            elements = corrected_digital_interface(mask, dx, fitting)
            records = candidate_face_flux_audit(elements, np.eye(3), float(alpha))
            face_total = math.fsum(float(row["face_integrated_weight_m2"]) for row in records)
            adapter_total = math.fsum(float(row["adapter_weight_m2"]) for row in records)
            discrepancies = np.asarray([row["absolute_relative_discrepancy"] for row in records])
            multiface = np.asarray([int(row["element_count"]) > 1 for row in records])
            rows.append({
                "geometry": geometry,
                "alpha": alpha,
                "candidate_count": len(records),
                "multiface_candidate_count": int(np.count_nonzero(multiface)),
                "multiface_candidate_fraction": float(np.mean(multiface)),
                "face_integrated_total_weight_m2": face_total,
                "adapter_total_weight_m2": adapter_total,
                "total_signed_relative_discrepancy": adapter_total / face_total - 1.0,
                "candidate_mean_abs_relative_discrepancy": float(np.mean(discrepancies)),
                "candidate_p95_abs_relative_discrepancy": float(np.percentile(discrepancies, 95)),
                "candidate_max_abs_relative_discrepancy": float(np.max(discrepancies)),
                "multiface_max_abs_relative_discrepancy": (
                    float(np.max(discrepancies[multiface])) if np.any(multiface) else 0.0
                ),
            })
    return rows, counter


def build_planar_system(
    *,
    domain_length_um: float,
    front_radius_um: float,
    voxel_size_um: float,
    fitting_radius_um: float,
    minimum_alignment: float,
    density: float,
    alpha: float,
    rng_seed: int,
) -> DynamicInterfaceFluxSystem:
    dx = voxel_size_um * 1e-6
    cells = int(round(domain_length_um / voxel_size_um))
    if not math.isclose(cells * voxel_size_um, domain_length_um, abs_tol=1e-12):
        raise ValueError("domain length must be divisible by voxel size")
    fronts = [
        FrontRegion(
            name=name,
            solid_mask=central_plane_mask(cells, normal),
            reference_normal=normal,
            centre_m=np.zeros(3),
            radius_m=front_radius_um * 1e-6,
            minimum_alignment=minimum_alignment,
        )
        for name, normal in zip(FRONT_NAMES, NORMALS, strict=True)
    ]
    return DynamicInterfaceFluxSystem(
        fronts,
        voxel_size_m=dx,
        fitting_radius_m=fitting_radius_um * 1e-6,
        solid_density_kg_m3=density,
        grain_to_lab=np.eye(3),
        alpha=alpha,
        rng_seed=rng_seed,
    )


def dynamic_case_worker(case: dict[str, object]) -> dict[str, object]:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    system = build_planar_system(
        domain_length_um=float(case["domain_length_um"]),
        front_radius_um=float(case["front_radius_um"]),
        voxel_size_um=float(case["voxel_size_um"]),
        fitting_radius_um=float(case["fitting_radius_um"]),
        minimum_alignment=float(case["minimum_normal_alignment"]),
        density=float(case["solid_density_kg_m3"]),
        alpha=float(case["alpha"]),
        rng_seed=int(case["rng_seed"]),
    )
    initial = {name: system.front_snapshot(name) for name in FRONT_NAMES}
    radius = float(case["front_radius_um"]) * 1e-6
    total_volume = (
        len(FRONT_NAMES)
        * math.pi
        * radius**2
        * float(case["reference_advance_um"])
        * 1e-6
    )
    tracemalloc.start()
    update = system.advance(total_volume * float(case["solid_density_kg_m3"]))
    _, peak_memory = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    kappa = np.exp(float(case["alpha"]) * (np.asarray([1.0, 0.5, 1.0 / 3.0]) - 1.0 / 3.0))
    expected = (
        len(FRONT_NAMES)
        * float(case["reference_advance_um"])
        * kappa
        / np.sum(kappa)
    )
    rows = []
    consistency = {}
    for index, (name, family) in enumerate(zip(FRONT_NAMES, FAMILIES, strict=True)):
        final = system.front_snapshot(name)
        actual = (final.area_weighted_position_m - initial[name].area_weighted_position_m) * 1e6
        partial = final.cumulative_allocated_volume_m3 / initial[name].reconstructed_area_m2 * 1e6
        realized_equivalent = (
            final.realized_voxels * system.voxel_volume_m3
            / initial[name].reconstructed_area_m2 * 1e6
        )
        rows.append({
            "alpha": case["alpha"],
            "voxel_size_um": case["voxel_size_um"],
            "orientation_family": family,
            "expected_continuum_advance_um": expected[index],
            "actual_reconstructed_front_advance_um": actual,
            "partial_volume_equivalent_advance_um": partial,
            "realized_volume_equivalent_advance_um": realized_equivalent,
            "actual_relative_error": abs(actual / expected[index] - 1.0),
            "partial_relative_error": abs(partial / expected[index] - 1.0),
            "realized_relative_error": abs(realized_equivalent / expected[index] - 1.0),
            "initial_reconstructed_area_m2": initial[name].reconstructed_area_m2,
            "final_reconstructed_area_m2": final.reconstructed_area_m2,
            "initial_area_relative_error": initial[name].reconstructed_area_m2 / (math.pi * radius**2) - 1.0,
            "final_to_initial_area_ratio": final.reconstructed_area_m2 / initial[name].reconstructed_area_m2,
            "realized_voxels": final.realized_voxels,
            "active_fractional_volume_m3": final.active_fractional_volume_m3,
            "inactive_fractional_volume_m3": final.inactive_fractional_volume_m3,
            "cumulative_allocated_volume_m3": final.cumulative_allocated_volume_m3,
            "global_closure_error_m3": update.closure_error_m3,
            "event_batches": update.event_batches,
            "full_reconstruction_count": update.reconstruction_count,
            "local_patch_update_count": update.local_patch_update_count,
            "local_recomputed_elements": update.local_recomputed_elements,
            "runtime_s": update.runtime_s,
            "peak_tracemalloc_bytes": peak_memory,
        })
        consistency[name] = system.local_cache_consistency(name)
    return {"rows": rows, "consistency": consistency}


def reservoir_difference(first: DynamicInterfaceFluxSystem, second: DynamicInterfaceFluxSystem) -> float:
    owners = set(first.local_reservoirs_m3) | set(second.local_reservoirs_m3)
    return max((
        abs(first.local_reservoirs_m3.get(owner, 0.0) - second.local_reservoirs_m3.get(owner, 0.0))
        for owner in owners
    ), default=0.0)


def timestep_partition_check(section: dict[str, object], split_counts: Sequence[int]) -> list[dict[str, object]]:
    common = dict(
        domain_length_um=float(section["domain_length_um"]),
        front_radius_um=float(section["front_radius_um"]),
        voxel_size_um=float(section["voxel_size_um"]),
        fitting_radius_um=float(section["fitting_radius_um"]),
        minimum_alignment=float(section["minimum_normal_alignment"]),
        density=float(section["solid_density_kg_m3"]),
        alpha=0.0,
        rng_seed=int(section["rng_seed"]),
    )
    radius = common["front_radius_um"] * 1e-6
    total_volume = 3.0 * math.pi * radius**2 * float(section["reference_advance_um"]) * 1e-6
    total_mass = total_volume * common["density"]
    reference = build_planar_system(**common)
    reference.advance(total_mass)
    rows = []
    for split_count in split_counts:
        split = build_planar_system(**common)
        for _ in range(int(split_count)):
            split.advance(total_mass / int(split_count))
        masks_exact = all(
            np.array_equal(reference.fronts[name].solid_mask, split.fronts[name].solid_mask)
            for name in FRONT_NAMES
        )
        rows.append({
            "split_count": split_count,
            "masks_exact": masks_exact,
            "maximum_reservoir_abs_difference_m3": reservoir_difference(reference, split),
            "realized_voxel_difference": split.cumulative_realized_voxels - reference.cumulative_realized_voxels,
            "interface_allocated_abs_difference_m3": abs(
                split.cumulative_interface_allocated_volume_m3
                - reference.cumulative_interface_allocated_volume_m3
            ),
            "reference_closure_error_m3": reference.closure_error_m3(),
            "split_closure_error_m3": split.closure_error_m3(),
            "reference_event_batches": sum(update.event_batches for update in reference.history),
            "split_event_batches": sum(update.event_batches for update in split.history),
        })
    return rows


def rollback_check(section: dict[str, object]) -> dict[str, object]:
    system = build_planar_system(
        domain_length_um=float(section["domain_length_um"]),
        front_radius_um=float(section["front_radius_um"]),
        voxel_size_um=float(section["voxel_size_um"]),
        fitting_radius_um=float(section["fitting_radius_um"]),
        minimum_alignment=float(section["minimum_normal_alignment"]),
        density=float(section["solid_density_kg_m3"]),
        alpha=0.0,
        rng_seed=int(section["rng_seed"]),
    )
    radius = float(section["front_radius_um"]) * 1e-6
    volume = 3.0 * math.pi * radius**2 * float(section["reference_advance_um"]) * 1e-6
    before = copy.deepcopy(system.snapshot_state())
    caught = False
    try:
        system.advance(
            volume * float(section["solid_density_kg_m3"]),
            inject_failure_after_event_batches=1,
        )
    except RuntimeError as error:
        caught = "injected" in str(error)
    after = system.snapshot_state()
    geometry_exact = all(
        np.array_equal(before["fronts"][name].solid_mask, after["fronts"][name].solid_mask)
        for name in FRONT_NAMES
    )
    return {
        "failure_caught": caught,
        "geometry_exact": geometry_exact,
        "reservoirs_exact": before["local_reservoirs_m3"] == after["local_reservoirs_m3"],
        "input_ledger_exact": before["cumulative_input_volume_m3"] == after["cumulative_input_volume_m3"],
        "allocated_ledger_exact": before["cumulative_interface_allocated_volume_m3"] == after["cumulative_interface_allocated_volume_m3"],
        "deferred_ledger_exact": before["cumulative_deferred_volume_m3"] == after["cumulative_deferred_volume_m3"],
        "rng_exact": before["rng_state"] == after["rng_state"],
        "cache_keys_exact": before["surface_cache"].keys() == after["surface_cache"].keys(),
    }


def ownership_edge_cases(config: dict[str, object]) -> dict[str, object]:
    audit_rows, _ = adapter_audit(config)
    concave = [row for row in audit_rows if row["geometry"] == "concave_corner"]
    convex = [row for row in audit_rows if row["geometry"] == "convex_edge"]
    smoke = config["dynamic_smoke"] if "dynamic_smoke" in config else config["planar"]
    system = build_planar_system(
        domain_length_um=float(smoke["domain_length_um"]),
        front_radius_um=float(smoke["front_radius_um"]),
        voxel_size_um=float(
            smoke.get(
                "voxel_size_um",
                config.get("timestep_partition", {}).get("voxel_size_um"),
            )
        ),
        fitting_radius_um=float(smoke["fitting_radius_um"]),
        minimum_alignment=float(smoke["minimum_normal_alignment"]),
        density=float(smoke["solid_density_kg_m3"]),
        alpha=0.0,
        rng_seed=int(smoke["rng_seed"]),
    )
    active = system.active_owners()
    shape = next(iter(system.fronts.values())).solid_mask.shape
    domain_candidates_valid = all(
        all(0 <= coordinate[axis] < shape[axis] for axis in range(3))
        for _, coordinate in active
    )
    return {
        "concave_multiface_candidates_present": any(int(row["multiface_candidate_count"]) > 0 for row in concave),
        "concave_max_adapter_discrepancy": max(float(row["multiface_max_abs_relative_discrepancy"]) for row in concave),
        "convex_results_retained": len(convex),
        "candidate_collision_rule": "all distinct faces owned by one coordinate are integrated before one reservoir is updated",
        "patch_disappearance_rule": "reservoir ownership is coordinate-persistent; inactive deposited fractional volume is retained at that owner",
        "patch_creation_rule": "a newly active coordinate starts at zero unless it already owns a retained reservoir",
        "domain_boundary_rule": "only internal salt/empty faces create candidates; no exterior-domain candidate is synthesized",
        "domain_candidates_valid": domain_candidates_valid,
    }


def development(config: dict[str, object], output: Path) -> dict[str, object]:
    audit_rows, counter_rows = adapter_audit(config)
    section = config["dynamic_smoke"]
    case = {**section, "alpha": 0.0}
    smoke_result = dynamic_case_worker(case)
    timestep_rows = timestep_partition_check(section, [int(section["split_count"])])
    rollback = rollback_check(section)
    ownership = ownership_edge_cases(config)
    write_csv(output / "adapter_discrepancy.csv", audit_rows)
    write_csv(output / "perpendicular_face_counterexample.csv", counter_rows)
    write_csv(output / "dynamic_smoke.csv", smoke_result["rows"])
    write_csv(output / "timestep_partition.csv", timestep_rows)
    write_json(output / "rollback.json", rollback)
    write_json(output / "local_cache_consistency.json", smoke_result["consistency"])
    write_json(output / "ownership_edge_cases.json", ownership)
    acceptance = config["acceptance"]
    counter_error = max(abs(float(row["adapter_to_face_ratio"]) / float(row["exact_ratio_exp_minus_alpha_over_2"]) - 1.0) for row in counter_rows)
    concave_max = max(
        float(row["multiface_max_abs_relative_discrepancy"])
        for row in audit_rows if row["geometry"] == "concave_corner"
    )
    cache_rows = smoke_result["consistency"].values()
    gates = {
        "analytical_counterexample_pass": counter_error <= float(acceptance["counterexample_relative_error_max"]),
        "analytical_counterexample_max_relative_error": counter_error,
        "adapter_nonequivalence_detected": concave_max >= float(acceptance["concave_multiface_discrepancy_min"]),
        "concave_multiface_max_relative_discrepancy": concave_max,
        "smoke_mass_closure_pass": max(abs(float(row["global_closure_error_m3"])) for row in smoke_result["rows"]) <= float(acceptance["volume_closure_abs_m3_max"]),
        "timestep_partition_pass": all(row["masks_exact"] and abs(float(row["maximum_reservoir_abs_difference_m3"])) <= float(acceptance["timestep_reservoir_abs_m3_max"]) and int(row["realized_voxel_difference"]) == 0 for row in timestep_rows),
        "rollback_pass": all(bool(value) for value in rollback.values()),
        "local_cache_pass": all(bool(row["keys_exact"]) and bool(row["valid_flags_exact"]) and float(row["maximum_area_abs_error_m2"]) <= float(acceptance["local_area_abs_error_m2_max"]) and float(row["maximum_normal_vector_error"]) <= float(acceptance["local_normal_vector_error_max"]) for row in cache_rows),
        "ownership_edge_cases_pass": bool(ownership["concave_multiface_candidates_present"]) and bool(ownership["domain_candidates_valid"]),
    }
    gates["overall_pass"] = all(bool(value) for key, value in gates.items() if key.endswith("_pass") or key == "adapter_nonequivalence_detected")
    gates["decision"] = "PROCEED_TO_LOCKED" if gates["overall_pass"] else "STOP_NO_GO"
    write_json(output / "gate_summary.json", gates)
    return gates


def locked(config: dict[str, object], output: Path) -> dict[str, object]:
    section = config["planar"]
    cases = [
        {**section, "voxel_size_um": voxel, "alpha": alpha, "rng_seed": int(section["rng_seed"]) + index}
        for index, (voxel, alpha) in enumerate(
            (pair for voxel in section["voxel_sizes_um"] for pair in [(voxel, alpha) for alpha in config["alpha_values"]])
        )
    ]
    workers = min(int(config["resource_limits"]["maximum_workers"]), len(cases))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(dynamic_case_worker, cases))
    rows = [row for result in results for row in result["rows"]]
    cache_rows = []
    for case, result in zip(cases, results, strict=True):
        for name, record in result["consistency"].items():
            cache_rows.append({
                "alpha": case["alpha"],
                "voxel_size_um": case["voxel_size_um"],
                "front": name,
                **record,
            })
    write_csv(output / "evolving_planar_fronts.csv", rows)
    write_csv(output / "local_cache_consistency.csv", cache_rows)
    timestep_rows = timestep_partition_check(
        {**section, "voxel_size_um": config["timestep_partition"]["voxel_size_um"]},
        config["timestep_partition"]["split_counts"],
    )
    write_csv(output / "timestep_partition.csv", timestep_rows)
    rollback = rollback_check({**section, "voxel_size_um": config["timestep_partition"]["voxel_size_um"]})
    write_json(output / "rollback.json", rollback)
    ownership = ownership_edge_cases(config)
    write_json(output / "ownership_edge_cases.json", ownership)

    acceptance = config["acceptance"]
    alpha_zero = [row for row in rows if float(row["alpha"]) == 0.0]
    sizes = sorted({float(row["voxel_size_um"]) for row in rows}, reverse=True)
    finest = sizes[-1]
    alpha_zero_metrics = {}
    for size in sizes:
        selected = [row for row in alpha_zero if float(row["voxel_size_um"]) == size]
        actual = np.asarray([float(row["actual_reconstructed_front_advance_um"]) for row in selected])
        partial = np.asarray([float(row["partial_volume_equivalent_advance_um"]) for row in selected])
        alpha_zero_metrics[str(size)] = {
            "actual_orientation_spread_over_mean": float((np.max(actual) - np.min(actual)) / np.mean(actual)),
            "partial_orientation_spread_over_mean": float((np.max(partial) - np.min(partial)) / np.mean(partial)),
            "maximum_actual_relative_error": max(float(row["actual_relative_error"]) for row in selected),
            "maximum_partial_relative_error": max(float(row["partial_relative_error"]) for row in selected),
        }
    finest_rows = [row for row in rows if float(row["voxel_size_um"]) == finest]
    maximum_actual_finest = max(float(row["actual_relative_error"]) for row in finest_rows)
    maximum_partial_finest = max(float(row["partial_relative_error"]) for row in finest_rows)
    coarsest_actual = max(float(row["actual_relative_error"]) for row in rows if float(row["voxel_size_um"]) == sizes[0])
    gates = {
        "mass_closure_pass": max(abs(float(row["global_closure_error_m3"])) for row in rows) <= float(acceptance["volume_closure_abs_m3_max"]),
        "alpha_zero_actual_orientation_pass": alpha_zero_metrics[str(finest)]["actual_orientation_spread_over_mean"] <= float(acceptance["alpha_zero_actual_orientation_spread_max"]),
        "alpha_zero_partial_orientation_pass": alpha_zero_metrics[str(finest)]["partial_orientation_spread_over_mean"] <= float(acceptance["alpha_zero_partial_orientation_spread_max"]),
        "finest_actual_error_pass": maximum_actual_finest <= float(acceptance["finest_actual_relative_error_max"]),
        "finest_partial_error_pass": maximum_partial_finest <= float(acceptance["finest_partial_relative_error_max"]),
        "actual_refinement_pass": maximum_actual_finest < coarsest_actual,
        "timestep_partition_pass": all(row["masks_exact"] and int(row["realized_voxel_difference"]) == 0 and float(row["maximum_reservoir_abs_difference_m3"]) <= float(acceptance["timestep_reservoir_abs_m3_max"]) for row in timestep_rows),
        "rollback_pass": all(bool(value) for value in rollback.values()),
        "local_cache_pass": all(bool(row["keys_exact"]) and bool(row["valid_flags_exact"]) and float(row["maximum_area_abs_error_m2"]) <= float(acceptance["local_area_abs_error_m2_max"]) and float(row["maximum_normal_vector_error"]) <= float(acceptance["local_normal_vector_error_max"]) for row in cache_rows),
        "ownership_edge_cases_pass": bool(ownership["concave_multiface_candidates_present"]) and bool(ownership["domain_candidates_valid"]),
        "alpha_zero_metrics": alpha_zero_metrics,
        "maximum_finest_actual_relative_error": maximum_actual_finest,
        "maximum_finest_partial_relative_error": maximum_partial_finest,
    }
    pass_keys = [key for key in gates if key.endswith("_pass")]
    gates["planar_overall_pass"] = all(bool(gates[key]) for key in pass_keys)
    gates["single_grain_executed"] = False
    gates["decision"] = (
        "PLANAR_GO_SINGLE_GRAIN_REQUIRED"
        if gates["planar_overall_pass"]
        else "NO_GO_DYNAMIC_PLANAR_FRONT"
    )
    write_json(output / "gate_summary.json", gates)
    return gates


def write_manifest(config_path: Path, output: Path, gates: dict[str, object]) -> None:
    write_json(output / "manifest.json", {
        "phase": "Phase 2A.3b-4",
        "config": str(config_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
        "config_sha256": sha256(config_path),
        "base_commit": "920aed2b8933801eb44ac65efada513e8cbc243f",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "cpu_count_visible": os.cpu_count(),
        "gpu_used": False,
        "decision": gates["decision"],
        "scope_exclusions": ["multigrain", "PNM", "permeability", "film hydraulics", "vapour transport", "NaCl calibration"],
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    output = args.output.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output.mkdir(parents=True, exist_ok=True)
    if config["stage"] == "method_development":
        gates = development(config, output)
    elif config["stage"] == "locked_verification":
        gates = locked(config, output)
    else:
        raise ValueError("unknown stage")
    write_manifest(config_path, output, gates)
    write_checksums(output)
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
