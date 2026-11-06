"""Run Phase 2A.3b-3 interface-area closure development or locked checks."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
from time import perf_counter
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import run_phase2a3b1_validation as phase2a3b1  # noqa: E402
from salt_phase2 import (  # noqa: E402
    AreaCandidateWeightModel,
    ConservativePatchLedger,
    Material,
    OrientationCandidateWeightModel,
)
from salt_phase2.ca import SaltGrowthCA  # noqa: E402
from salt_phase2.interface_flux import (  # noqa: E402
    allocate_interface_volume,
    bundle_elements_by_candidate,
    corrected_digital_interface,
    frozen_front_realization,
    triangulated_interface,
)
from salt_phase2.orientation_verification import (  # noqa: E402
    ReferenceGeometry,
    central_plane_box_area,
)
from salt_phase2.surface_orientation import (  # noqa: E402
    angular_errors_deg,
    cubic_attachment_weight,
    rotation_matrix_xyz_degrees,
)
from salt_phase2.anisotropy_diagnosis import support_advances_m  # noqa: E402


FAMILIES = ("<100>", "<110>", "<111>")
NORMALS_CRYSTAL = np.asarray(
    [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 1.0]], dtype=float
)
NORMALS_CRYSTAL /= np.linalg.norm(NORMALS_CRYSTAL, axis=1)[:, None]


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


def analytical_rows(config: dict[str, object]) -> list[dict[str, object]]:
    side = float(config["analytical_patch_side_um"]) * 1.0e-6
    areas = np.asarray(
        [side**2, math.sqrt(2.0) * side**2, 3.0 * math.sqrt(3.0) * side**2 / 4.0]
    )
    total_volume = float(config["analytical_total_volume_m3"])
    rotations = {
        "identity": np.eye(3),
        "joint_xyz_19_-31_7": rotation_matrix_xyz_degrees(19.0, -31.0, 7.0),
    }
    rows = []
    for rotation_name, rotation in rotations.items():
        normals_lab = NORMALS_CRYSTAL @ rotation.T
        for alpha in config["alpha_values"]:
            kappa = cubic_attachment_weight(normals_lab, rotation, float(alpha))
            allocation = allocate_interface_volume(total_volume, areas, kappa)
            reference = allocation.normal_advances_m[-1]
            for candidate_count in config["candidate_count_variants"]:
                for index, family in enumerate(FAMILIES):
                    rows.append({
                        "rotation": rotation_name,
                        "alpha": alpha,
                        "candidate_count_metadata_only": candidate_count,
                        "orientation_family": family,
                        "area_m2": areas[index],
                        "kappa": kappa[index],
                        "allocated_volume_m3": allocation.patch_volumes_m3[index],
                        "normal_advance_m": allocation.normal_advances_m[index],
                        "advance_ratio_to_111": allocation.normal_advances_m[index] / reference,
                        "kappa_ratio_to_111": kappa[index] / kappa[-1],
                        "ratio_relative_error": abs(
                            allocation.normal_advances_m[index] / reference - kappa[index] / kappa[-1]
                        ) / (kappa[index] / kappa[-1]),
                        "volume_closure_error_m3": allocation.closure_error_m3,
                    })
    return rows


def geometry_set(section: dict[str, object], *, locked: bool) -> list[ReferenceGeometry]:
    length = float(section["domain_length_um"]) * 1.0e-6
    plane = np.asarray(section["plane_normal"], dtype=float)
    plane /= np.linalg.norm(plane)
    centre = (17e-6, -11e-6, 9e-6) if locked else (0.0, 0.0, 0.0)
    return [
        ReferenceGeometry("plane", "plane", length, normal=tuple(plane)),
        ReferenceGeometry(
            "sphere", "sphere", length, centre_m=centre,
            radius_m=float(section["sphere_radius_um"]) * 1.0e-6,
        ),
        ReferenceGeometry(
            "rotated_cube", "cube", length,
            half_width_m=float(section["cube_half_width_um"]) * 1.0e-6,
            rotation_xyz_deg=tuple(float(v) for v in section["cube_rotation_xyz_deg"]),
        ),
    ]


def surface_rows(config: dict[str, object], *, locked: bool) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    section = config["surface_benchmarks"]
    rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    for geometry in geometry_set(section, locked=locked):
        for voxel_um in section["voxel_sizes_um"]:
            dx = float(voxel_um) * 1.0e-6
            mask = geometry.voxelize(dx)
            fitting = max(float(section["fitting_radius_um"]) * 1.0e-6, 1.01 * dx)
            methods = {
                "triangulated": triangulated_interface(mask, dx),
                "corrected_digital": corrected_digital_interface(mask, dx, fitting),
            }
            for method, elements in methods.items():
                reference, classes = geometry.reference_normals(
                    elements.centroids_m, max(2.0 * dx, fitting)
                )
                errors = angular_errors_deg(elements.normals, reference)
                unique = np.isin(classes, ["smooth", "face"]) & elements.normal_valid & np.isfinite(errors)
                bundles = bundle_elements_by_candidate(elements)
                rows.append({
                    "geometry": geometry.name,
                    "method": method,
                    "voxel_size_um": voxel_um,
                    "cells_per_domain": mask.shape[0],
                    "element_count": len(elements.areas_m2),
                    "candidate_bundle_count": len(bundles),
                    "multiface_bundle_fraction": float(np.mean([p.element_count > 1 for p in bundles])),
                    "valid_normal_fraction": float(np.mean(elements.normal_valid)),
                    "physical_edge_element_fraction": float(np.mean(np.isin(classes, ["edge", "corner"]))),
                    "estimated_area_m2": elements.total_area_m2,
                    "reference_area_m2": geometry.exact_area_m2,
                    "signed_area_relative_error": elements.total_area_m2 / geometry.exact_area_m2 - 1.0,
                    "absolute_area_relative_error": abs(elements.total_area_m2 / geometry.exact_area_m2 - 1.0),
                    "unique_normal_median_error_deg": float(np.median(errors[unique])) if np.any(unique) else float("nan"),
                    "unique_normal_p90_error_deg": float(np.percentile(errors[unique], 90)) if np.any(unique) else float("nan"),
                    "runtime_s": elements.runtime_s,
                })
                if geometry.kind == "cube" and method == "corrected_digital":
                    face_kappa = cubic_attachment_weight(
                        elements.normals[elements.normal_valid], geometry.rotation, 1.5 * math.log(2.0)
                    )
                    face_flux = float(np.sum(elements.areas_m2[elements.normal_valid] * face_kappa))
                    bundle_valid = [patch for patch in bundles if patch.normal_valid]
                    bundle_normals = np.asarray([patch.normal for patch in bundle_valid])
                    bundle_kappa = cubic_attachment_weight(bundle_normals, geometry.rotation, 1.5 * math.log(2.0))
                    bundle_flux = float(np.sum(
                        np.asarray([patch.area_m2 for patch in bundle_valid]) * bundle_kappa
                    ))
                    edge_rows.append({
                        "voxel_size_um": voxel_um,
                        "rule_a": "unique_face_partition_then_kappa",
                        "rule_b": "candidate_bundle_then_kappa",
                        "face_partition_weight_m2": face_flux,
                        "candidate_bundle_weight_m2": bundle_flux,
                        "relative_sensitivity": abs(bundle_flux / face_flux - 1.0),
                        "digital_multiface_fraction": float(np.mean([p.element_count > 1 for p in bundles])),
                        "physical_edge_element_fraction": float(np.mean(np.isin(classes, ["edge", "corner"]))),
                        "selection_basis": "predeclared geometric partition; neither rule selected by appearance",
                    })
    return rows, edge_rows


def select_surface_method(rows: Sequence[dict[str, object]], config: dict[str, object]) -> dict[str, object]:
    sizes = sorted({float(row["voxel_size_um"]) for row in rows}, reverse=True)
    coarsest, finest = sizes[0], sizes[-1]
    candidates = []
    for method in ("triangulated", "corrected_digital"):
        method_rows = [row for row in rows if row["method"] == method]
        per_geometry = {}
        monotone_endpoint = True
        for geometry in ("plane", "sphere", "rotated_cube"):
            coarse = next(row for row in method_rows if row["geometry"] == geometry and float(row["voxel_size_um"]) == coarsest)
            fine = next(row for row in method_rows if row["geometry"] == geometry and float(row["voxel_size_um"]) == finest)
            per_geometry[geometry] = float(fine["absolute_area_relative_error"])
            monotone_endpoint &= float(fine["absolute_area_relative_error"]) <= float(coarse["absolute_area_relative_error"]) + 1.0e-15
        candidates.append({
            "method": method,
            "endpoint_refinement_all_geometries": monotone_endpoint,
            "finest_errors": per_geometry,
            "maximum_finest_error": max(per_geometry.values()),
            "within_predeclared_error": max(per_geometry.values()) <= float(config["selection_rule"]["maximum_finest_grid_area_error"]),
        })
    eligible = [row for row in candidates if row["endpoint_refinement_all_geometries"] and row["within_predeclared_error"]]
    selected = min(eligible, key=lambda row: row["maximum_finest_error"])["method"] if eligible else None
    return {"candidates": candidates, "selected_method": selected, "selection_rule": config["selection_rule"]}


def central_plane_mask(cells: int, normal: np.ndarray) -> np.ndarray:
    axes = [np.arange(cells, dtype=float) - 0.5 * (cells - 1)] * 3
    x, y, z = np.meshgrid(*axes, indexing="ij", sparse=True)
    return np.asarray(normal[0] * x + normal[1] * y + normal[2] * z <= 0.0)


def locked_discrete_rows(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["frozen_front"]
    side = float(section["physical_patch_side_um"]) * 1.0e-6
    areas = np.asarray([
        side**2, math.sqrt(2.0) * side**2, 3.0 * math.sqrt(3.0) * side**2 / 4.0
    ])
    total_volume = math.fsum(areas) * float(section["isotropic_reference_advance_um"]) * 1.0e-6
    rows = []
    for alpha in config["alpha_values"]:
        kappa = cubic_attachment_weight(NORMALS_CRYSTAL, np.eye(3), float(alpha))
        exact = allocate_interface_volume(total_volume, areas, kappa)
        for voxel_um in section["voxel_sizes_um"]:
            started = perf_counter()
            dx = float(voxel_um) * 1.0e-6
            cells = int(round(side / dx))
            counts = np.asarray([
                len(phase2a3b1.candidate_coordinates_from_mask(central_plane_mask(cells, normal)))
                for normal in NORMALS_CRYSTAL
            ], dtype=int)
            new = frozen_front_realization(total_volume, areas, kappa, dx, counts)
            old_prob = counts * kappa / np.sum(counts * kappa)
            old_advance = total_volume * old_prob / areas
            for index, family in enumerate(FAMILIES):
                rows.append({
                    "alpha": alpha,
                    "voxel_size_um": voxel_um,
                    "orientation_family": family,
                    "candidate_count": counts[index],
                    "candidate_density_m_minus2": counts[index] / areas[index],
                    "exact_allocated_advance_um": exact.normal_advances_m[index] * 1.0e6,
                    "new_voxel_only_advance_um": new["voxel_only_advances_m"][index] * 1.0e6,
                    "new_relative_advance_error": abs(new["voxel_only_advances_m"][index] / exact.normal_advances_m[index] - 1.0),
                    "new_subvoxel_volume_m3": new["reservoir_volumes_m3"][index],
                    "new_retained_volume_m3": new["retained_volumes_m3"][index],
                    "old_site_weighted_advance_um": old_advance[index] * 1.0e6,
                    "old_site_weighted_relative_error": abs(old_advance[index] / exact.normal_advances_m[index] - 1.0),
                    "volume_closure_error_m3": new["closure_error_m3"],
                    "capacity_limited_patches": new["capacity_limited_patches"],
                    "runtime_s": perf_counter() - started,
                })
    return rows


def ledger_and_local_checks(config: dict[str, object]) -> tuple[dict[str, object], list[dict[str, object]]]:
    ledger = ConservativePatchLedger(2160.0, 10e-6)
    voxel_mass = ledger.voxel_mass_kg
    first = ledger.update(
        3.5 * voxel_mass, ("p0", "p1"), (1.0, 2.0), (1.0, 1.0),
        (("shared", "a", "b"), ("shared", "c", "d")),
    )
    before = copy.deepcopy(ledger.snapshot_state())
    rollback = False
    try:
        ledger.update(
            2.0 * voxel_mass, ("p0",), (1.0,), (1.0,), (("e", "f"),),
            inject_failure_after_voxels=1,
        )
    except RuntimeError:
        rollback = ledger.snapshot_state() == before
    ledger_row = {
        "first_realized_voxels": first.realized_voxels,
        "unique_consumed_candidates": len(ledger.consumed_candidates),
        "double_consumption_absent": len(ledger.consumed_candidates) == first.realized_voxels,
        "mass_closure_error_kg": ledger.closure_error_kg(),
        "rollback_exact": rollback,
    }

    section = config["local_consistency"]
    dx = float(section["voxel_size_um"]) * 1e-6
    geometry = ReferenceGeometry(
        "local_cube", "cube", float(section["domain_length_um"]) * 1e-6,
        half_width_m=float(section["cube_half_width_um"]) * 1e-6,
        rotation_xyz_deg=tuple(section["cube_rotation_xyz_deg"]),
    )
    mask = geometry.voxelize(dx)
    fitting = float(section["fitting_radius_um"]) * 1e-6
    before_bundles = {p.coordinate: p for p in bundle_elements_by_candidate(corrected_digital_interface(mask, dx, fitting))}
    insertion = min(before_bundles, key=lambda c: np.linalg.norm(np.asarray(c) - 0.5 * np.asarray(mask.shape)))
    updated = mask.copy()
    updated[insertion] = True
    full_after = {p.coordinate: p for p in bundle_elements_by_candidate(corrected_digital_interface(updated, dx, fitting))}
    radius_cells = fitting / dx + math.sqrt(3.0) + 1.0
    centre = np.asarray(insertion, dtype=float)
    affected = {
        coordinate for coordinate in set(before_bundles) | set(full_after)
        if np.linalg.norm(np.asarray(coordinate, dtype=float) - centre) <= radius_cells
    }
    local_after = {coordinate: patch for coordinate, patch in before_bundles.items() if coordinate not in affected}
    local_after.update({coordinate: patch for coordinate, patch in full_after.items() if coordinate in affected})
    keys_exact = set(local_after) == set(full_after)
    area_error = max(
        (abs(local_after[key].area_m2 - full_after[key].area_m2) for key in full_after), default=0.0
    ) if keys_exact else float("inf")
    normal_error = max((
        np.linalg.norm(np.asarray(local_after[key].normal) - np.asarray(full_after[key].normal))
        for key in full_after if local_after[key].normal_valid and full_after[key].normal_valid
    ), default=0.0) if keys_exact else float("inf")
    local_rows = [{
        "insertion_coordinate": json.dumps(insertion),
        "invalidation_radius_cells": radius_cells,
        "affected_candidate_count": len(affected),
        "candidate_keys_exact": keys_exact,
        "maximum_area_abs_error_m2": area_error,
        "maximum_normal_vector_error": normal_error,
        "note": "fresh full reconstruction was used as the validation oracle; only affected patch records were replaced",
    }]
    return ledger_row, local_rows


def single_grain_rows(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["single_grain"]
    voxel_m = float(section["voxel_size_um"]) * 1e-6
    grain_m = float(section["grain_size_um"]) * 1e-6
    rotation = rotation_matrix_xyz_degrees(*section["rotation_xyz_deg"])
    directions = np.asarray(section["directions_crystal"], dtype=float)
    rows = []
    for seed in section["rng_seeds"]:
        for method in ("old_site_weighted", "new_area_weighted"):
            state = phase2a3b1.continuous_cube_state(
                shape=tuple(section["shape"]), substrate_layers=1,
                grain_size_m=grain_m, voxel_size_m=voxel_m, rotation=rotation,
            )
            initial = state.material == Material.SALT
            model_class = OrientationCandidateWeightModel if method == "old_site_weighted" else AreaCandidateWeightModel
            model = model_class(
                voxel_size_m=voxel_m,
                fitting_radius_m=float(section["fitting_radius_um"]) * 1e-6,
                grain_to_lab=rotation,
                alpha=float(section["alpha"]),
                edge_policy="smooth_average",
                update_strategy="local",
                quality_minimum=float(section["quality_minimum"]),
            )
            ca = SaltGrowthCA(
                state, solid_density_kg_m3=2160.0, voxel_size_m=voxel_m,
                growth_mode="stochastic", update_mode="sequential",
                rng_seed=int(seed), candidate_weight_model=model,
            )
            started = perf_counter()
            result = ca.grow(np.nextafter(int(section["events"]) * ca.voxel_mass_kg, math.inf))
            final = state.material == Material.SALT
            advances = support_advances_m(
                initial, final, voxel_size_m=voxel_m,
                grain_to_lab=rotation, directions_crystal=directions,
            )
            row = {
                "method": method,
                "rng_seed": seed,
                "events": result.transaction.realized_solid_voxels,
                "committed_volume_m3": result.transaction.realized_solid_voxels * voxel_m**3,
                "mass_closure_error_kg": ca.allocator.closure_error_kg(),
                "runtime_s": perf_counter() - started,
            }
            for index, direction in enumerate(section["direction_labels"]):
                row[f"advance_{direction}_um"] = advances[index] * 1e6
            rows.append(row)
    return rows


def development(config: dict[str, object], output: Path) -> dict[str, object]:
    analytical = analytical_rows(config)
    surfaces, edges = surface_rows(config, locked=False)
    selection = select_surface_method(surfaces, config)
    max_analytical = max(float(row["ratio_relative_error"]) for row in analytical)
    gates = {
        "analytical_allocation_pass": max_analytical <= float(config["selection_rule"]["maximum_analytical_relative_error"]),
        "maximum_analytical_ratio_relative_error": max_analytical,
        "surface_selection": selection,
        "overall_pass": max_analytical <= float(config["selection_rule"]["maximum_analytical_relative_error"]) and selection["selected_method"] is not None,
    }
    gates["decision"] = "PROCEED_TO_LOCKED" if gates["overall_pass"] else "STOP_NO_GO"
    write_csv(output / "analytical_allocation.csv", analytical)
    write_csv(output / "surface_benchmarks.csv", surfaces)
    write_csv(output / "edge_rule_sensitivity.csv", edges)
    write_json(output / "method_selection.json", selection)
    write_json(output / "gate_summary.json", gates)
    return gates


def locked(config: dict[str, object], output: Path) -> dict[str, object]:
    analytical = analytical_rows(config)
    surfaces, edges = surface_rows(config, locked=True)
    discrete = locked_discrete_rows(config)
    ledger, local = ledger_and_local_checks(config)
    write_csv(output / "analytical_allocation.csv", analytical)
    write_csv(output / "surface_benchmarks.csv", surfaces)
    write_csv(output / "edge_rule_sensitivity.csv", edges)
    write_csv(output / "frozen_front_realization.csv", discrete)
    write_json(output / "mass_ledger_rollback.json", ledger)
    write_csv(output / "local_patch_consistency.csv", local)

    acceptance = config["acceptance"]
    selected = config["selected_surface_method"]
    selected_rows = [row for row in surfaces if row["method"] == selected]
    finest = min(float(row["voxel_size_um"]) for row in selected_rows)
    surface_finest = max(
        float(row["absolute_area_relative_error"])
        for row in selected_rows if float(row["voxel_size_um"]) == finest
    )
    alpha_grid = sorted({float(row["alpha"]) for row in discrete})
    sizes = sorted({float(row["voxel_size_um"]) for row in discrete}, reverse=True)
    refinement = {}
    refinement_slopes = {}
    for alpha in alpha_grid:
        maxima = [max(
            float(row["new_relative_advance_error"]) for row in discrete
            if float(row["alpha"]) == alpha and float(row["voxel_size_um"]) == size
        ) for size in sizes]
        refinement[str(alpha)] = maxima
        refinement_slopes[str(alpha)] = float(
            np.polyfit(np.log(np.asarray(sizes)), np.log(np.maximum(maxima, 1.0e-30)), 1)[0]
        )
    gates = {
        "analytical_ratio_pass": max(float(row["ratio_relative_error"]) for row in analytical) <= float(acceptance["analytical_ratio_relative_error_max"]),
        "mass_closure_pass": max(abs(float(row["volume_closure_error_m3"])) for row in discrete) <= float(acceptance["volume_closure_abs_m3_max"]),
        "surface_area_pass": surface_finest <= float(acceptance["finest_surface_area_relative_error_max"]),
        "surface_area_finest_max_relative_error": surface_finest,
        "discrete_refinement_by_alpha": refinement,
        "discrete_refinement_loglog_slopes": refinement_slopes,
        "discrete_refinement_pass": all(
            values[-1] < float(acceptance["finest_to_coarsest_error_ratio_max"]) * values[0]
            and refinement_slopes[key] >= float(acceptance["refinement_order_min"])
            for key, values in refinement.items()
        ),
        "discrete_finest_error_pass": max(values[-1] for values in refinement.values()) <= float(acceptance["finest_discrete_advance_relative_error_max"]),
        "capacity_pass": all(int(row["capacity_limited_patches"]) == 0 for row in discrete),
        "ledger_closure_pass": abs(float(ledger["mass_closure_error_kg"])) <= float(acceptance["mass_closure_abs_kg_max"]),
        "rollback_pass": bool(ledger["rollback_exact"]),
        "no_double_consumption_pass": bool(ledger["double_consumption_absent"]),
        "local_patch_consistency_pass": all(row["candidate_keys_exact"] and float(row["maximum_area_abs_error_m2"]) <= float(acceptance["local_area_abs_error_m2_max"]) and float(row["maximum_normal_vector_error"]) <= float(acceptance["local_normal_vector_error_max"]) for row in local),
    }
    pass_keys = [key for key in gates if key.endswith("_pass")]
    gates["overall_pass"] = all(bool(gates[key]) for key in pass_keys)
    if gates["overall_pass"]:
        single = single_grain_rows(config)
        write_csv(output / "single_grain_comparison.csv", single)
        gates["single_grain_executed"] = True
        gates["decision"] = "GO_PHASE2A3B3_NUMERICAL_CLOSURE_ONLY"
    else:
        gates["single_grain_executed"] = False
        gates["decision"] = "NO_GO_CONSERVATIVE_DISCRETE_REALIZATION"
    write_json(output / "gate_summary.json", gates)
    return gates


def manifest(config_path: Path, output: Path, gates: dict[str, object]) -> None:
    write_json(output / "manifest.json", {
        "phase": "Phase 2A.3b-3",
        "config": str(config_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
        "config_sha256": sha256(config_path),
        "base_commit": "58c8db851d8c5de91c6499e7b255ef37256bb530",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "cpu_count_visible": os.cpu_count(),
        "gpu_used": False,
        "decision": gates["decision"],
        "scope_exclusions": ["multigrain", "PNM", "film hydraulics", "permeability", "full vapour transport", "NaCl kinetic calibration"],
    })


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if config["stage"] == "method_development":
        gates = development(config, output)
    elif config["stage"] == "locked_verification":
        gates = locked(config, output)
    else:
        raise ValueError("unknown stage")
    manifest(config_path, output, gates)
    write_checksums(output)
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
