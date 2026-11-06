"""Run Phase 2A.3b-1 mechanism development or locked verification."""

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
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy import ndimage


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase1 import EvaporationParameters, compute_mass_budget  # noqa: E402
from salt_phase2 import (  # noqa: E402
    Material,
    OrientationCandidateWeightModel,
    PhysicalState,
    SaltGrowthCA,
    air_exposed_face_area_m2,
    classify_crystal_directions,
    closure_tolerance_kg,
    crystallographic_candidate_diagnostics,
    directional_growth_extents,
    fixed_front_probabilities,
    map_candidate_normals,
    rotation_matrix_xyz_degrees,
    sample_fixed_front,
    symmetric_chamfer_distance,
)


FAMILY_NAMES = ("<100>", "<110>", "<111>", "other")
FAMILY_COLUMNS = {"<100>": "f100", "<110>": "f110", "<111>": "f111", "other": "other"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"no rows for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def normalize_svg(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    path.write_text(
        "\n".join(line.rstrip() for line in text.splitlines()) + "\n",
        encoding="utf-8",
    )


def write_checksums(output: Path) -> None:
    lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (output / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def rotation_from_angles(angles: Sequence[float]) -> np.ndarray:
    return rotation_matrix_xyz_degrees(*(float(value) for value in angles))


def continuous_cube_state(
    *,
    shape: tuple[int, int, int],
    substrate_layers: int,
    grain_size_m: float,
    voxel_size_m: float,
    rotation: np.ndarray,
    anchored: bool = False,
) -> PhysicalState:
    coordinates = [
        (-0.5 * count + 0.5 + np.arange(count)) * voxel_size_m
        for count in shape
    ]
    x, y, z = np.meshgrid(*coordinates, indexing="ij", sparse=True)
    local = [
        sum(component * rotation[j, i] for j, component in enumerate((x, y, z)))
        for i in range(3)
    ]
    centred = np.asarray(
        np.maximum.reduce([np.abs(component) for component in local])
        <= 0.5 * grain_size_m,
        dtype=bool,
    )
    seed = centred
    if anchored:
        occupied = np.argwhere(centred)
        shift = substrate_layers - int(np.min(occupied[:, 2]))
        shifted = occupied.copy()
        shifted[:, 2] += shift
        if np.min(shifted) < 0 or np.any(shifted >= np.asarray(shape)):
            raise RuntimeError("anchored seed shift leaves the domain")
        seed = np.zeros(shape, dtype=bool)
        seed[tuple(shifted.T)] = True
    material = np.full(shape, Material.EMPTY, dtype=np.uint8)
    material[:, :, :substrate_layers] = Material.SUBSTRATE
    if np.any(seed & (material == Material.SUBSTRATE)):
        raise RuntimeError("seed overlaps substrate")
    material[seed] = Material.SALT
    if anchored:
        contacts = np.argwhere(seed & (np.indices(shape)[2] == substrate_layers))
        centre_xy = 0.5 * (np.asarray(shape[:2]) - 1)
        distances = np.linalg.norm(contacts[:, :2] - centre_xy, axis=1)
        anchor_coordinate = tuple(int(value) for value in contacts[int(np.argmin(distances))])
        liquid = np.zeros(shape, dtype=bool)
        liquid[anchor_coordinate] = True
    else:
        liquid = seed.copy()
    state = PhysicalState(material, liquid, seed.copy(), substrate_layers)
    state.validate()
    return state


def sphere_mask(shape: tuple[int, int, int], radius_m: float, voxel_size_m: float) -> np.ndarray:
    axes = [(-0.5 * count + 0.5 + np.arange(count)) * voxel_size_m for count in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij", sparse=True)
    return np.asarray(x**2 + y**2 + z**2 <= radius_m**2, dtype=bool)


def candidate_coordinates_from_material(material: np.ndarray) -> list[tuple[int, int, int]]:
    salt = material == Material.SALT
    kernel = np.zeros((3, 3, 3), dtype=np.uint8)
    kernel[0, 1, 1] = kernel[2, 1, 1] = 1
    kernel[1, 0, 1] = kernel[1, 2, 1] = 1
    kernel[1, 1, 0] = kernel[1, 1, 2] = 1
    adjacent = ndimage.convolve(salt.astype(np.uint8), kernel, mode="constant") > 0
    eligible = adjacent & (material == Material.EMPTY)
    return [tuple(int(value) for value in row) for row in np.argwhere(eligible)]


def candidate_coordinates_from_mask(salt: np.ndarray) -> list[tuple[int, int, int]]:
    material = np.full(salt.shape, Material.EMPTY, dtype=np.uint8)
    material[salt] = Material.SALT
    return candidate_coordinates_from_material(material)


def fixed_front_geometry(config: dict[str, object], geometry_name: str) -> tuple[np.ndarray, np.ndarray]:
    fixed = config["fixed_front"]
    shape = tuple(int(value) for value in fixed["shape"])
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    if geometry_name == "sphere":
        return sphere_mask(shape, float(fixed["sphere_radius_um"]) * 1.0e-6, voxel_m), np.eye(3)
    if geometry_name == "rotated_cube":
        rotation = rotation_from_angles(fixed["cube_rotation_xyz_deg"])
        state = continuous_cube_state(
            shape=shape,
            substrate_layers=1,
            grain_size_m=float(config["grain_size_um"]) * 1.0e-6,
            voxel_size_m=voxel_m,
            rotation=rotation,
        )
        return state.material == Material.SALT, rotation
    raise ValueError(geometry_name)


def fixed_front_analysis(
    config: dict[str, object],
    *,
    sample_counts: Sequence[int],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fixed = config["fixed_front"]
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    fitting_m = float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"]) * 1.0e-6
    candidate_rows: list[dict[str, object]] = []
    group_rows: list[dict[str, object]] = []
    largest_sample = max(int(value) for value in sample_counts)
    distribution_index = 0
    for geometry_name in fixed["geometries"]:
        salt, rotation = fixed_front_geometry(config, geometry_name)
        coordinates = candidate_coordinates_from_mask(salt)
        for policy in config["edge_policies"]:
            for alpha in config["alpha_values"]:
                records = crystallographic_candidate_diagnostics(
                    salt,
                    coordinates,
                    voxel_size_m=voxel_m,
                    fitting_radius_m=fitting_m,
                    grain_to_lab=rotation,
                    alpha=float(alpha),
                    edge_policy=str(policy),
                    quality_minimum=float(config["quality_minimum"]),
                    near_tolerance_deg=float(config["near_family_tolerance_deg"]),
                )
                probabilities = fixed_front_probabilities(records)
                for sample_count in sample_counts:
                    counts = sample_fixed_front(
                        probabilities,
                        sample_count=int(sample_count),
                        rng_seed=int(fixed["rng_seed"]) + distribution_index * 101 + int(sample_count),
                    )
                    for family in FAMILY_NAMES:
                        selector = np.asarray(
                            [record.orientation_family == family for record in records]
                        )
                        expected = float(np.sum(probabilities[selector]))
                        observed_count = int(np.sum(counts[selector]))
                        observed = observed_count / int(sample_count)
                        standard_error = math.sqrt(
                            max(expected * (1.0 - expected) / int(sample_count), 0.0)
                        )
                        standardized = (
                            (observed - expected) / standard_error
                            if standard_error > 0.0
                            else 0.0
                        )
                        group_rows.append(
                            {
                                "geometry": geometry_name,
                                "edge_policy": policy,
                                "alpha": alpha,
                                "sample_count": sample_count,
                                "orientation_family": family,
                                "candidate_count": int(np.count_nonzero(selector)),
                                "represented_area_m2": float(
                                    np.sum(
                                        [
                                            record.represented_area_m2
                                            for record, selected in zip(records, selector, strict=True)
                                            if selected
                                        ]
                                    )
                                ),
                                "expected_probability": expected,
                                "empirical_probability": observed,
                                "absolute_error": abs(observed - expected),
                                "standard_error": standard_error,
                                "standardized_error": standardized,
                                "ci95_lower": max(observed - 1.96 * standard_error, 0.0),
                                "ci95_upper": min(observed + 1.96 * standard_error, 1.0),
                            }
                        )
                    if int(sample_count) == largest_sample:
                        for record, probability, count in zip(records, probabilities, counts, strict=True):
                            candidate_rows.append(
                                {
                                    "geometry": geometry_name,
                                    "edge_policy": policy,
                                    "alpha": alpha,
                                    "coordinate": json.dumps(record.coordinate),
                                    "site_class": record.site_class,
                                    "normal_lab": json.dumps(record.normal_lab),
                                    "normal_crystal": json.dumps(record.normal_crystal),
                                    "q4": record.q4,
                                    "orientation_family": record.orientation_family,
                                    "family_angular_error_deg": record.family_angular_error_deg,
                                    "represented_area_m2": record.represented_area_m2,
                                    "normal_quality": record.normal_quality,
                                    "anisotropy_eligible": record.anisotropy_eligible,
                                    "weight": record.anisotropy_weight,
                                    "expected_probability": probability,
                                    "empirical_count": int(count),
                                    "empirical_probability": int(count) / largest_sample,
                                }
                            )
                distribution_index += 1
    return candidate_rows, group_rows


def run_ca_once(
    *,
    shape: tuple[int, int, int],
    grain_size_m: float,
    voxel_size_m: float,
    fitting_radius_m: float,
    rotation: np.ndarray,
    alpha: float,
    edge_policy: str,
    strategy: str,
    quality_minimum: float,
    rng_seed: int,
    events: int,
) -> tuple[np.ndarray, tuple[tuple[int, int, int], ...], OrientationCandidateWeightModel, SaltGrowthCA]:
    state = continuous_cube_state(
        shape=shape,
        substrate_layers=1,
        grain_size_m=grain_size_m,
        voxel_size_m=voxel_size_m,
        rotation=rotation,
    )
    model = OrientationCandidateWeightModel(
        voxel_size_m=voxel_size_m,
        fitting_radius_m=fitting_radius_m,
        grain_to_lab=rotation,
        alpha=alpha,
        edge_policy=edge_policy,
        update_strategy=strategy,
        quality_minimum=quality_minimum,
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=2160.0,
        voxel_size_m=voxel_size_m,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=rng_seed,
        candidate_weight_model=model,
    )
    result = ca.grow(np.nextafter(events * ca.voxel_mass_kg, math.inf))
    if result.transaction.realized_solid_voxels != events:
        raise RuntimeError("cache case did not realize requested events")
    return (
        state.material == Material.SALT,
        result.selected_coordinates,
        model,
        ca,
    )


def cache_case_worker(case: dict[str, object]) -> dict[str, object]:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    rotation = rotation_from_angles(case["rotation_xyz_deg"])
    common = {
        "shape": tuple(case["shape"]),
        "grain_size_m": float(case["grain_size_um"]) * 1.0e-6,
        "voxel_size_m": float(case["voxel_size_um"]) * 1.0e-6,
        "fitting_radius_m": float(case["fitting_radius_um"]) * 1.0e-6,
        "rotation": rotation,
        "alpha": float(case["alpha"]),
        "edge_policy": str(case["edge_policy"]),
        "quality_minimum": float(case["quality_minimum"]),
        "rng_seed": int(case["rng_seed"]),
        "events": int(case["events"]),
    }
    started = perf_counter()
    full_mask, full_sequence, full_model, full_ca = run_ca_once(strategy="full", **common)
    full_runtime = perf_counter() - started
    started = perf_counter()
    local_mask, local_sequence, local_model, local_ca = run_ca_once(strategy="local", **common)
    local_runtime = perf_counter() - started
    weight_error = float("nan")
    if full_sequence == local_sequence:
        weight_error = float(
            np.max(
                np.abs(
                    np.asarray([record.anisotropy_weight for record in full_model.selection_history])
                    - np.asarray([record.anisotropy_weight for record in local_model.selection_history])
                )
            )
        )
    return {
        "rotation_index": case["rotation_index"],
        "rotation_xyz_deg": json.dumps(case["rotation_xyz_deg"]),
        "edge_policy": case["edge_policy"],
        "rng_seed": case["rng_seed"],
        "events": case["events"],
        "selected_sequence_exact": full_sequence == local_sequence,
        "final_mask_exact": bool(np.array_equal(full_mask, local_mask)),
        "selected_weight_max_abs_error": weight_error,
        "full_runtime_s": full_runtime,
        "local_runtime_s": local_runtime,
        "speedup_full_over_local": full_runtime / local_runtime,
        "full_recomputed_candidates": full_model.recomputed_candidates,
        "local_recomputed_candidates": local_model.recomputed_candidates,
        "full_mass_closure_error_kg": full_ca.allocator.closure_error_kg(),
        "local_mass_closure_error_kg": local_ca.allocator.closure_error_kg(),
    }


def cache_case_matrix(config: dict[str, object], section_name: str) -> list[dict[str, object]]:
    section = config[section_name]
    fitting_um = float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"])
    cases = []
    for rotation_index, angles in enumerate(section["rotations_xyz_deg"]):
        for policy in config["edge_policies"]:
            for seed in section["rng_seeds"]:
                cases.append(
                    {
                        "rotation_index": rotation_index,
                        "rotation_xyz_deg": angles,
                        "edge_policy": policy,
                        "rng_seed": seed,
                        "events": section["events"],
                        "alpha": section["alpha"],
                        "shape": section["shape"],
                        "grain_size_um": config["grain_size_um"],
                        "voxel_size_um": config["voxel_size_um"],
                        "fitting_radius_um": fitting_um,
                        "quality_minimum": config["quality_minimum"],
                    }
                )
    with ProcessPoolExecutor(max_workers=min(int(section["maximum_workers"]), len(cases))) as executor:
        return list(executor.map(cache_case_worker, cases))


def targeted_cache_scenarios(config: dict[str, object]) -> list[dict[str, object]]:
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    fitting_m = float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"]) * 1.0e-6
    shape = (24, 24, 24)
    scenarios: dict[str, np.ndarray] = {}
    convex = np.zeros(shape, dtype=bool)
    convex[7:17, 7:17, 7:17] = True
    scenarios["convex_corner"] = convex
    concave = np.zeros(shape, dtype=bool)
    concave[5:15, 5:19, 5:19] = True
    concave[5:19, 5:15, 5:19] = True
    scenarios["concave_reentrant"] = concave
    three_face = np.zeros(shape, dtype=bool)
    three_face[:12, :, :] = True
    three_face[:, :12, :] = True
    three_face[:, :, :12] = True
    scenarios["three_face_corner"] = three_face
    boundary = np.zeros(shape, dtype=bool)
    boundary[0:10, 7:17, 7:17] = True
    scenarios["domain_boundary"] = boundary
    rows = []
    for name, salt in scenarios.items():
        coordinates_before = candidate_coordinates_from_mask(salt)
        initial_records = map_candidate_normals(
            salt,
            coordinates_before,
            voxel_size_m=voxel_m,
            fitting_radius_m=fitting_m,
            grain_to_lab=np.eye(3),
            alpha=float(config["alpha_values"][0]),
            edge_policy="smooth_average",
            quality_minimum=float(config["quality_minimum"]),
        )
        preferred_classes = ("corner_like", "edge_like", "face_like")
        insertion_record = next(
            record
            for site_class in preferred_classes
            for record in initial_records
            if record.site_class == site_class
        )
        model = OrientationCandidateWeightModel(
            voxel_size_m=voxel_m,
            fitting_radius_m=fitting_m,
            grain_to_lab=np.eye(3),
            alpha=float(config["alpha_values"][0]),
            edge_policy="smooth_average",
            update_strategy="local",
            quality_minimum=float(config["quality_minimum"]),
        )
        model.diagnostics_for_candidates(salt, coordinates_before)
        updated = salt.copy()
        updated[insertion_record.coordinate] = True
        model.notify_insertion(insertion_record.coordinate)
        coordinates_after = candidate_coordinates_from_mask(updated)
        local = model.diagnostics_for_candidates(updated, coordinates_after)
        full = map_candidate_normals(
            updated,
            coordinates_after,
            voxel_size_m=voxel_m,
            fitting_radius_m=fitting_m,
            grain_to_lab=np.eye(3),
            alpha=float(config["alpha_values"][0]),
            edge_policy="smooth_average",
            quality_minimum=float(config["quality_minimum"]),
        )
        weight_error = float(
            np.max(
                np.abs(
                    np.asarray([record.anisotropy_weight for record in local])
                    - np.asarray([record.anisotropy_weight for record in full])
                )
            )
        )
        class_mismatch = sum(
            first.site_class != second.site_class
            for first, second in zip(local, full, strict=True)
        )
        rows.append(
            {
                "scenario": name,
                "insertion_coordinate": json.dumps(insertion_record.coordinate),
                "insertion_site_class": insertion_record.site_class,
                "candidate_count_after": len(coordinates_after),
                "maximum_weight_error": weight_error,
                "site_class_mismatch_count": class_mismatch,
                "pass": weight_error <= 1.0e-12 and class_mismatch == 0,
            }
        )
    return rows


def development_stage(config: dict[str, object], output: Path) -> dict[str, object]:
    candidate_rows, group_rows = fixed_front_analysis(
        config, sample_counts=config["fixed_front"]["sample_counts"]
    )
    write_csv(output / "fixed_front_candidates.csv", candidate_rows)
    write_csv(output / "fixed_front_probability_groups.csv", group_rows)
    cache_rows = cache_case_matrix(config, "cache_development")
    write_csv(output / "cache_full_vs_local.csv", cache_rows)
    targeted_rows = targeted_cache_scenarios(config)
    write_csv(output / "cache_targeted_insertions.csv", targeted_rows)
    minimum_samples = int(config["fixed_front"]["minimum_locked_sample_count"])
    locked_rows = [row for row in group_rows if int(row["sample_count"]) >= minimum_samples]
    probability_pass = all(
        float(row["absolute_error"])
        <= float(config["fixed_front"]["maximum_group_absolute_probability_error"])
        and abs(float(row["standardized_error"]))
        <= float(config["fixed_front"]["maximum_group_standardized_error"])
        for row in locked_rows
    )
    cache_pass = all(
        row["selected_sequence_exact"]
        and row["final_mask_exact"]
        and float(row["selected_weight_max_abs_error"])
        <= float(config["cache_development"]["maximum_selected_weight_error"])
        for row in cache_rows
    )
    targeted_pass = all(row["pass"] for row in targeted_rows)
    gates = {
        "fixed_front_probability_agreement": probability_pass,
        "cache_100_event_exactness": cache_pass,
        "targeted_cache_insertions": targeted_pass,
        "selected_locked_fixed_front_sample_count": minimum_samples,
        "overall_pass": probability_pass and cache_pass and targeted_pass,
        "decision": "PROCEED_TO_LOCKED" if probability_pass and cache_pass and targeted_pass else "STOP",
    }
    write_json(output / "gate_summary.json", gates)
    return gates


def cube_attribution(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["cube_attribution"]
    grain_m = float(config["grain_size_um"]) * 1.0e-6
    rotation = rotation_from_angles(section["rotation_xyz_deg"])
    rows = []
    for cells in section["cells_across"]:
        voxel_m = grain_m / int(cells)
        fitting_m = float(config["fitting_radius_over_grain_size"]) * grain_m
        shape = (3 * int(cells),) * 3
        state = continuous_cube_state(
            shape=shape,
            substrate_layers=1,
            grain_size_m=grain_m,
            voxel_size_m=voxel_m,
            rotation=rotation,
        )
        salt = state.material == Material.SALT
        coordinates = candidate_coordinates_from_material(state.material)
        baseline = crystallographic_candidate_diagnostics(
            salt,
            coordinates,
            voxel_size_m=voxel_m,
            fitting_radius_m=fitting_m,
            grain_to_lab=rotation,
            alpha=0.0,
            edge_policy="smooth_average",
            quality_minimum=float(config["quality_minimum"]),
            near_tolerance_deg=float(config["near_family_tolerance_deg"]),
        )
        p0 = fixed_front_probabilities(baseline)
        for policy in config["edge_policies"]:
            for alpha in config["alpha_values"]:
                records = crystallographic_candidate_diagnostics(
                    salt,
                    coordinates,
                    voxel_size_m=voxel_m,
                    fitting_radius_m=fitting_m,
                    grain_to_lab=rotation,
                    alpha=float(alpha),
                    edge_policy=str(policy),
                    quality_minimum=float(config["quality_minimum"]),
                    near_tolerance_deg=float(config["near_family_tolerance_deg"]),
                )
                probabilities = fixed_front_probabilities(records)
                total_variation = 0.5 * float(np.sum(np.abs(probabilities - p0)))
                for site_class in ("face_like", "edge_like", "corner_like"):
                    selector = np.asarray([record.site_class == site_class for record in records])
                    class_shift = 0.5 * float(np.sum(np.abs(probabilities[selector] - p0[selector])))
                    weights = np.asarray(
                        [record.anisotropy_weight for record, chosen in zip(records, selector, strict=True) if chosen]
                    )
                    rows.append(
                        {
                            "cells_across_grain": cells,
                            "voxel_size_um": voxel_m * 1.0e6,
                            "edge_policy": policy,
                            "alpha": alpha,
                            "site_class": site_class,
                            "candidate_count": int(np.count_nonzero(selector)),
                            "candidate_fraction": float(np.mean(selector)),
                            "represented_area_m2": float(
                                sum(record.represented_area_m2 for record, chosen in zip(records, selector, strict=True) if chosen)
                            ),
                            "mean_weight": float(np.mean(weights)) if weights.size else float("nan"),
                            "weight_coefficient_of_variation": (
                                float(np.std(weights) / np.mean(weights)) if weights.size else float("nan")
                            ),
                            "expected_probability": float(np.sum(probabilities[selector])),
                            "total_variation_contribution": class_shift,
                            "fraction_of_total_anisotropy_shift": (
                                class_shift / total_variation if total_variation > 0.0 else 0.0
                            ),
                        }
                    )
    return rows


def full_cache_check(
    model: OrientationCandidateWeightModel,
    state: PhysicalState,
    coordinates: list[tuple[int, int, int]],
) -> dict[str, object]:
    salt = state.material == Material.SALT
    local = model.diagnostics_for_candidates(salt, coordinates)
    full = map_candidate_normals(
        salt,
        coordinates,
        voxel_size_m=model.voxel_size_m,
        fitting_radius_m=model.fitting_radius_m,
        grain_to_lab=model.grain_to_lab,
        alpha=model.alpha,
        edge_policy=model.edge_policy,
        quality_minimum=model.quality_minimum,
    )
    local_weights = np.asarray([record.anisotropy_weight for record in local])
    full_weights = np.asarray([record.anisotropy_weight for record in full])
    angles = []
    for first, second in zip(local, full, strict=True):
        first_normal = np.asarray(first.combined_normal, dtype=float)
        second_normal = np.asarray(second.combined_normal, dtype=float)
        first_valid = np.all(np.isfinite(first_normal))
        second_valid = np.all(np.isfinite(second_normal))
        if not first_valid and not second_valid:
            angles.append(0.0)
        elif first_valid != second_valid:
            angles.append(float("inf"))
        else:
            dot = float(np.dot(first_normal, second_normal))
            angles.append(math.degrees(math.acos(np.clip(dot, -1.0, 1.0))))
    return {
        "maximum_weight_error": float(np.max(np.abs(local_weights - full_weights))),
        "maximum_normal_angle_error_deg": float(np.max(angles)),
        "site_class_mismatch_count": sum(
            first.site_class != second.site_class for first, second in zip(local, full, strict=True)
        ),
        "anisotropy_eligibility_mismatch_count": sum(
            first.anisotropy_eligible != second.anisotropy_eligible
            for first, second in zip(local, full, strict=True)
        ),
    }


def family_surface_summary(records) -> dict[str, float]:
    totals = {
        family: float(sum(record.represented_area_m2 for record in records if record.orientation_family == family))
        for family in FAMILY_NAMES
    }
    total = sum(totals.values())
    result = {f"surface_area_{FAMILY_COLUMNS[family]}_m2": totals[family] for family in FAMILY_NAMES}
    result.update(
        {
            f"surface_area_fraction_{FAMILY_COLUMNS[family]}": totals[family] / total if total else 0.0
            for family in FAMILY_NAMES
        }
    )
    result["surface_area_total_m2"] = total
    return result


def selected_event_rows(
    *,
    case_id: str,
    stage_fraction: float,
    stage_event_start: int,
    records,
    rotation: np.ndarray,
    shape: tuple[int, int, int],
    voxel_m: float,
    near_tolerance_deg: float,
) -> list[dict[str, object]]:
    lab_normals = np.asarray([record.combined_normal for record in records], dtype=float)
    crystal_normals = lab_normals @ rotation
    families, errors = classify_crystal_directions(
        crystal_normals, near_tolerance_deg=near_tolerance_deg
    )
    q4 = np.sum(crystal_normals**4, axis=1)
    centre = 0.5 * (np.asarray(shape, dtype=float) - 1.0)
    rows = []
    for offset, (record, normal_crystal, family, error, q4_value) in enumerate(
        zip(records, crystal_normals, families, errors, q4, strict=True)
    ):
        point_lab = (np.asarray(record.coordinate, dtype=float) - centre) * voxel_m
        point_crystal = point_lab @ rotation
        rows.append(
            {
                "case_id": case_id,
                "stage_fraction": stage_fraction,
                "event_index": stage_event_start + offset,
                "coordinate": json.dumps(record.coordinate),
                "site_class": record.site_class,
                "anisotropy_eligible": record.anisotropy_eligible,
                "normal_quality": record.normal_quality,
                "weight": record.anisotropy_weight,
                "normal_crystal": json.dumps(tuple(float(value) for value in normal_crystal)),
                "q4": float(q4_value),
                "orientation_family": str(family),
                "family_angular_error_deg": float(error),
                "represented_area_m2": record.adjacent_salt_faces * voxel_m**2,
                "point_crystal_x_m": point_crystal[0],
                "point_crystal_y_m": point_crystal[1],
                "point_crystal_z_m": point_crystal[2],
            }
        )
    return rows


def locked_growth_worker(case: dict[str, object]) -> dict[str, object]:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    shape = tuple(int(value) for value in case["shape"])
    voxel_m = float(case["voxel_size_um"]) * 1.0e-6
    grain_m = float(case["grain_size_um"]) * 1.0e-6
    rotation = rotation_from_angles(case["rotation_xyz_deg"])
    state = continuous_cube_state(
        shape=shape,
        substrate_layers=1,
        grain_size_m=grain_m,
        voxel_size_m=voxel_m,
        rotation=rotation,
    )
    initial_mask = state.seed_mask.copy()
    model = OrientationCandidateWeightModel(
        voxel_size_m=voxel_m,
        fitting_radius_m=float(case["fitting_radius_um"]) * 1.0e-6,
        grain_to_lab=rotation,
        alpha=float(case["alpha"]),
        edge_policy=str(case["edge_policy"]),
        update_strategy="local",
        quality_minimum=float(case["quality_minimum"]),
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=float(case["solid_density_kg_m3"]),
        voxel_size_m=voxel_m,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=int(case["rng_seed"]),
        candidate_weight_model=model,
    )
    previous_target = 0
    previous_surface = None
    event_rows: list[dict[str, object]] = []
    stage_rows: list[dict[str, object]] = []
    started = perf_counter()
    for stage_fraction, target_events in zip(case["stage_fractions"], case["stage_events"], strict=True):
        delta = int(target_events) - previous_target
        history_start = len(model.selection_history)
        result = ca.grow(np.nextafter(delta * ca.voxel_mass_kg, math.inf))
        if result.transaction.realized_solid_voxels != delta:
            raise RuntimeError("locked growth did not realize a stage target")
        selected = model.selection_history[history_start:]
        event_rows.extend(
            selected_event_rows(
                case_id=str(case["case_id"]),
                stage_fraction=float(stage_fraction),
                stage_event_start=previous_target,
                records=selected,
                rotation=rotation,
                shape=shape,
                voxel_m=voxel_m,
                near_tolerance_deg=float(case["near_family_tolerance_deg"]),
            )
        )
        coordinates = candidate_coordinates_from_material(state.material)
        cache_check = full_cache_check(model, state, coordinates)
        surface_records = crystallographic_candidate_diagnostics(
            state.material == Material.SALT,
            coordinates,
            voxel_size_m=voxel_m,
            fitting_radius_m=float(case["fitting_radius_um"]) * 1.0e-6,
            grain_to_lab=rotation,
            alpha=float(case["alpha"]),
            edge_policy=str(case["edge_policy"]),
            quality_minimum=float(case["quality_minimum"]),
            near_tolerance_deg=float(case["near_family_tolerance_deg"]),
        )
        surface = family_surface_summary(surface_records)
        if previous_surface is None:
            initial_state = continuous_cube_state(
                shape=shape,
                substrate_layers=1,
                grain_size_m=grain_m,
                voxel_size_m=voxel_m,
                rotation=rotation,
            )
            initial_coordinates = candidate_coordinates_from_material(initial_state.material)
            initial_records = crystallographic_candidate_diagnostics(
                initial_state.material == Material.SALT,
                initial_coordinates,
                voxel_size_m=voxel_m,
                fitting_radius_m=float(case["fitting_radius_um"]) * 1.0e-6,
                grain_to_lab=rotation,
                alpha=float(case["alpha"]),
                edge_policy=str(case["edge_policy"]),
                quality_minimum=float(case["quality_minimum"]),
                near_tolerance_deg=float(case["near_family_tolerance_deg"]),
            )
            previous_surface = family_surface_summary(initial_records)
        stage_event_subset = event_rows[-delta:]
        family_counts = {
            family: sum(row["orientation_family"] == family for row in stage_event_subset)
            for family in FAMILY_NAMES
        }
        final_mask = state.material == Material.SALT
        air_area = air_exposed_face_area_m2(
            final_mask,
            state.material == Material.EMPTY,
            voxel_m,
        )
        extents = directional_growth_extents(
            initial_mask,
            final_mask,
            voxel_size_m=voxel_m,
            grain_to_lab=rotation,
        )
        row: dict[str, object] = {
            "case_id": case["case_id"],
            "rotation_index": case["rotation_index"],
            "rotation_xyz_deg": json.dumps(case["rotation_xyz_deg"]),
            "edge_policy": case["edge_policy"],
            "alpha": case["alpha"],
            "rng_seed": case["rng_seed"],
            "stage_fraction": stage_fraction,
            "target_new_voxels": target_events,
            "realized_fraction_of_initial_voxels": int(target_events) / int(np.count_nonzero(initial_mask)),
            "initial_seed_voxels": int(np.count_nonzero(initial_mask)),
            "current_salt_voxels": int(np.count_nonzero(final_mask)),
            "cumulative_input_mass_kg": ca.allocator.cumulative_input_mass_kg,
            "cumulative_committed_mass_kg": ca.allocator.cumulative_committed_solid_mass_kg,
            "deferred_mass_kg": ca.allocator.deferred_precipitation_mass_kg,
            "mass_closure_error_kg": ca.allocator.closure_error_kg(),
            "air_exposed_area_m2": air_area,
            "surface_to_volume_ratio_m_inv": air_area / (np.count_nonzero(final_mask) * voxel_m**3),
            "runtime_cumulative_s": perf_counter() - started,
            **cache_check,
            **surface,
            **extents,
        }
        for family in FAMILY_NAMES:
            column = FAMILY_COLUMNS[family]
            count = family_counts[family]
            mean_area = 0.5 * (
                float(previous_surface[f"surface_area_{column}_m2"])
                + float(surface[f"surface_area_{column}_m2"])
            )
            row[f"new_events_{column}"] = count
            row[f"new_event_fraction_{column}"] = count / delta if delta else 0.0
            row[f"area_normalized_events_{column}_m_minus2"] = count / mean_area if mean_area > 0.0 else 0.0
        stage_rows.append(row)
        previous_surface = surface
        previous_target = int(target_events)
    ca.assert_invariants()
    final_mask = state.material == Material.SALT
    return {
        "case_id": case["case_id"],
        "stage_rows": stage_rows,
        "event_rows": event_rows,
        "final_mask": final_mask,
    }


def locked_growth_matrix(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["growth_campaign"]
    fitting_um = float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"])
    cases = []
    for rotation_index, angles in enumerate(section["rotations_xyz_deg"]):
        for policy in config["edge_policies"]:
            for alpha in config["alpha_values_with_zero"]:
                for seed in section["rng_seeds"]:
                    cases.append(
                        {
                            "case_id": f"R{rotation_index}_a{float(alpha):+.12g}_{policy}_seed{seed}",
                            "rotation_index": rotation_index,
                            "rotation_xyz_deg": angles,
                            "edge_policy": policy,
                            "alpha": alpha,
                            "rng_seed": seed,
                            "shape": section["shape"],
                            "stage_fractions": section["stage_fractions"],
                            "stage_events": section["stage_events"],
                            "voxel_size_um": config["voxel_size_um"],
                            "grain_size_um": config["grain_size_um"],
                            "fitting_radius_um": fitting_um,
                            "quality_minimum": config["quality_minimum"],
                            "near_family_tolerance_deg": config["near_family_tolerance_deg"],
                            "solid_density_kg_m3": section["solid_density_kg_m3"],
                        }
                    )
    with ProcessPoolExecutor(max_workers=min(int(section["maximum_workers"]), len(cases))) as executor:
        return list(executor.map(locked_growth_worker, cases))


def event_points(rows: Sequence[dict[str, object]], case_id: str, stage_fraction: float) -> np.ndarray:
    selected = [
        row
        for row in rows
        if row["case_id"] == case_id and float(row["stage_fraction"]) <= stage_fraction
    ]
    return np.asarray(
        [[row["point_crystal_x_m"], row["point_crystal_y_m"], row["point_crystal_z_m"]] for row in selected],
        dtype=float,
    )


def family_fraction_vector(rows: Sequence[dict[str, object]], case_id: str, stage_fraction: float) -> np.ndarray:
    selected = [
        row
        for row in rows
        if row["case_id"] == case_id and float(row["stage_fraction"]) <= stage_fraction
    ]
    counts = np.asarray(
        [sum(row["orientation_family"] == family for row in selected) for family in FAMILY_NAMES],
        dtype=float,
    )
    return counts / np.sum(counts)


def rotation_and_signal_analysis(
    config: dict[str, object],
    stage_rows: list[dict[str, object]],
    event_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    section = config["growth_campaign"]
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    stage_index = {(row["case_id"], float(row["stage_fraction"])): row for row in stage_rows}
    rotation_rows = []
    for policy in config["edge_policies"]:
        for alpha in config["alpha_values_with_zero"]:
            for seed in section["rng_seeds"]:
                for stage in section["stage_fractions"]:
                    ids = [f"R{rotation}_a{float(alpha):+.12g}_{policy}_seed{seed}" for rotation in (0, 1)]
                    points = [event_points(event_rows, case_id, float(stage)) for case_id in ids]
                    fractions = [family_fraction_vector(event_rows, case_id, float(stage)) for case_id in ids]
                    extents = []
                    for case_id in ids:
                        row = stage_index[(case_id, float(stage))]
                        extents.append(
                            np.asarray([row[f"growth_extent_{axis}_total_m"] for axis in ("x", "y", "z")], dtype=float)
                        )
                    rotation_rows.append(
                        {
                            "edge_policy": policy,
                            "alpha": alpha,
                            "rng_seed": seed,
                            "stage_fraction": stage,
                            "growth_point_chamfer_voxels": symmetric_chamfer_distance(points[0], points[1]) / voxel_m,
                            "event_family_fraction_l1": float(np.sum(np.abs(fractions[0] - fractions[1]))),
                            "directional_extent_relative_difference": float(
                                np.linalg.norm(extents[0] - extents[1])
                                / max(np.linalg.norm(extents[0]), voxel_m)
                            ),
                        }
                    )
    signal_rows = []
    alpha_magnitude = float(config["alpha_magnitude"])
    for rotation in (0, 1):
        for policy in config["edge_policies"]:
            for alpha in (-alpha_magnitude, alpha_magnitude):
                for seed in section["rng_seeds"]:
                    for stage in section["stage_fractions"]:
                        candidate_id = f"R{rotation}_a{alpha:+.12g}_{policy}_seed{seed}"
                        baseline_id = f"R{rotation}_a{0.0:+.12g}_{policy}_seed{seed}"
                        candidate_points = event_points(event_rows, candidate_id, float(stage))
                        baseline_points = event_points(event_rows, baseline_id, float(stage))
                        candidate_fraction = family_fraction_vector(event_rows, candidate_id, float(stage))
                        baseline_fraction = family_fraction_vector(event_rows, baseline_id, float(stage))
                        signal_rows.append(
                            {
                                "rotation_index": rotation,
                                "edge_policy": policy,
                                "alpha": alpha,
                                "rng_seed": seed,
                                "stage_fraction": stage,
                                "growth_point_chamfer_voxels": symmetric_chamfer_distance(candidate_points, baseline_points) / voxel_m,
                                "event_family_fraction_l1": float(np.sum(np.abs(candidate_fraction - baseline_fraction))),
                            }
                        )
    stochastic_rows = []
    seeds = list(section["rng_seeds"])
    for rotation in (0, 1):
        for policy in config["edge_policies"]:
            for stage in section["stage_fractions"]:
                case_ids = [f"R{rotation}_a{0.0:+.12g}_{policy}_seed{seed}" for seed in seeds]
                for first_index in range(len(case_ids)):
                    for second_index in range(first_index + 1, len(case_ids)):
                        first = event_points(event_rows, case_ids[first_index], float(stage))
                        second = event_points(event_rows, case_ids[second_index], float(stage))
                        first_fraction = family_fraction_vector(event_rows, case_ids[first_index], float(stage))
                        second_fraction = family_fraction_vector(event_rows, case_ids[second_index], float(stage))
                        stochastic_rows.append(
                            {
                                "rotation_index": rotation,
                                "edge_policy": policy,
                                "stage_fraction": stage,
                                "seed_pair": f"{seeds[first_index]}-{seeds[second_index]}",
                                "growth_point_chamfer_voxels": symmetric_chamfer_distance(first, second) / voxel_m,
                                "event_family_fraction_l1": float(np.sum(np.abs(first_fraction - second_fraction))),
                            }
                        )
    alpha_zero_lookup = {
        (row["edge_policy"], row["rng_seed"], float(row["stage_fraction"])): row
        for row in rotation_rows
        if np.isclose(float(row["alpha"]), 0.0)
    }
    for row in rotation_rows:
        baseline = alpha_zero_lookup[(row["edge_policy"], row["rng_seed"], float(row["stage_fraction"]))]
        row["chamfer_excess_over_alpha_zero"] = float(row["growth_point_chamfer_voxels"]) - float(baseline["growth_point_chamfer_voxels"])
        row["family_l1_excess_over_alpha_zero"] = float(row["event_family_fraction_l1"]) - float(baseline["event_family_fraction_l1"])
    return rotation_rows, signal_rows, stochastic_rows


def anchored_control(config: dict[str, object]) -> tuple[dict[str, object], np.ndarray]:
    section = config["anchored_control"]
    shape = tuple(int(value) for value in section["shape"])
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    grain_m = float(config["grain_size_um"]) * 1.0e-6
    rotation = rotation_from_angles(section["rotation_xyz_deg"])
    state = continuous_cube_state(
        shape=shape,
        substrate_layers=int(section["substrate_layers"]),
        grain_size_m=grain_m,
        voxel_size_m=voxel_m,
        rotation=rotation,
        anchored=True,
    )
    anchor_coordinates = np.argwhere(state.liquid_accessible)
    mass_parameters = EvaporationParameters(
        temperature_k=float(section["temperature_k"]),
        water_activity=float(section["water_activity"]),
        relative_humidity=float(section["relative_humidity"]),
        mass_transfer_coefficient_m_s=float(section["mass_transfer_coefficient_m_s"]),
        evaporation_area_m2=float(section["evaporation_area_m2"]),
        timestep_s=float(section["duration_s"]),
        salt_mass_fraction=float(section["salt_mass_fraction"]),
        solid_density_kg_m3=float(section["solid_density_kg_m3"]),
    )
    budget = compute_mass_budget(mass_parameters)
    model = OrientationCandidateWeightModel(
        voxel_size_m=voxel_m,
        fitting_radius_m=float(config["fitting_radius_over_grain_size"]) * grain_m,
        grain_to_lab=rotation,
        alpha=float(section["alpha"]),
        edge_policy=str(section["edge_policy"]),
        update_strategy="local",
        quality_minimum=float(config["quality_minimum"]),
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=float(section["solid_density_kg_m3"]),
        voxel_size_m=voxel_m,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=int(section["rng_seed"]),
        candidate_weight_model=model,
    )
    initial_mask = state.seed_mask.copy()
    initial_empty_mask = state.material == Material.EMPTY
    initial_air_area = air_exposed_face_area_m2(initial_mask, initial_empty_mask, voxel_m)
    started = perf_counter()
    result = ca.grow(budget.precipitated_salt_mass_kg)
    runtime = perf_counter() - started
    final_mask = state.material == Material.SALT
    salt_coordinates = np.argwhere(initial_mask)
    contact_count = int(np.count_nonzero(initial_mask[:, :, int(section["substrate_layers"])]))
    air_area = air_exposed_face_area_m2(final_mask, state.material == Material.EMPTY, voxel_m)
    return (
        {
            "scenario": "substrate_anchored_single_supply",
            "rotation_xyz_deg": json.dumps(section["rotation_xyz_deg"]),
            "substrate_layers": section["substrate_layers"],
            "substrate_contact_voxels": contact_count,
            "minimum_seed_z_index": int(np.min(salt_coordinates[:, 2])),
            "liquid_supply_anchor_count": int(anchor_coordinates.shape[0]),
            "liquid_supply_anchor_coordinate": json.dumps(tuple(int(value) for value in anchor_coordinates[0])),
            "air_exposed_area_initial_m2": initial_air_area,
            "air_exposed_area_final_m2": air_area,
            "temperature_k": section["temperature_k"],
            "water_activity": section["water_activity"],
            "relative_humidity": section["relative_humidity"],
            "duration_s": section["duration_s"],
            "evaporation_area_m2": section["evaporation_area_m2"],
            "evaporation_regime": budget.regime,
            "delivered_salt_mass_kg": budget.precipitated_salt_mass_kg,
            "requested_solid_voxels": result.transaction.requested_solid_voxels,
            "realized_solid_voxels": result.transaction.realized_solid_voxels,
            "realized_salt_mass_kg": result.transaction.committed_solid_mass_kg,
            "deferred_mass_kg": result.transaction.deferred_precipitation_mass_kg,
            "mass_closure_error_kg": ca.allocator.closure_error_kg(),
            "supply_connected_salt_voxels": result.morphology.supply_connected_salt_voxels,
            "total_salt_voxels": result.morphology.total_salt_voxels,
            "runtime_s": runtime,
            "film_hydraulics_solved": False,
        },
        final_mask,
    )


def locked_stage(config: dict[str, object], output: Path) -> dict[str, object]:
    candidate_rows, probability_rows = fixed_front_analysis(
        config, sample_counts=[int(config["fixed_front"]["sample_count"])]
    )
    write_csv(output / "fixed_front_candidates.csv", candidate_rows)
    write_csv(output / "fixed_front_probability_groups.csv", probability_rows)
    attribution_rows = cube_attribution(config)
    write_csv(output / "cube_site_attribution.csv", attribution_rows)
    cache_rows = cache_case_matrix(config, "cache_verification")
    write_csv(output / "cache_full_vs_local.csv", cache_rows)
    targeted_rows = targeted_cache_scenarios(config)
    write_csv(output / "cache_targeted_insertions.csv", targeted_rows)

    growth_results = locked_growth_matrix(config)
    stage_rows = [row for result in growth_results for row in result["stage_rows"]]
    event_rows = [row for result in growth_results for row in result["event_rows"]]
    write_csv(output / "growth_stages.csv", stage_rows)
    write_csv(output / "growth_events.csv", event_rows)
    rotation_rows, signal_rows, stochastic_rows = rotation_and_signal_analysis(
        config, stage_rows, event_rows
    )
    write_csv(output / "growth_rotation_covariance.csv", rotation_rows)
    write_csv(output / "anisotropy_signal.csv", signal_rows)
    write_csv(output / "alpha_zero_stochastic_variation.csv", stochastic_rows)

    anchored_row, anchored_mask = anchored_control(config)
    write_csv(output / "anchored_control.csv", [anchored_row])
    final_masks = {result["case_id"]: result["final_mask"] for result in growth_results}
    representative_ids = config["figures"]["representative_case_ids"]
    np.savez_compressed(
        output / "representative_final_masks.npz",
        **{case_id: final_masks[case_id] for case_id in representative_ids},
        anchored_control=anchored_mask,
    )

    probability_pass = all(
        float(row["absolute_error"])
        <= float(config["acceptance"]["fixed_front_max_absolute_probability_error"])
        and abs(float(row["standardized_error"]))
        <= float(config["acceptance"]["fixed_front_max_standardized_error"])
        for row in probability_rows
    )
    cache_pass = all(
        row["selected_sequence_exact"]
        and row["final_mask_exact"]
        and float(row["selected_weight_max_abs_error"])
        <= float(config["acceptance"]["cache_weight_tolerance"])
        for row in cache_rows
    ) and all(row["pass"] for row in targeted_rows)
    periodic_cache_pass = all(
        float(row["maximum_weight_error"]) <= float(config["acceptance"]["cache_weight_tolerance"])
        and int(row["site_class_mismatch_count"]) == 0
        and int(row["anisotropy_eligibility_mismatch_count"]) == 0
        for row in stage_rows
    )
    mass_pass = all(
        abs(float(row["mass_closure_error_kg"]))
        <= closure_tolerance_kg(
            float(row["cumulative_input_mass_kg"]),
            float(config["growth_campaign"]["solid_density_kg_m3"])
            * (float(config["voxel_size_um"]) * 1.0e-6) ** 3,
        )
        for row in stage_rows
    ) and float(anchored_row["mass_closure_error_kg"]) == 0.0

    final_stage = max(float(value) for value in config["growth_campaign"]["stage_fractions"])
    alpha_zero_rotation = [
        float(row["growth_point_chamfer_voxels"])
        for row in rotation_rows
        if np.isclose(float(row["alpha"]), 0.0) and np.isclose(float(row["stage_fraction"]), final_stage)
    ]
    stochastic = [
        float(row["growth_point_chamfer_voxels"])
        for row in stochastic_rows
        if np.isclose(float(row["stage_fraction"]), final_stage)
    ]
    baseline = max(float(np.median(alpha_zero_rotation)), float(np.median(stochastic)))
    signal_by_condition = {}
    signal_pass = True
    for policy in config["edge_policies"]:
        for alpha in (-float(config["alpha_magnitude"]), float(config["alpha_magnitude"])):
            values = [
                float(row["growth_point_chamfer_voxels"])
                for row in signal_rows
                if row["edge_policy"] == policy
                and np.isclose(float(row["alpha"]), alpha)
                and np.isclose(float(row["stage_fraction"]), final_stage)
            ]
            median_signal = float(np.median(values))
            ratio = median_signal / baseline if baseline > 0.0 else float("inf")
            signal_by_condition[f"{policy}|{alpha:+.12g}"] = {
                "median_chamfer_voxels": median_signal,
                "baseline_chamfer_voxels": baseline,
                "signal_to_baseline_ratio": ratio,
            }
            signal_pass &= ratio >= float(config["acceptance"]["morphology_signal_to_baseline_ratio_min"])
    rotation_pass = all(
        float(row["chamfer_excess_over_alpha_zero"])
        <= float(config["acceptance"]["rotation_chamfer_excess_max_voxels"])
        for row in rotation_rows
        if not np.isclose(float(row["alpha"]), 0.0)
    )
    anchored_pass = (
        int(anchored_row["substrate_contact_voxels"]) > 0
        and int(anchored_row["liquid_supply_anchor_count"]) == 1
        and int(anchored_row["realized_solid_voxels"]) > 0
    )
    gates = {
        "fixed_front_probability_agreement": probability_pass,
        "cache_long_sequence_exactness": cache_pass,
        "cache_periodic_full_checks": periodic_cache_pass,
        "mass_closure": mass_pass,
        "rotation_covariance_growth_increment": rotation_pass,
        "morphology_signal_beyond_alpha_zero_grid_and_stochastic_baseline": signal_pass,
        "anchored_supply_control": anchored_pass,
        "morphology_signal_summary": signal_by_condition,
    }
    gates["overall_pass"] = all(
        value for key, value in gates.items() if key != "morphology_signal_summary"
    )
    gates["decision"] = (
        "GO for complete single-grain numerical morphology model"
        if gates["overall_pass"]
        else "NO-GO pending failed Phase 2A.3b-1 gates"
    )
    write_json(output / "gate_summary.json", gates)
    plot_locked(config, output, probability_rows, attribution_rows, stage_rows, signal_rows, final_masks, anchored_mask)
    return gates


def plot_locked(
    config,
    output,
    probability_rows,
    attribution_rows,
    stage_rows,
    signal_rows,
    final_masks,
    anchored_mask,
):
    mpl.rcParams.update(
        {"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"}
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.25), layout="constrained")
    selected = [
        row for row in probability_rows
        if row["geometry"] == "sphere"
        and row["edge_policy"] == "smooth_average"
        and float(row["alpha"]) > 0.0
    ]
    x = np.arange(len(FAMILY_NAMES))
    expected = [next(float(row["expected_probability"]) for row in selected if row["orientation_family"] == family) for family in FAMILY_NAMES]
    observed = [next(float(row["empirical_probability"]) for row in selected if row["orientation_family"] == family) for family in FAMILY_NAMES]
    axes[0].bar(x - 0.18, expected, width=0.36, label="expected")
    axes[0].bar(x + 0.18, observed, width=0.36, label="sampled")
    axes[0].set(xticks=x, xticklabels=FAMILY_NAMES, ylabel="Fixed-front selection probability")
    axes[0].legend(frameon=False)
    selected_attr = [
        row for row in attribution_rows
        if int(row["cells_across_grain"]) == max(config["cube_attribution"]["cells_across"])
        and row["edge_policy"] == "smooth_average"
        and float(row["alpha"]) > 0.0
    ]
    axes[1].bar(
        [row["site_class"].replace("_like", "") for row in selected_attr],
        [float(row["fraction_of_total_anisotropy_shift"]) for row in selected_attr],
        color=["#0072B2", "#E69F00", "#D55E00"],
    )
    axes[1].set(ylabel="Fraction of probability shift", xlabel="Digital candidate class")
    final_stage = max(config["growth_campaign"]["stage_fractions"])
    for policy, marker in zip(config["edge_policies"], ("o", "s"), strict=True):
        values = []
        labels = []
        for alpha in (-config["alpha_magnitude"], config["alpha_magnitude"]):
            subset = [
                float(row["growth_point_chamfer_voxels"])
                for row in signal_rows
                if row["edge_policy"] == policy
                and np.isclose(float(row["alpha"]), float(alpha))
                and np.isclose(float(row["stage_fraction"]), float(final_stage))
            ]
            values.append(float(np.median(subset)))
            labels.append("negative" if alpha < 0 else "positive")
        axes[2].plot(labels, values, marker=marker, label=policy)
    axes[2].set(ylabel="Median growth-only Chamfer vs alpha zero (voxels)", xlabel="Anisotropy sign")
    axes[2].legend(frameon=False, fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        path = output / f"mechanism_summary.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white", metadata={"Date": None} if suffix == "svg" else None)
        if suffix == "svg":
            normalize_svg(path)
    plt.close(fig)

    ids = config["figures"]["representative_case_ids"]
    masks = [final_masks[case_id] for case_id in ids] + [anchored_mask]
    free_initial_masks = []
    for case_id in ids:
        rotation_index = int(case_id.split("_", maxsplit=1)[0][1:])
        free_initial_masks.append(
            continuous_cube_state(
                shape=tuple(config["growth_campaign"]["shape"]),
                substrate_layers=1,
                grain_size_m=float(config["grain_size_um"]) * 1.0e-6,
                voxel_size_m=float(config["voxel_size_um"]) * 1.0e-6,
                rotation=rotation_from_angles(
                    config["growth_campaign"]["rotations_xyz_deg"][rotation_index]
                ),
            ).seed_mask
        )
    anchored_initial = continuous_cube_state(
        shape=tuple(config["anchored_control"]["shape"]),
        substrate_layers=int(config["anchored_control"]["substrate_layers"]),
        grain_size_m=float(config["grain_size_um"]) * 1.0e-6,
        voxel_size_m=float(config["voxel_size_um"]) * 1.0e-6,
        rotation=rotation_from_angles(config["anchored_control"]["rotation_xyz_deg"]),
        anchored=True,
    ).seed_mask
    initial_masks = free_initial_masks + [anchored_initial]
    titles = [r"$\alpha<0$", r"$\alpha=0$", r"$\alpha>0$", "anchored control"]
    increment_cmap = mpl.colors.ListedColormap(["white", "#BDBDBD", "#D55E00"])
    increment_norm = mpl.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5], increment_cmap.N)
    fig, axes = plt.subplots(1, len(masks), figsize=(3.0 * len(masks), 3.0), layout="constrained")
    for axis, initial_mask, mask, title in zip(axes, initial_masks, masks, titles, strict=True):
        initial_projection = np.max(initial_mask, axis=1).T
        growth_projection = np.max(mask & ~initial_mask, axis=1).T
        display = initial_projection.astype(np.uint8)
        display[growth_projection] = 2
        axis.imshow(
            display,
            origin="lower",
            cmap=increment_cmap,
            norm=increment_norm,
            interpolation="nearest",
        )
        axis.set(title=title, xlabel="x voxel", ylabel="z voxel")
    axes[0].legend(
        handles=[
            Patch(facecolor="#BDBDBD", label="initial seed"),
            Patch(facecolor="#D55E00", label="new deposit"),
        ],
        frameon=False,
        fontsize=7,
        loc="upper left",
    )
    for suffix in ("png", "svg"):
        path = output / f"long_growth_morphologies.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white", metadata={"Date": None} if suffix == "svg" else None)
        if suffix == "svg":
            normalize_svg(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("development", "locked"), required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    expected_stage_labels = {
        "development": {"development", "method_development"},
        "locked": {"locked", "locked_verification"},
    }
    if config["stage"] not in expected_stage_labels[args.stage]:
        raise ValueError("config stage does not match requested stage")
    started = perf_counter()
    gates = (
        development_stage(config, args.output)
        if args.stage == "development"
        else locked_stage(config, args.output)
    )
    manifest = {
        "stage": args.stage,
        "base_commit": config["base_commit"],
        "config": args.config.resolve().relative_to(REPOSITORY_ROOT.resolve()).as_posix(),
        "config_sha256": sha256(args.config),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "runtime_s": perf_counter() - started,
        "maximum_workers": (
            config["cache_development"]["maximum_workers"]
            if args.stage == "development"
            else config["growth_campaign"]["maximum_workers"]
        ),
        "physical_validation_claimed": False,
    }
    write_json(args.output / "manifest.json", manifest)
    write_checksums(args.output)
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
