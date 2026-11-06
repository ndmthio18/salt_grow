"""Run the bounded Phase 2A.3b-0 orientation/CA interface validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from time import perf_counter

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase2 import (  # noqa: E402
    Material,
    OrientationCandidateWeightModel,
    PhysicalState,
    SaltGrowthCA,
    map_candidate_normals,
    proper_cubic_rotations,
    rotation_matrix_xyz_degrees,
)
from salt_phase2.orientation_verification import (  # noqa: E402
    ReferenceGeometry,
    evaluate_normal_case,
)
from salt_phase2.surface_orientation import cubic_attachment_weight, cubic_q4  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def continuous_cube_state(
    *,
    shape: tuple[int, int, int],
    substrate_layers: int,
    grain_size_m: float,
    voxel_size_m: float,
    rotation: np.ndarray,
) -> PhysicalState:
    """Voxelize one continuous rotated cube without changing its physical size."""

    coordinates = [
        (-0.5 * count + 0.5 + np.arange(count)) * voxel_size_m
        for count in shape
    ]
    x, y, z = np.meshgrid(*coordinates, indexing="ij", sparse=True)
    local = [
        sum(component * rotation[j, i] for j, component in enumerate((x, y, z)))
        for i in range(3)
    ]
    seed = np.asarray(
        np.maximum.reduce([np.abs(component) for component in local])
        <= 0.5 * grain_size_m,
        dtype=bool,
    )
    material = np.full(shape, Material.EMPTY, dtype=np.uint8)
    material[:, :, :substrate_layers] = Material.SUBSTRATE
    if np.any(seed & (material == Material.SUBSTRATE)):
        raise RuntimeError("continuous seed overlaps the prescribed substrate")
    material[seed] = Material.SALT
    state = PhysicalState(
        material=material,
        liquid_accessible=seed.copy(),
        seed_mask=seed.copy(),
        substrate_layers=substrate_layers,
    )
    state.validate()
    return state


def candidate_coordinates(state: PhysicalState) -> list[tuple[int, int, int]]:
    salt = state.material == Material.SALT
    kernel = np.zeros((3, 3, 3), dtype=np.uint8)
    kernel[0, 1, 1] = kernel[2, 1, 1] = 1
    kernel[1, 0, 1] = kernel[1, 2, 1] = 1
    kernel[1, 1, 0] = kernel[1, 1, 2] = 1
    adjacent = ndimage.convolve(salt.astype(np.uint8), kernel, mode="constant") > 0
    mask = adjacent & (state.material == Material.EMPTY)
    return [tuple(int(value) for value in row) for row in np.argwhere(mask)]


def scale_transfer_matrix(config: dict[str, object]) -> list[dict[str, object]]:
    settings = config["physical_scale_transfer"]
    rows: list[dict[str, object]] = []
    for grain_um in settings["grain_sizes_um"]:
        grain_m = float(grain_um) * 1.0e-6
        for cells_across in settings["cells_across"]:
            voxel_m = grain_m / int(cells_across)
            for lambda_value in settings["fitting_radius_over_grain_size"]:
                fitting_m = float(lambda_value) * grain_m
                for rotation_index, angles in enumerate(settings["rotations_xyz_deg"]):
                    geometry = ReferenceGeometry(
                        name=f"cube_L{int(grain_um)}_R{rotation_index}",
                        kind="cube",
                        domain_length_m=3.0 * grain_m,
                        half_width_m=0.5 * grain_m,
                        rotation_xyz_deg=tuple(float(value) for value in angles),
                    )
                    record, _, classes = evaluate_normal_case(
                        geometry,
                        voxel_m,
                        "weighted_plane_fit",
                        fitting_m,
                        max_points=int(settings["maximum_sample_points"]),
                    )
                    sample_count = int(record["sample_count"])
                    rows.append(
                        {
                            "grain_size_um": grain_um,
                            "cells_across_grain": cells_across,
                            "voxel_size_um": voxel_m * 1.0e6,
                            "lambda_rfit_over_L": lambda_value,
                            "fitting_radius_um": fitting_m * 1.0e6,
                            "fitting_radius_cells": fitting_m / voxel_m,
                            "rotation_index": rotation_index,
                            "rotation_xyz_deg": json.dumps(angles),
                            "face_interior_samples": int(record["unique_normal_count"]),
                            "ambiguous_edge_corner_fraction": float(
                                np.mean(np.isin(classes, ["edge", "corner"]))
                            ),
                            "median_face_error_deg": record["median_angular_error_deg"],
                            "p90_face_error_deg": record["p90_angular_error_deg"],
                            "maximum_face_error_deg": record["max_angular_error_deg"],
                            "surface_area_relative_error": record[
                                "surface_area_relative_error"
                            ],
                            "volume_relative_error": record["volume_relative_error"],
                            "neighbourhood_spans_entire_grain": bool(fitting_m >= grain_m),
                            "sample_count": sample_count,
                            "runtime_s": record["runtime_s"],
                        }
                    )
    return rows


def corrected_resolution_gate(output: Path) -> list[dict[str, object]]:
    source = ROOT / "outputs" / "phase2a3a_locked" / "resolution_gate.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        original = list(csv.DictReader(handle))
    rows: list[dict[str, object]] = []
    for row in original:
        p90 = row["p90_face_normal_error_deg"]
        corrected = bool(
            int(row["cells_across_seed"]) > 2
            and row["normal_estimator_status"] == "evaluated"
            and float(row["valid_normal_fraction"]) >= 0.90
            and int(row["face_interior_samples"]) >= 100
            and p90 not in {"", "None"}
            and float(p90) <= 10.0
        )
        rows.append(
            {
                **row,
                "original_eligible_flag_preserved": row[
                    "eligible_as_facet_convergence_evidence"
                ],
                "corrected_eligible_flag": corrected,
                "correction_basis": (
                    "requires evaluated normals, valid fraction >=0.90, >=100 face samples, "
                    "and p90 <=10 deg; original locked thresholds/data unchanged"
                ),
            }
        )
    write_csv(output / "corrected_resolution_gate.csv", rows)
    metadata = {
        "original_locked_file": source.relative_to(REPOSITORY_ROOT).as_posix(),
        "original_locked_sha256": sha256(source),
        "original_locked_file_modified": False,
        "one_alignment_volume_result": "preserved only as the original locked diagnostic",
        "robust_alignment_rotation_audit": "remains a separately labelled post-lock diagnostic",
        "four_voxel_original_flag": next(
            row["original_eligible_flag_preserved"]
            for row in rows
            if int(row["cells_across_seed"]) == 4
        ),
        "four_voxel_corrected_flag": next(
            row["corrected_eligible_flag"]
            for row in rows
            if int(row["cells_across_seed"]) == 4
        ),
    }
    (output / "corrected_resolution_gate_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    return rows


def run_growth(
    config: dict[str, object],
    *,
    rotation_index: int,
    alpha: float,
    edge_policy: str,
    update_strategy: str,
    rng_seed: int,
    with_orientation_model: bool = True,
) -> tuple[dict[str, object], np.ndarray, tuple[tuple[int, int, int], ...], object | None]:
    campaign = config["growth_campaign"]
    shape = tuple(int(value) for value in campaign["shape"])
    voxel_m = float(campaign["voxel_size_um"]) * 1.0e-6
    grain_m = float(campaign["grain_size_um"]) * 1.0e-6
    fitting_m = float(campaign["fitting_radius_over_grain_size"]) * grain_m
    rotation = rotation_matrix_xyz_degrees(
        *campaign["rotations_xyz_deg"][rotation_index]
    )
    state = continuous_cube_state(
        shape=shape,
        substrate_layers=int(campaign["substrate_layers"]),
        grain_size_m=grain_m,
        voxel_size_m=voxel_m,
        rotation=rotation,
    )
    model = (
        OrientationCandidateWeightModel(
            voxel_size_m=voxel_m,
            fitting_radius_m=fitting_m,
            grain_to_lab=rotation,
            alpha=alpha,
            edge_policy=edge_policy,
            update_strategy=update_strategy,
            quality_minimum=float(config["normal_quality_minimum"]),
        )
        if with_orientation_model
        else None
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=float(campaign["solid_density_kg_m3"]),
        voxel_size_m=voxel_m,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=rng_seed,
        candidate_weight_model=model,
    )
    events = int(campaign["new_voxel_events"])
    started = perf_counter()
    result = ca.grow((events + 0.25) * ca.voxel_mass_kg)
    runtime = perf_counter() - started
    if result.transaction.realized_solid_voxels != events:
        raise RuntimeError("growth case did not realize the predeclared event count")
    final_mask = state.material == Material.SALT
    points_lab = np.argwhere(final_mask).astype(float) * voxel_m
    centred = points_lab - np.mean(points_lab, axis=0)
    covariance = centred.T @ centred / points_lab.shape[0]
    eigenvalues = np.sort(np.linalg.eigvalsh(covariance))[::-1]
    points_crystal = centred @ rotation
    extents = np.ptp(points_crystal, axis=0)
    selected_history = model.selection_history if model is not None else []
    class_counts = {
        name: sum(record.site_class == name for record in selected_history)
        for name in ("face_like", "edge_like", "corner_like")
    }
    candidate_summary = model.performance_summary() if model is not None else {}
    candidate_evaluations = int(candidate_summary.get("candidate_evaluations", 0))
    row = {
        "case_id": (
            f"R{rotation_index}_a{alpha:+.12g}_{edge_policy}_seed{rng_seed}"
            if with_orientation_model
            else f"R{rotation_index}_isotropic_reference_seed{rng_seed}"
        ),
        "rotation_index": rotation_index,
        "rotation_xyz_deg": json.dumps(campaign["rotations_xyz_deg"][rotation_index]),
        "alpha": alpha if with_orientation_model else 0.0,
        "edge_policy": edge_policy if with_orientation_model else "no_orientation_model",
        "update_strategy": update_strategy if with_orientation_model else "none",
        "rng_seed": rng_seed,
        "seed_voxels": state.initial_seed_voxel_count,
        "new_voxels": result.transaction.realized_solid_voxels,
        "total_voxels": int(np.count_nonzero(final_mask)),
        "committed_mass_kg": result.transaction.committed_solid_mass_kg,
        "retained_mass_kg": result.transaction.unallocated_mass_kg,
        "mass_closure_error_kg": ca.allocator.closure_error_kg(),
        "gyration_eigenvalue_1_m2": eigenvalues[0],
        "gyration_eigenvalue_2_m2": eigenvalues[1],
        "gyration_eigenvalue_3_m2": eigenvalues[2],
        "crystal_extent_x_m": extents[0],
        "crystal_extent_y_m": extents[1],
        "crystal_extent_z_m": extents[2],
        "selected_face_like": class_counts["face_like"],
        "selected_edge_like": class_counts["edge_like"],
        "selected_corner_like": class_counts["corner_like"],
        "selected_anisotropy_eligible": sum(
            record.anisotropy_eligible for record in selected_history
        ),
        "candidate_evaluations": candidate_evaluations,
        "candidate_face_like_fraction": (
            float(candidate_summary.get("face_like_evaluations", 0))
            / candidate_evaluations
            if candidate_evaluations
            else 0.0
        ),
        "candidate_edge_like_fraction": (
            float(candidate_summary.get("edge_like_evaluations", 0))
            / candidate_evaluations
            if candidate_evaluations
            else 0.0
        ),
        "candidate_corner_like_fraction": (
            float(candidate_summary.get("corner_like_evaluations", 0))
            / candidate_evaluations
            if candidate_evaluations
            else 0.0
        ),
        "candidate_anisotropy_eligible_fraction": (
            float(candidate_summary.get("anisotropy_eligible_evaluations", 0))
            / candidate_evaluations
            if candidate_evaluations
            else 0.0
        ),
        "normal_recomputed_candidates": int(
            candidate_summary.get("recomputed_candidates", 0)
        ),
        "normal_runtime_s": float(candidate_summary.get("normal_runtime_s", 0.0)),
        "total_runtime_s": runtime,
    }
    return row, final_mask.copy(), result.selected_coordinates, model


def performance_benchmark(config: dict[str, object]) -> tuple[list[dict[str, object]], str]:
    performance = config["performance"]
    adjusted = json.loads(json.dumps(config))
    copy_config["growth_campaign"].update(
        {
            "grain_size_um": performance["grain_size_um"],
            "voxel_size_um": performance["voxel_size_um"],
            "fitting_radius_over_grain_size": performance[
                "fitting_radius_over_grain_size"
            ],
            "rotations_xyz_deg": [performance["rotation_xyz_deg"]],
            "new_voxel_events": performance["events"],
        }
    )
    runs = {}
    rows = []
    for strategy in performance["strategies"]:
        row, mask, sequence, model = run_growth(
            adjusted,
            rotation_index=0,
            alpha=float(performance["alpha"]),
            edge_policy=str(performance["edge_policy"]),
            update_strategy=strategy,
            rng_seed=int(performance["rng_seed"]),
        )
        runs[strategy] = (row, mask, sequence, model)
    reference = runs["full"]
    for strategy in performance["strategies"]:
        row, mask, sequence, model = runs[strategy]
        selected_weight_error = float("nan")
        if sequence == reference[2]:
            selected_weight_error = float(
                np.max(
                    np.abs(
                        np.asarray(
                            [record.anisotropy_weight for record in model.selection_history]
                        )
                        - np.asarray(
                            [
                                record.anisotropy_weight
                                for record in reference[3].selection_history
                            ]
                        )
                    )
                )
            )
        rows.append(
            {
                "strategy": strategy,
                "total_runtime_s": row["total_runtime_s"],
                "normal_runtime_s": row["normal_runtime_s"],
                "recomputed_candidates": row["normal_recomputed_candidates"],
                "selected_sequence_equals_full": sequence == reference[2],
                "final_mask_equals_full": bool(np.array_equal(mask, reference[1])),
                "selected_weight_max_abs_error_vs_full": selected_weight_error,
                "eligible_for_selection": strategy in {"full", "local"},
            }
        )
    tolerance = float(config["acceptance"]["local_strategy_weight_tolerance"])
    eligible = [
        row
        for row in rows
        if row["eligible_for_selection"]
        and row["selected_sequence_equals_full"]
        and row["final_mask_equals_full"]
        and float(row["selected_weight_max_abs_error_vs_full"]) <= tolerance
    ]
    selected = min(eligible, key=lambda row: float(row["total_runtime_s"]))["strategy"]
    return rows, str(selected)


def candidate_snapshot_rows(config: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object]]:
    campaign = config["growth_campaign"]
    voxel_m = float(campaign["voxel_size_um"]) * 1.0e-6
    grain_m = float(campaign["grain_size_um"]) * 1.0e-6
    fitting_m = float(campaign["fitting_radius_over_grain_size"]) * grain_m
    rows: list[dict[str, object]] = []
    symmetry_errors: list[float] = []
    symmetry = proper_cubic_rotations()[7]
    for rotation_index, angles in enumerate(campaign["rotations_xyz_deg"]):
        rotation = rotation_matrix_xyz_degrees(*angles)
        state = continuous_cube_state(
            shape=tuple(campaign["shape"]),
            substrate_layers=int(campaign["substrate_layers"]),
            grain_size_m=grain_m,
            voxel_size_m=voxel_m,
            rotation=rotation,
        )
        coordinates = candidate_coordinates(state)
        salt = state.material == Material.SALT
        for policy in config["edge_policies"]:
            for alpha in config["alpha_values"]:
                records = map_candidate_normals(
                    salt,
                    coordinates,
                    voxel_size_m=voxel_m,
                    fitting_radius_m=fitting_m,
                    grain_to_lab=rotation,
                    alpha=float(alpha),
                    edge_policy=policy,
                    quality_minimum=float(config["normal_quality_minimum"]),
                )
                equivalent = map_candidate_normals(
                    salt,
                    coordinates,
                    voxel_size_m=voxel_m,
                    fitting_radius_m=fitting_m,
                    grain_to_lab=rotation @ symmetry,
                    alpha=float(alpha),
                    edge_policy=policy,
                    quality_minimum=float(config["normal_quality_minimum"]),
                )
                symmetry_errors.append(
                    float(
                        np.max(
                            np.abs(
                                np.asarray([r.anisotropy_weight for r in records])
                                - np.asarray([r.anisotropy_weight for r in equivalent])
                            )
                        )
                    )
                )
                for record in records:
                    rows.append(
                        {
                            "rotation_index": rotation_index,
                            "alpha": alpha,
                            "edge_policy": policy,
                            "coordinate": json.dumps(record.coordinate),
                            "site_class": record.site_class,
                            "adjacent_salt_faces": record.adjacent_salt_faces,
                            "valid_face_normals": record.valid_face_normals,
                            "normal_set": json.dumps(record.normal_set),
                            "combined_normal": json.dumps(record.combined_normal),
                            "normal_quality": record.normal_quality,
                            "anisotropy_eligible": record.anisotropy_eligible,
                            "anisotropy_weight": record.anisotropy_weight,
                        }
                    )
    directions = np.asarray([[1, 0, 0], [1, 1, 0], [1, 1, 1]], dtype=float)
    alpha = float(config["alpha_magnitude"])
    invariants = {
        "q4_100_110_111": cubic_q4(directions).tolist(),
        "positive_alpha_weights_100_110_111": cubic_attachment_weight(
            directions, np.eye(3), alpha
        ).tolist(),
        "negative_alpha_weights_100_110_111": cubic_attachment_weight(
            directions, np.eye(3), -alpha
        ).tolist(),
        "maximum_candidate_weight_error_under_cubic_symmetry": max(symmetry_errors),
    }
    return rows, invariants


def comparison_tables(
    config: dict[str, object],
    rows: list[dict[str, object]],
    masks: dict[str, np.ndarray],
    sequences: dict[str, tuple[tuple[int, int, int], ...]],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    index = {row["case_id"]: row for row in rows}
    alpha_zero: list[dict[str, object]] = []
    campaign = config["growth_campaign"]
    for rotation_index in range(len(campaign["rotations_xyz_deg"])):
        for seed in campaign["rng_seeds"]:
            reference_id = f"R{rotation_index}_isotropic_reference_seed{seed}"
            for policy in config["edge_policies"]:
                candidate_id = f"R{rotation_index}_a{0.0:+.12g}_{policy}_seed{seed}"
                alpha_zero.append(
                    {
                        "rotation_index": rotation_index,
                        "rng_seed": seed,
                        "edge_policy": policy,
                        "selected_sequence_exact": sequences[candidate_id]
                        == sequences[reference_id],
                        "final_mask_exact": bool(
                            np.array_equal(masks[candidate_id], masks[reference_id])
                        ),
                        "mass_exact": index[candidate_id]["committed_mass_kg"]
                        == index[reference_id]["committed_mass_kg"],
                    }
                )
    covariance: list[dict[str, object]] = []
    for policy in config["edge_policies"]:
        for alpha in config["alpha_values"]:
            for seed in campaign["rng_seeds"]:
                ids = [
                    f"R{rotation}_a{float(alpha):+.12g}_{policy}_seed{seed}"
                    for rotation in range(2)
                ]
                vectors = []
                for case_id in ids:
                    vector = np.asarray(
                        [
                            index[case_id][f"gyration_eigenvalue_{axis}_m2"]
                            for axis in (1, 2, 3)
                        ],
                        dtype=float,
                    )
                    vectors.append(vector / np.sum(vector))
                difference = float(
                    np.linalg.norm(vectors[0] - vectors[1])
                    / max(np.linalg.norm(vectors[0]), 1.0e-30)
                )
                covariance.append(
                    {
                        "edge_policy": policy,
                        "alpha": alpha,
                        "rng_seed": seed,
                        "joint_rotation_normalized_gyration_relative_difference": difference,
                        "threshold": config["acceptance"][
                            "joint_rotation_gyration_relative_difference_max"
                        ],
                        "pass": difference
                        <= float(
                            config["acceptance"][
                                "joint_rotation_gyration_relative_difference_max"
                            ]
                        ),
                    }
                )
    def paired_metrics(first_id: str, second_id: str) -> dict[str, object]:
        first_mask = masks[first_id]
        second_mask = masks[second_id]
        intersection = int(np.count_nonzero(first_mask & second_mask))
        union = int(np.count_nonzero(first_mask | second_mask))
        first_sequence = sequences[first_id]
        second_sequence = sequences[second_id]
        first_gyration = np.asarray(
            [index[first_id][f"gyration_eigenvalue_{axis}_m2"] for axis in (1, 2, 3)],
            dtype=float,
        )
        second_gyration = np.asarray(
            [index[second_id][f"gyration_eigenvalue_{axis}_m2"] for axis in (1, 2, 3)],
            dtype=float,
        )
        return {
            "first_case_id": first_id,
            "second_case_id": second_id,
            "selected_sequence_exact": first_sequence == second_sequence,
            "selected_coordinate_set_overlap_fraction": len(
                set(first_sequence) & set(second_sequence)
            )
            / len(first_sequence),
            "final_mask_jaccard": intersection / union,
            "final_mask_different_voxels": int(np.count_nonzero(first_mask != second_mask)),
            "gyration_relative_difference": float(
                np.linalg.norm(first_gyration - second_gyration)
                / max(np.linalg.norm(first_gyration), 1.0e-30)
            ),
        }

    policy_sensitivity: list[dict[str, object]] = []
    for rotation_index in range(len(campaign["rotations_xyz_deg"])):
        for alpha in config["alpha_values"]:
            for seed in campaign["rng_seeds"]:
                first_id = (
                    f"R{rotation_index}_a{float(alpha):+.12g}_smooth_average_seed{seed}"
                )
                second_id = (
                    f"R{rotation_index}_a{float(alpha):+.12g}_isotropic_ambiguous_seed{seed}"
                )
                policy_sensitivity.append(
                    {
                        "rotation_index": rotation_index,
                        "alpha": alpha,
                        "rng_seed": seed,
                        **paired_metrics(first_id, second_id),
                    }
                )

    anisotropy_effects: list[dict[str, object]] = []
    alpha_magnitude = float(config["alpha_magnitude"])
    for rotation_index in range(len(campaign["rotations_xyz_deg"])):
        for policy in config["edge_policies"]:
            for seed in campaign["rng_seeds"]:
                ids = {
                    "negative": f"R{rotation_index}_a{-alpha_magnitude:+.12g}_{policy}_seed{seed}",
                    "zero": f"R{rotation_index}_a{0.0:+.12g}_{policy}_seed{seed}",
                    "positive": f"R{rotation_index}_a{alpha_magnitude:+.12g}_{policy}_seed{seed}",
                }
                for comparison, first_key, second_key in (
                    ("negative_vs_zero", "negative", "zero"),
                    ("positive_vs_zero", "positive", "zero"),
                    ("positive_vs_negative", "positive", "negative"),
                ):
                    anisotropy_effects.append(
                        {
                            "rotation_index": rotation_index,
                            "edge_policy": policy,
                            "rng_seed": seed,
                            "comparison": comparison,
                            **paired_metrics(ids[first_key], ids[second_key]),
                        }
                    )
    return alpha_zero, covariance, policy_sensitivity, anisotropy_effects


def plot_results(
    output: Path,
    scale_rows: list[dict[str, object]],
    performance_rows: list[dict[str, object]],
    growth_rows: list[dict[str, object]],
    representative_masks: dict[str, np.ndarray],
) -> None:
    mpl.rcParams.update(
        {
            "font.size": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.2), layout="constrained")
    for grain in sorted({float(row["grain_size_um"]) for row in scale_rows}):
        for lambda_value in sorted(
            {float(row["lambda_rfit_over_L"]) for row in scale_rows}
        ):
            subset = [
                row
                for row in scale_rows
                if float(row["grain_size_um"]) == grain
                and float(row["lambda_rfit_over_L"]) == lambda_value
            ]
            x = sorted({int(row["cells_across_grain"]) for row in subset})
            y = [
                np.mean(
                    [
                        float(row["p90_face_error_deg"])
                        for row in subset
                        if int(row["cells_across_grain"]) == value
                    ]
                )
                for value in x
            ]
            axes[0].plot(x, y, "o-", label=f"L={grain:.0f} µm, λ={lambda_value:g}")
    axes[0].axhline(10.0, color="black", linestyle=":", label="10° gate")
    axes[0].set(xlabel="Cells across grain", ylabel="Mean face-normal p90 error (°)")
    axes[0].legend(frameon=False, fontsize=6)
    names = [row["strategy"] for row in performance_rows]
    axes[1].bar(
        names,
        [float(row["total_runtime_s"]) for row in performance_rows],
        color=["#0072B2", "#009E73", "#D55E00"],
    )
    axes[1].set(ylabel="Runtime for 8 growth events (s)", xlabel="Normal update strategy")
    modeled = [row for row in growth_rows if row["edge_policy"] != "no_orientation_model"]
    policies = list(dict.fromkeys(str(row["edge_policy"]) for row in modeled))
    for policy, marker in zip(policies, ("o", "s"), strict=True):
        means = []
        for alpha in (-float(np.log(2) * 1.5), 0.0, float(np.log(2) * 1.5)):
            selected = [
                float(row["candidate_anisotropy_eligible_fraction"])
                for row in modeled
                if row["edge_policy"] == policy and np.isclose(float(row["alpha"]), alpha)
            ]
            means.append(np.mean(selected))
        axes[2].plot([-1, 0, 1], means, marker=marker, label=policy)
    axes[2].set(
        xlabel="Anisotropy sign",
        ylabel="Anisotropy-eligible candidate fraction",
        xticks=[-1, 0, 1],
        xticklabels=["negative", "zero", "positive"],
    )
    axes[2].legend(frameon=False, fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"phase2a3b0_diagnostics.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)

    keys = list(representative_masks)
    fig, axes = plt.subplots(1, len(keys), figsize=(3.0 * len(keys), 3.0), layout="constrained")
    if len(keys) == 1:
        axes = [axes]
    for axis, key in zip(axes, keys, strict=True):
        mask = representative_masks[key]
        axis.imshow(np.max(mask, axis=1).T, origin="lower", cmap="viridis", interpolation="nearest")
        axis.set(title=key, xlabel="x voxel", ylabel="z voxel")
    for suffix in ("png", "svg"):
        fig.savefig(output / f"representative_growth.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "configs" / "phase2a3b0_protocol.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs" / "phase2a3b0",
    )
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))

    scale_rows = scale_transfer_matrix(config)
    write_csv(args.output / "scale_transfer.csv", scale_rows)
    corrected_resolution_gate(args.output)
    performance_rows, selected_strategy = performance_benchmark(config)
    write_csv(args.output / "performance.csv", performance_rows)

    snapshot_rows, invariants = candidate_snapshot_rows(config)
    write_csv(args.output / "candidate_snapshot_diagnostics.csv", snapshot_rows)
    (args.output / "anisotropy_invariants.json").write_text(
        json.dumps(invariants, indent=2), encoding="utf-8"
    )

    campaign = config["growth_campaign"]
    growth_rows: list[dict[str, object]] = []
    masks: dict[str, np.ndarray] = {}
    sequences: dict[str, tuple[tuple[int, int, int], ...]] = {}
    for rotation_index in range(len(campaign["rotations_xyz_deg"])):
        for rng_seed in campaign["rng_seeds"]:
            row, mask, sequence, _ = run_growth(
                config,
                rotation_index=rotation_index,
                alpha=0.0,
                edge_policy="no_orientation_model",
                update_strategy="none",
                rng_seed=int(rng_seed),
                with_orientation_model=False,
            )
            growth_rows.append(row)
            masks[row["case_id"]] = mask
            sequences[row["case_id"]] = sequence
            for policy in config["edge_policies"]:
                for alpha in config["alpha_values"]:
                    row, mask, sequence, _ = run_growth(
                        config,
                        rotation_index=rotation_index,
                        alpha=float(alpha),
                        edge_policy=policy,
                        update_strategy=selected_strategy,
                        rng_seed=int(rng_seed),
                    )
                    growth_rows.append(row)
                    masks[row["case_id"]] = mask
                    sequences[row["case_id"]] = sequence
    write_csv(args.output / "growth_cases.csv", growth_rows)
    alpha_zero, covariance, policy_sensitivity, anisotropy_effects = comparison_tables(
        config, growth_rows, masks, sequences
    )
    write_csv(args.output / "alpha_zero_comparisons.csv", alpha_zero)
    write_csv(args.output / "orientation_covariance.csv", covariance)
    write_csv(args.output / "edge_policy_sensitivity.csv", policy_sensitivity)
    write_csv(args.output / "anisotropy_effects.csv", anisotropy_effects)

    reproducibility_rows = []
    for rotation_index in range(len(campaign["rotations_xyz_deg"])):
        for policy in config["edge_policies"]:
            row, mask, sequence, _ = run_growth(
                config,
                rotation_index=rotation_index,
                alpha=float(config["alpha_magnitude"]),
                edge_policy=policy,
                update_strategy=selected_strategy,
                rng_seed=int(campaign["rng_seeds"][0]),
            )
            case_id = row["case_id"]
            reproducibility_rows.append(
                {
                    "case_id": case_id,
                    "selected_sequence_exact": sequence == sequences[case_id],
                    "final_mask_exact": bool(np.array_equal(mask, masks[case_id])),
                }
            )
    write_csv(args.output / "pinned_reproducibility.csv", reproducibility_rows)

    alpha = float(config["alpha_magnitude"])
    representative_ids = [
        f"R0_a{value:+.12g}_smooth_average_seed{campaign['rng_seeds'][0]}"
        for value in (-alpha, 0.0, alpha)
    ]
    representative_masks = {
        label: masks[case_id]
        for label, case_id in zip(("negative alpha", "alpha zero", "positive alpha"), representative_ids, strict=True)
    }
    np.savez_compressed(args.output / "representative_masks.npz", **representative_masks)
    plot_results(args.output, scale_rows, performance_rows, growth_rows, representative_masks)

    mass_ok = all(
        abs(float(row["mass_closure_error_kg"]))
        <= 32.0
        * math.ulp(
            max(abs(float(row["committed_mass_kg"])), np.finfo(float).tiny)
        )
        for row in growth_rows
    )
    scale_selected = [
        row
        for row in scale_rows
        if int(row["cells_across_grain"]) >= 16
        and np.isclose(
            float(row["lambda_rfit_over_L"]),
            float(campaign["fitting_radius_over_grain_size"]),
        )
    ]
    gates = {
        "mass_closure": mass_ok,
        "alpha_zero_exact": all(
            row["selected_sequence_exact"] and row["final_mask_exact"] and row["mass_exact"]
            for row in alpha_zero
        ),
        "pinned_reproducibility": all(
            row["selected_sequence_exact"] and row["final_mask_exact"]
            for row in reproducibility_rows
        ),
        "cubic_symmetry": invariants[
            "maximum_candidate_weight_error_under_cubic_symmetry"
        ]
        <= float(config["acceptance"]["cubic_symmetry_weight_tolerance"]),
        "scale_transfer": all(
            float(row["p90_face_error_deg"])
            <= float(config["acceptance"]["scale_transfer_face_p90_deg_max"])
            for row in scale_selected
        ),
        "joint_rotation_covariance": all(row["pass"] for row in covariance),
        "performance_strategy": selected_strategy,
        "four_voxel_corrected_ineligible": not next(
            row["corrected_eligible_flag"]
            for row in corrected_resolution_gate(args.output)
            if int(row["cells_across_seed"]) == 4
        ),
    }
    gates["overall_pass"] = all(
        value for key, value in gates.items() if key != "performance_strategy"
    )
    gates["decision"] = (
        "GO for full Phase 2A.3b single-grain numerical pilot"
        if gates["overall_pass"]
        else "NO-GO pending failed Phase 2A.3b-0 gates"
    )
    (args.output / "gate_summary.json").write_text(
        json.dumps(gates, indent=2), encoding="utf-8"
    )
    manifest = {
        "stage": config["stage"],
        "base_commit": config["base_commit"],
        "config": args.config.resolve().relative_to(REPOSITORY_ROOT.resolve()).as_posix(),
        "config_sha256": sha256(args.config),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "growth_case_count": len(growth_rows),
        "candidate_snapshot_row_count": len(snapshot_rows),
        "physical_validation_claimed": False,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    checksum_lines = []
    for path in sorted(args.output.iterdir()):
        if path.name != "SHA256SUMS.txt":
            checksum_lines.append(f"{sha256(path)}  {path.name}")
    (args.output / "SHA256SUMS.txt").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
