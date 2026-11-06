"""Run the pre-registered Phase 2A.3a development or locked verification set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT.parent / "phase1_mass_budget" / "src"))

from salt_phase2.orientation_verification import (  # noqa: E402
    evaluate_normal_case,
    geometry_catalogue,
)
from salt_phase2.surface_orientation import (  # noqa: E402
    angular_errors_deg,
    cubic_attachment_weight,
    cubic_q4,
    estimate_plane_fit_normals,
    estimate_sdf_normals,
    proper_cubic_rotations,
    rotation_matrix_xyz_degrees,
)


COLORS = {
    "signed_distance_gradient": "#0072B2",
    "weighted_plane_fit": "#D55E00",
}
MARKERS = {"signed_distance_gradient": "o", "weighted_plane_fit": "s"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot_development(rows: list[dict[str, object]], output: Path) -> None:
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), layout="constrained")
    for method in COLORS:
        subset = [row for row in rows if row["method"] == method]
        lengths = sorted({float(row["physical_length_um"]) for row in subset})
        medians = []
        p90s = []
        for length in lengths:
            selected = [row for row in subset if float(row["physical_length_um"]) == length]
            medians.append(np.mean([float(row["median_angular_error_deg"]) for row in selected]))
            p90s.append(np.mean([float(row["p90_angular_error_deg"]) for row in selected]))
        label = method.replace("_", " ")
        axes[0].plot(lengths, medians, marker=MARKERS[method], color=COLORS[method], label=label)
        axes[1].plot(lengths, p90s, marker=MARKERS[method], color=COLORS[method], label=label)
    axes[0].set(xlabel="Physical smoothing/fitting length (µm)", ylabel="Mean median angular error (°)")
    axes[1].set(xlabel="Physical smoothing/fitting length (µm)", ylabel="Mean 90th-percentile error (°)")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"method_development.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def _plot_locked(rows: list[dict[str, object]], output: Path) -> None:
    cases = list(dict.fromkeys(str(row["case"]) for row in rows))
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.35), layout="constrained")
    palette = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
    markers = ["o", "s", "^", "D"]
    for case, color, marker in zip(cases, palette, markers, strict=True):
        subset = sorted(
            [row for row in rows if row["case"] == case],
            key=lambda row: float(row["voxel_size_um"]),
            reverse=True,
        )
        x = [float(row["voxel_size_um"]) for row in subset]
        axes[0].plot(x, [float(row["median_angular_error_deg"]) for row in subset], marker=marker, color=color, label=case)
        axes[0].plot(x, [float(row["p90_angular_error_deg"]) for row in subset], marker=marker, color=color, linestyle="--")
        axes[1].plot(x, [100 * float(row["surface_area_relative_error"]) for row in subset], marker=marker, color=color)
        axes[2].plot(x, [100 * float(row["volume_relative_error"]) for row in subset], marker=marker, color=color)
    axes[0].set(xlabel="Voxel size (µm; finer to the right)", ylabel="Angular error (°)")
    axes[1].set(xlabel="Voxel size (µm; finer to the right)", ylabel="Surface-area error (%)")
    axes[2].set(xlabel="Voxel size (µm; finer to the right)", ylabel="Volume error (%)")
    for axis in axes:
        axis.invert_xaxis()
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=7)
    axes[0].text(0.02, 0.98, "solid: median\ndashed: p90", transform=axes[0].transAxes, va="top", fontsize=7)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"locked_refinement.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)


def _summarize_development(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, float], list[dict[str, object]]] = {}
    for row in rows:
        key = (str(row["method"]), float(row["physical_length_um"]))
        groups.setdefault(key, []).append(row)
    summary = []
    for (method, length), values in groups.items():
        summary.append(
            {
                "method": method,
                "physical_length_um": length,
                "mean_median_error_deg": float(np.mean([float(v["median_angular_error_deg"]) for v in values])),
                "mean_p90_error_deg": float(np.mean([float(v["p90_angular_error_deg"]) for v in values])),
                "mean_area_relative_error": float(np.mean([float(v["surface_area_relative_error"]) for v in values])),
                "total_runtime_s": float(np.sum([float(v["runtime_s"]) for v in values])),
            }
        )
    return sorted(summary, key=lambda row: (float(row["mean_p90_error_deg"]), float(row["mean_median_error_deg"]), float(row["total_runtime_s"])))


def _estimate(geometry, voxel_m: float, selection: dict[str, object], max_points: int):
    mask = geometry.voxelize(voxel_m)
    length_m = float(selection["physical_length_um"]) * 1.0e-6
    if selection["method"] == "weighted_plane_fit":
        estimate = estimate_plane_fit_normals(
            mask, voxel_m, length_m, max_points=max_points
        )
    else:
        estimate = estimate_sdf_normals(
            mask, voxel_m, length_m, max_points=max_points
        )
    reference, classes = geometry.reference_normals(
        estimate.points_m, max(2.0 * voxel_m, length_m)
    )
    return mask, estimate, reference, classes


def _anisotropy_invariants(config: dict[str, object]) -> dict[str, object]:
    directions = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]], dtype=float)
    q4 = cubic_q4(directions)
    alpha_values = [float(value) for value in config["alpha_values"]]
    weights = {
        f"{alpha:.16g}": cubic_attachment_weight(directions, np.eye(3), alpha).tolist()
        for alpha in alpha_values
    }
    probe = np.array([[0.2, -0.4, 0.7], [1.0, 2.0, -3.0]])
    alpha = max(alpha_values)
    baseline = cubic_attachment_weight(probe, np.eye(3), alpha)
    symmetry_error = max(
        float(np.max(np.abs(cubic_attachment_weight(probe, symmetry, alpha) - baseline)))
        for symmetry in proper_cubic_rotations()
    )
    return {
        "q4_100": float(q4[0]),
        "q4_110": float(q4[1]),
        "q4_111": float(q4[2]),
        "weights_by_alpha_100_110_111": weights,
        "proper_cubic_rotation_count": len(proper_cubic_rotations()),
        "maximum_cubic_symmetry_weight_error": symmetry_error,
        "interpretation": "Demonstration attachment weights only; not facet-growth velocities or calibrated NaCl kinetics.",
    }


def _rotation_covariance(
    config: dict[str, object], catalogue: dict[str, object]
) -> list[dict[str, object]]:
    selection = config["selected_estimator"]
    maximum = int(config["maximum_sample_points_per_case"])
    reference_geometry = catalogue[config["rotation_covariance"]["reference_geometry"]]
    rotated_geometry = catalogue[config["rotation_covariance"]["rotated_geometry"]]
    joint = rotated_geometry.rotation
    grain = rotation_matrix_xyz_degrees(*config["rotation_covariance"]["grain_rotation_xyz_deg"])
    alpha = max(float(value) for value in config["alpha_values"])
    rows = []
    for voxel_um in config["voxel_sizes_um"]:
        voxel_m = float(voxel_um) * 1.0e-6
        _, ref_estimate, ref_truth, ref_classes = _estimate(
            reference_geometry, voxel_m, selection, maximum
        )
        _, rot_estimate, _, rot_classes = _estimate(
            rotated_geometry, voxel_m, selection, maximum
        )
        ref_valid = ref_estimate.valid & (ref_classes == "face")
        rot_valid = rot_estimate.valid & (rot_classes == "face")
        tree = cKDTree(ref_estimate.points_m[ref_valid])
        ref_indices = np.flatnonzero(ref_valid)
        rotated_back_points = rot_estimate.points_m[rot_valid] @ joint
        distances, local_indices = tree.query(rotated_back_points, k=1)
        keep = distances <= float(config["rotation_covariance"]["matching_radius_voxels"]) * voxel_m
        matched_ref = ref_indices[local_indices[keep]]
        rotated_indices = np.flatnonzero(rot_valid)[keep]
        rotated_back_normals = rot_estimate.normals[rotated_indices] @ joint
        normal_errors = angular_errors_deg(
            rotated_back_normals, ref_estimate.normals[matched_ref]
        )
        isotropic_ref = cubic_attachment_weight(
            ref_estimate.normals[matched_ref], grain, 0.0
        )
        isotropic_rot = cubic_attachment_weight(
            rot_estimate.normals[rotated_indices], joint @ grain, 0.0
        )
        weight_ref = cubic_attachment_weight(
            ref_estimate.normals[matched_ref], grain, alpha
        )
        weight_fixed_rotated_crystal = cubic_attachment_weight(
            ref_estimate.normals[matched_ref], joint @ grain, alpha
        )
        weight_joint = cubic_attachment_weight(
            rot_estimate.normals[rotated_indices], joint @ grain, alpha
        )
        weight_difference = np.abs(weight_joint - weight_ref)
        exact_ref_weight = cubic_attachment_weight(ref_truth[matched_ref], grain, alpha)
        rows.append(
            {
                "voxel_size_um": voxel_um,
                "reference_face_samples": int(np.count_nonzero(ref_valid)),
                "rotated_face_samples": int(np.count_nonzero(rot_valid)),
                "matched_samples": int(np.count_nonzero(keep)),
                "matched_fraction_of_rotated_faces": float(np.mean(keep)),
                "matching_distance_p90_um": float(np.percentile(distances[keep], 90) * 1.0e6),
                "A_normal_median_error_deg": float(np.nanmedian(normal_errors)),
                "A_normal_p90_error_deg": float(np.nanpercentile(normal_errors, 90)),
                "A_isotropic_max_weight_error": float(np.max(np.abs(isotropic_rot - isotropic_ref))),
                "B_fixed_geometry_crystal_rotation_rms_change": float(
                    np.sqrt(np.mean((weight_fixed_rotated_crystal - weight_ref) ** 2))
                ),
                "B_reference_normal_weight_rmse": float(
                    np.sqrt(np.mean((weight_ref - exact_ref_weight) ** 2))
                ),
                "C_joint_rotation_weight_rmse": float(np.sqrt(np.mean(weight_difference**2))),
                "C_joint_rotation_weight_p90_abs_error": float(np.percentile(weight_difference, 90)),
            }
        )
    return rows


def _resolution_gate(config: dict[str, object], catalogue: dict[str, object]):
    selection = config["selected_estimator"]
    geometry = catalogue[config["resolution_gate"]["geometry"]]
    seed_width_m = 2.0 * geometry.half_width_m
    geometry = replace(geometry, domain_length_m=3.0 * seed_width_m)
    rows = []
    for cells_across in config["resolution_gate"]["cells_across_seed"]:
        voxel_m = seed_width_m / int(cells_across)
        domain_cells = int(round(geometry.domain_length_m / voxel_m))
        if not np.isclose(domain_cells * voxel_m, geometry.domain_length_m):
            raise RuntimeError("resolution-gate domain is not integral")
        mask = geometry.voxelize(voxel_m)
        volume = np.count_nonzero(mask) * voxel_m**3
        row = {
            "cells_across_seed": cells_across,
            "voxel_size_um": voxel_m * 1.0e6,
            "volume_relative_error": abs(volume / geometry.exact_volume_m3 - 1.0),
            "normal_estimator_status": "not_run",
            "valid_normal_fraction": 0.0,
            "face_interior_samples": 0,
            "p90_face_normal_error_deg": None,
            "eligible_as_facet_convergence_evidence": bool(cells_across > 2),
        }
        try:
            _, estimate, truth, classes = _estimate(
                geometry,
                voxel_m,
                selection,
                int(config["resolution_gate"]["maximum_sample_points"]),
            )
            errors = angular_errors_deg(estimate.normals, truth)
            face = estimate.valid & (classes == "face") & np.isfinite(errors)
            row.update(
                {
                    "normal_estimator_status": "evaluated",
                    "valid_normal_fraction": float(np.mean(estimate.valid)),
                    "face_interior_samples": int(np.count_nonzero(face)),
                    "p90_face_normal_error_deg": (
                        float(np.percentile(errors[face], 90)) if np.any(face) else None
                    ),
                }
            )
        except ValueError as error:
            row["normal_estimator_status"] = f"unresolved: {error}"
        rows.append(row)
    criteria = config["resolution_gate"]["criteria"]
    def minimum(predicate):
        matches = [int(row["cells_across_seed"]) for row in rows if predicate(row)]
        return min(matches) if matches else None
    requirements = {
        "minimum_cells_for_volume": minimum(
            lambda row: float(row["volume_relative_error"]) <= float(criteria["volume_relative_error_max"])
        ),
        "minimum_cells_for_surface_normals": minimum(
            lambda row: row["p90_face_normal_error_deg"] is not None
            and float(row["valid_normal_fraction"]) >= float(criteria["valid_normal_fraction_min"])
            and float(row["p90_face_normal_error_deg"]) <= float(criteria["surface_normal_p90_deg_max"])
        ),
        "minimum_cells_for_facet_orientation": minimum(
            lambda row: bool(row["eligible_as_facet_convergence_evidence"])
            and int(row["face_interior_samples"]) >= int(criteria["facet_face_samples_min"])
            and row["p90_face_normal_error_deg"] is not None
            and float(row["p90_face_normal_error_deg"]) <= float(criteria["facet_p90_deg_max"])
        ),
        "minimum_cells_across_narrow_gap": int(criteria["gap_cells_min"]),
        "gap_basis": "A one-cell boundary-location uncertainty gives <=20% relative width uncertainty only at five or more cells across; this is a geometry gate, not a flow-resolution claim.",
        "one_and_two_voxel_policy": "Reported as unresolved diagnostics only and excluded from facet-convergence evidence.",
    }
    return rows, requirements


def _locked_gates(
    rows: list[dict[str, object]],
    covariance: list[dict[str, object]],
    invariants: dict[str, object],
    config: dict[str, object],
) -> dict[str, object]:
    thresholds = config["tolerances"]
    finest = min(float(value) for value in config["voxel_sizes_um"])
    coarse = max(float(value) for value in config["voxel_sizes_um"])
    checks = []
    for case in config["locked_cases"]:
        fine = next(row for row in rows if row["case"] == case and float(row["voxel_size_um"]) == finest)
        coarse_row = next(row for row in rows if row["case"] == case and float(row["voxel_size_um"]) == coarse)
        case_checks = {
            "median": float(fine["median_angular_error_deg"]) <= float(thresholds["finest_median_angular_error_deg_max"]),
            "p90": float(fine["p90_angular_error_deg"]) <= float(thresholds["finest_p90_angular_error_deg_max"]),
            "maximum": float(fine["max_angular_error_deg"]) <= float(thresholds["finest_max_angular_error_deg_max"]),
            "area": float(fine["surface_area_relative_error"]) <= float(thresholds["finest_surface_area_relative_error_max"]),
            "volume": float(fine["volume_relative_error"]) <= float(thresholds["finest_volume_relative_error_max"]),
            "p90_refinement": float(fine["p90_angular_error_deg"]) <= float(coarse_row["p90_angular_error_deg"]) + float(thresholds["refinement_slack_deg"]),
            "area_refinement": float(fine["surface_area_relative_error"]) <= float(coarse_row["surface_area_relative_error"]),
        }
        checks.append({"case": case, "checks": case_checks, "pass": all(case_checks.values())})
    covariance_fine = next(row for row in covariance if float(row["voxel_size_um"]) == finest)
    covariance_checks = {
        "matched_fraction": float(covariance_fine["matched_fraction_of_rotated_faces"]) >= float(thresholds["rotation_matched_fraction_min"]),
        "normal_p90": float(covariance_fine["A_normal_p90_error_deg"]) <= float(thresholds["rotation_normal_p90_deg_max"]),
        "isotropic_exact": float(covariance_fine["A_isotropic_max_weight_error"]) <= float(thresholds["exact_invariant_tolerance"]),
        "joint_weight_rmse": float(covariance_fine["C_joint_rotation_weight_rmse"]) <= float(thresholds["joint_weight_rmse_max"]),
        "joint_weight_p90": float(covariance_fine["C_joint_rotation_weight_p90_abs_error"]) <= float(thresholds["joint_weight_p90_abs_error_max"]),
    }
    invariant_checks = {
        "q4_100": abs(float(invariants["q4_100"]) - 1.0) <= float(thresholds["exact_invariant_tolerance"]),
        "q4_110": abs(float(invariants["q4_110"]) - 0.5) <= float(thresholds["exact_invariant_tolerance"]),
        "q4_111": abs(float(invariants["q4_111"]) - 1.0 / 3.0) <= float(thresholds["exact_invariant_tolerance"]),
        "cubic_symmetry": float(invariants["maximum_cubic_symmetry_weight_error"]) <= float(thresholds["exact_invariant_tolerance"]),
    }
    overall = all(item["pass"] for item in checks) and all(covariance_checks.values()) and all(invariant_checks.values())
    return {
        "case_checks": checks,
        "rotation_covariance_checks": covariance_checks,
        "analytical_invariant_checks": invariant_checks,
        "overall_pass": overall,
        "decision": "GO for a bounded single-grain CA orientation pilot" if overall else "NO-GO for single-grain CA orientation coupling",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("development", "locked"), required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    expected = "method_development" if args.stage == "development" else "locked_verification"
    if config.get("stage") != expected:
        raise ValueError(f"config stage must be {expected}")
    args.output.mkdir(parents=True, exist_ok=False)
    catalogue = geometry_catalogue()
    rows: list[dict[str, object]] = []
    if args.stage == "development":
        for candidate in config["candidates"]:
            for length_um in candidate["physical_lengths_um"]:
                for case in config["cases"]:
                    for voxel_um in config["voxel_sizes_um"]:
                        record, _, _ = evaluate_normal_case(
                            catalogue[case],
                            voxel_um * 1.0e-6,
                            candidate["method"],
                            length_um * 1.0e-6,
                            max_points=int(config["maximum_sample_points_per_case"]),
                        )
                        rows.append(record)
        _write_csv(args.output / "development_results.csv", rows)
        ranking = _summarize_development(rows)
        _write_csv(args.output / "method_ranking.csv", ranking)
        _plot_development(rows, args.output)
    else:
        selection = config["selected_estimator"]
        for case in config["locked_cases"]:
            for voxel_um in config["voxel_sizes_um"]:
                record, _, _ = evaluate_normal_case(
                    catalogue[case],
                    voxel_um * 1.0e-6,
                    selection["method"],
                    selection["physical_length_um"] * 1.0e-6,
                    max_points=int(config["maximum_sample_points_per_case"]),
                )
                rows.append(record)
        _write_csv(args.output / "locked_results.csv", rows)
        _plot_locked(rows, args.output)
        invariants = _anisotropy_invariants(config)
        (args.output / "anisotropy_invariants.json").write_text(
            json.dumps(invariants, indent=2), encoding="utf-8"
        )
        covariance = _rotation_covariance(config, catalogue)
        _write_csv(args.output / "rotation_covariance.csv", covariance)
        resolution_rows, requirements = _resolution_gate(config, catalogue)
        _write_csv(args.output / "resolution_gate.csv", resolution_rows)
        (args.output / "resolution_requirements.json").write_text(
            json.dumps(requirements, indent=2), encoding="utf-8"
        )
        gates = _locked_gates(rows, covariance, invariants, config)
        (args.output / "gate_summary.json").write_text(
            json.dumps(gates, indent=2), encoding="utf-8"
        )
    manifest = {
        "stage": args.stage,
        "config": str(args.config),
        "config_sha256": _sha256(args.config),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": mpl.__version__,
        "record_count": len(rows),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    checksums = []
    for path in sorted(args.output.iterdir()):
        if path.name != "SHA256SUMS.txt":
            checksums.append(f"{_sha256(path)}  {path.name}")
    (args.output / "SHA256SUMS.txt").write_text("\n".join(checksums) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
