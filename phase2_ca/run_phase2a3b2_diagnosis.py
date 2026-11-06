"""Run Phase 2A.3b-2 exploratory or controlled anisotropy diagnosis."""

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
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import run_phase2a3b1_validation as phase2a3b1  # noqa: E402
from salt_phase2 import (  # noqa: E402
    Material,
    closure_tolerance_kg,
    crystallographic_candidate_diagnostics,
    equivalent_displacement_m,
    group_allocation,
    ideal_cubic_direction_families,
    support_advances_m,
    symmetric_chamfer_distance,
    uniform_cube_facet_displacement_m,
    volume_fraction_for_facet_displacement,
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def normalize_svg(path: Path) -> None:
    path.write_text(
        "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n",
        encoding="utf-8",
    )


def write_checksums(output: Path) -> None:
    lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (output / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def case_id(alpha: float, policy: str, seed: int, rotation_index: int = 0) -> str:
    return f"R{rotation_index}_a{float(alpha):+.12g}_{policy}_seed{seed}"


def analytical_resolution(config: dict[str, object]) -> list[dict[str, object]]:
    grain_m = float(config["grain_size_um"]) * 1.0e-6
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    rows = []
    for fraction in config["analytical_growth_fractions"]:
        displacement = uniform_cube_facet_displacement_m(grain_m, float(fraction))
        rows.append(
            {
                "added_volume_fraction": fraction,
                "facet_displacement_um": displacement * 1.0e6,
                "facet_displacement_voxels": displacement / voxel_m,
                "is_at_least_one_voxel_layer": displacement >= voxel_m * (1.0 - 1.0e-12),
                "interpretation": "uniform cubic expansion, not nonuniform CA deposition",
            }
        )
    exact_fraction = volume_fraction_for_facet_displacement(grain_m, voxel_m)
    if not any(np.isclose(float(row["added_volume_fraction"]), exact_fraction) for row in rows):
        raise RuntimeError("analytical protocol must include the exact one-layer fraction")
    return rows


def central_plane_mask(shape: tuple[int, int, int], normal: np.ndarray) -> np.ndarray:
    axes = [np.arange(count, dtype=float) - 0.5 * (count - 1) for count in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij", sparse=True)
    field = normal[0] * x + normal[1] * y + normal[2] * z
    return np.asarray(field <= 0.0, dtype=bool)


def exact_central_section_area_m2(family: str, side_m: float) -> float:
    factors = {"<100>": 1.0, "<110>": math.sqrt(2.0), "<111>": 3.0 * math.sqrt(3.0) / 4.0}
    return factors[family] * side_m**2


def planar_probability_velocity_audit(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["planar_audit"]
    side_m = float(section["physical_patch_side_um"]) * 1.0e-6
    fitting_m = float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"]) * 1.0e-6
    total_volume = float(section["fixed_total_deposited_volume_m3"])
    normals = {
        "<100>": np.asarray([1.0, 0.0, 0.0]),
        "<110>": np.asarray([1.0, 1.0, 0.0]) / math.sqrt(2.0),
        "<111>": np.asarray([1.0, 1.0, 1.0]) / math.sqrt(3.0),
    }
    rows: list[dict[str, object]] = []
    for voxel_um in section["voxel_sizes_um"]:
        voxel_m = float(voxel_um) * 1.0e-6
        cells = int(round(side_m / voxel_m))
        if not math.isclose(cells * voxel_m, side_m, rel_tol=0.0, abs_tol=1.0e-15):
            raise ValueError("planar patch side must be divisible by voxel size")
        shape = (cells, cells, cells)
        for policy in section["edge_policies"]:
            for alpha in section["alpha_values"]:
                group_records = []
                for family, normal in normals.items():
                    salt = central_plane_mask(shape, normal)
                    coordinates = phase2a3b1.candidate_coordinates_from_mask(salt)
                    records = crystallographic_candidate_diagnostics(
                        salt,
                        coordinates,
                        voxel_size_m=voxel_m,
                        fitting_radius_m=fitting_m,
                        grain_to_lab=np.eye(3),
                        alpha=float(alpha),
                        edge_policy=str(policy),
                        quality_minimum=float(config["quality_minimum"]),
                        near_tolerance_deg=float(config["near_family_tolerance_deg"]),
                    )
                    area = exact_central_section_area_m2(family, side_m)
                    weights = np.asarray([record.anisotropy_weight for record in records], dtype=float)
                    digital_area = float(sum(record.represented_area_m2 for record in records))
                    group_records.append(
                        {
                            "orientation_family": family,
                            "candidate_count": len(records),
                            "physical_area_m2": area,
                            "digital_face_area_m2": digital_area,
                            "mean_weight": float(np.mean(weights)),
                            "weight_sum": float(np.sum(weights)),
                            "anisotropy_eligible_fraction": float(
                                np.mean([record.anisotropy_eligible for record in records])
                            ),
                            "face_like_fraction": float(np.mean([record.site_class == "face_like" for record in records])),
                            "edge_like_fraction": float(np.mean([record.site_class == "edge_like" for record in records])),
                            "corner_like_fraction": float(np.mean([record.site_class == "corner_like" for record in records])),
                        }
                    )
                counts = np.asarray([row["candidate_count"] for row in group_records], dtype=float)
                areas = np.asarray([row["physical_area_m2"] for row in group_records], dtype=float)
                weights = np.asarray([row["mean_weight"] for row in group_records], dtype=float)
                for allocation in ("site_weighted", "area_weighted"):
                    probabilities = group_allocation(counts, areas, weights, allocation=allocation)
                    expected_total_events = total_volume / voxel_m**3
                    displacements = total_volume * probabilities / areas
                    reference = float(displacements[0])
                    for base, probability, displacement in zip(
                        group_records, probabilities, displacements, strict=True
                    ):
                        rows.append(
                            {
                                "voxel_size_um": voxel_um,
                                "cells_per_patch_side": cells,
                                "edge_policy": policy,
                                "alpha": alpha,
                                "allocation": allocation,
                                **base,
                                "candidate_count_per_physical_area_m_minus2": float(base["candidate_count"]) / float(base["physical_area_m2"]),
                                "digital_to_physical_area_ratio": float(base["digital_face_area_m2"]) / float(base["physical_area_m2"]),
                                "group_probability": float(probability),
                                "expected_events": float(expected_total_events * probability),
                                "selected_events_per_physical_area_m_minus2": float(expected_total_events * probability / float(base["physical_area_m2"])),
                                "deposited_volume_per_area_m": float(displacement),
                                "effective_facet_displacement_um": float(displacement * 1.0e6),
                                "effective_velocity_ratio_to_f100": float(displacement / reference),
                                "same_total_deposited_volume_m3": total_volume,
                                "implemented_in_ca": allocation == "site_weighted",
                            }
                        )
    return rows


def stage_lookup(rows: Sequence[dict[str, str]]) -> dict[tuple[str, float], dict[str, str]]:
    return {(row["case_id"], float(row["stage_fraction"])): row for row in rows}


def facet_fraction_vector(row: dict[str, str]) -> np.ndarray:
    return np.asarray([float(row[f"surface_area_fraction_{FAMILY_COLUMNS[family]}"]) for family in FAMILY_NAMES])


def extent_vector(row: dict[str, str]) -> np.ndarray:
    return np.asarray([float(row[f"growth_extent_{axis}_total_m"]) for axis in ("x", "y", "z")])


def comparison_row(
    *,
    source: str,
    first_id: str,
    second_id: str,
    stage: float,
    event_rows: Sequence[dict[str, str]],
    stages: dict[tuple[str, float], dict[str, str]],
    voxel_m: float,
) -> dict[str, object]:
    first_points = phase2a3b1.event_points(event_rows, first_id, stage)
    second_points = phase2a3b1.event_points(event_rows, second_id, stage)
    first_family = phase2a3b1.family_fraction_vector(event_rows, first_id, stage)
    second_family = phase2a3b1.family_fraction_vector(event_rows, second_id, stage)
    first_stage = stages[(first_id, stage)]
    second_stage = stages[(second_id, stage)]
    first_extent = extent_vector(first_stage)
    second_extent = extent_vector(second_stage)
    return {
        "effect_source": source,
        "first_case_id": first_id,
        "second_case_id": second_id,
        "stage_fraction": stage,
        "growth_only_chamfer_voxels": symmetric_chamfer_distance(first_points, second_points) / voxel_m,
        "event_family_fraction_l1": float(np.sum(np.abs(first_family - second_family))),
        "facet_area_fraction_l1": float(np.sum(np.abs(facet_fraction_vector(first_stage) - facet_fraction_vector(second_stage)))),
        "directional_extent_relative_difference": float(
            np.linalg.norm(first_extent - second_extent) / max(np.linalg.norm(first_extent), voxel_m)
        ),
    }


def existing_phase2a3b1_audit(config: dict[str, object]) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    source = REPOSITORY_ROOT / str(config["existing_phase2a3b1_output"])
    events = read_csv(source / "growth_events.csv")
    stage_rows = read_csv(source / "growth_stages.csv")
    stages = stage_lookup(stage_rows)
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    rotations = (0, 1)
    policies = ("smooth_average", "isotropic_ambiguous")
    alphas = (-float(config["alpha_magnitude"]), 0.0, float(config["alpha_magnitude"]))
    seeds = (101, 202, 303)
    fractions = (0.01, 0.05, 0.10)
    comparisons: list[dict[str, object]] = []
    for stage in fractions:
        for rotation in rotations:
            for policy in policies:
                for seed in seeds:
                    baseline = case_id(0.0, policy, seed, rotation)
                    for alpha in (alphas[0], alphas[2]):
                        comparisons.append(
                            comparison_row(
                                source="alpha_response",
                                first_id=case_id(alpha, policy, seed, rotation),
                                second_id=baseline,
                                stage=stage,
                                event_rows=events,
                                stages=stages,
                                voxel_m=voxel_m,
                            )
                        )
            for alpha in alphas:
                for seed in seeds:
                    comparisons.append(
                        comparison_row(
                            source="edge_policy",
                            first_id=case_id(alpha, policies[0], seed, rotation),
                            second_id=case_id(alpha, policies[1], seed, rotation),
                            stage=stage,
                            event_rows=events,
                            stages=stages,
                            voxel_m=voxel_m,
                        )
                    )
        for policy in policies:
            for alpha in alphas:
                for seed in seeds:
                    comparisons.append(
                        comparison_row(
                            source="rotation_grid",
                            first_id=case_id(alpha, policy, seed, 0),
                            second_id=case_id(alpha, policy, seed, 1),
                            stage=stage,
                            event_rows=events,
                            stages=stages,
                            voxel_m=voxel_m,
                        )
                    )
        for rotation in rotations:
            for policy in policies:
                for first_seed, second_seed in ((101, 202), (101, 303), (202, 303)):
                    comparisons.append(
                        comparison_row(
                            source="alpha_zero_seed",
                            first_id=case_id(0.0, policy, first_seed, rotation),
                            second_id=case_id(0.0, policy, second_seed, rotation),
                            stage=stage,
                            event_rows=events,
                            stages=stages,
                            voxel_m=voxel_m,
                        )
                    )
    summary_rows: list[dict[str, object]] = []
    metrics = (
        "growth_only_chamfer_voxels",
        "event_family_fraction_l1",
        "facet_area_fraction_l1",
        "directional_extent_relative_difference",
    )
    for source_name in sorted({str(row["effect_source"]) for row in comparisons}):
        for stage in fractions:
            selected = [row for row in comparisons if row["effect_source"] == source_name and np.isclose(float(row["stage_fraction"]), stage)]
            for metric in metrics:
                values = np.asarray([float(row[metric]) for row in selected])
                summary_rows.append(
                    {
                        "effect_source": source_name,
                        "stage_fraction": stage,
                        "metric": metric,
                        "n": values.size,
                        "minimum": float(np.min(values)),
                        "q25": float(np.quantile(values, 0.25)),
                        "median": float(np.median(values)),
                        "q75": float(np.quantile(values, 0.75)),
                        "maximum": float(np.max(values)),
                    }
                )
    chamfer = np.asarray([float(row["growth_only_chamfer_voxels"]) for row in comparisons])
    family = np.asarray([float(row["event_family_fraction_l1"]) for row in comparisons])
    area = np.asarray([float(row["facet_area_fraction_l1"]) for row in comparisons])
    associations = {
        "n_comparisons": len(comparisons),
        "spearman_chamfer_vs_event_family_l1": float(spearmanr(chamfer, family).statistic),
        "spearman_chamfer_vs_facet_area_l1": float(spearmanr(chamfer, area).statistic),
        "interpretation": "Chamfer is complementary, not sufficient: it has no sign and does not directly measure facet-normal advance.",
    }
    return comparisons, summary_rows, associations


def make_growth_case(config: dict[str, object], section: dict[str, object], *, alpha: float, policy: str, seed: int) -> dict[str, object]:
    return {
        "case_id": case_id(alpha, policy, seed),
        "rotation_index": 0,
        "rotation_xyz_deg": section["rotation_xyz_deg"],
        "edge_policy": policy,
        "alpha": alpha,
        "rng_seed": seed,
        "shape": section["shape"],
        "stage_fractions": section["stage_fractions"],
        "stage_events": section["stage_events"],
        "voxel_size_um": config["voxel_size_um"],
        "grain_size_um": config["grain_size_um"],
        "fitting_radius_um": float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"]),
        "quality_minimum": config["quality_minimum"],
        "near_family_tolerance_deg": config["near_family_tolerance_deg"],
        "solid_density_kg_m3": section["solid_density_kg_m3"],
        "retain_stage_masks": True,
    }


def boundary_voxel_count(mask: np.ndarray) -> int:
    boundary = np.zeros(mask.shape, dtype=bool)
    for axis in range(3):
        first = [slice(None)] * 3
        last = [slice(None)] * 3
        first[axis] = 0
        last[axis] = -1
        boundary[tuple(first)] = True
        boundary[tuple(last)] = True
    return int(np.count_nonzero(mask & boundary))


def exploratory_stage(config: dict[str, object], output: Path) -> dict[str, object]:
    resolution_rows = analytical_resolution(config)
    write_csv(output / "analytical_growth_resolution.csv", resolution_rows)
    plane_rows = planar_probability_velocity_audit(config)
    write_csv(output / "planar_probability_velocity_audit.csv", plane_rows)
    comparisons, summaries, associations = existing_phase2a3b1_audit(config)
    write_csv(output / "existing_phase2a3b1_comparisons.csv", comparisons)
    write_csv(output / "existing_phase2a3b1_distribution_summary.csv", summaries)
    write_json(output / "existing_metric_association.json", associations)

    section = config["feasibility_case"]
    case = make_growth_case(
        config,
        section,
        alpha=float(section["alpha"]),
        policy=str(section["edge_policy"]),
        seed=int(section["rng_seed"]),
    )
    result = phase2a3b1.locked_growth_worker(case)
    write_csv(output / "feasibility_growth_stages.csv", result["stage_rows"])
    write_csv(output / "feasibility_growth_events.csv", result["event_rows"])
    np.savez_compressed(
        output / "feasibility_stage_masks.npz",
        **{f"stage_{index}": mask for index, mask in enumerate(result["stage_masks"])},
    )
    last = result["stage_rows"][-1]
    voxel_mass = float(section["solid_density_kg_m3"]) * (float(config["voxel_size_um"]) * 1.0e-6) ** 3
    cache_pass = all(
        float(row["maximum_weight_error"]) <= 1.0e-12
        and int(row["site_class_mismatch_count"]) == 0
        and int(row["anisotropy_eligibility_mismatch_count"]) == 0
        for row in result["stage_rows"]
    )
    mass_pass = all(
        abs(float(row["mass_closure_error_kg"]))
        <= closure_tolerance_kg(float(row["cumulative_input_mass_kg"]), voxel_mass)
        for row in result["stage_rows"]
    )
    boundary_count = boundary_voxel_count(result["final_mask"])
    gates = {
        "realized_all_target_events": int(last["target_new_voxels"]) == int(section["stage_events"][-1]),
        "domain_boundary_clear": boundary_count == 0,
        "domain_boundary_salt_voxels": boundary_count,
        "periodic_full_cache_exact": cache_pass,
        "mass_closure": mass_pass,
        "runtime_within_budget": float(last["runtime_cumulative_s"]) <= float(section["maximum_runtime_s"]),
        "final_runtime_s": float(last["runtime_cumulative_s"]),
        "controlled_campaign_maximum_workers": 6,
    }
    gates["overall_pass"] = all(value for key, value in gates.items() if isinstance(value, bool))
    gates["decision"] = "PROCEED_TO_CONTROLLED_PROTOCOL" if gates["overall_pass"] else "STOP"
    write_json(output / "gate_summary.json", gates)
    plot_exploratory(output, resolution_rows, summaries, plane_rows, result["stage_rows"])
    return gates


def plot_exploratory(output, resolution_rows, summaries, plane_rows, feasibility_rows) -> None:
    mpl.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.4), layout="constrained")
    fractions = np.asarray([float(row["added_volume_fraction"]) for row in resolution_rows])
    displacements = np.asarray([float(row["facet_displacement_voxels"]) for row in resolution_rows])
    axes[0, 0].plot(100.0 * fractions, displacements, marker="o")
    axes[0, 0].axhline(1.0, color="black", linestyle="--", linewidth=1)
    axes[0, 0].set(xlabel="Added solid volume (%)", ylabel="Uniform-face displacement (voxels)")
    final_summaries = [row for row in summaries if np.isclose(float(row["stage_fraction"]), 0.10) and row["metric"] == "growth_only_chamfer_voxels"]
    axes[0, 1].bar([str(row["effect_source"]).replace("_", "\n") for row in final_summaries], [float(row["median"]) for row in final_summaries])
    axes[0, 1].set(ylabel="Median growth-only Chamfer (voxels)", title="Existing 10% campaign")
    selected_planes = [row for row in plane_rows if float(row["voxel_size_um"]) == 10.0 and row["edge_policy"] == "smooth_average" and float(row["alpha"]) > 0.0]
    for allocation, marker in (("site_weighted", "o"), ("area_weighted", "s")):
        chosen = [row for row in selected_planes if row["allocation"] == allocation]
        axes[1, 0].plot([row["orientation_family"] for row in chosen], [float(row["effective_velocity_ratio_to_f100"]) for row in chosen], marker=marker, label=allocation)
    axes[1, 0].set(ylabel="Effective displacement ratio", xlabel="Frozen plane family")
    axes[1, 0].legend(frameon=False)
    axes[1, 1].plot([100.0 * float(row["stage_fraction"]) for row in feasibility_rows], [float(row["runtime_cumulative_s"]) for row in feasibility_rows], marker="o")
    axes[1, 1].set(xlabel="Added solid volume (%)", ylabel="Cumulative runtime (s)", title="Feasibility case")
    for axis in axes.ravel():
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        path = output / f"exploratory_diagnosis.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white", metadata={"Date": None} if suffix == "svg" else None)
        if suffix == "svg":
            normalize_svg(path)
    plt.close(fig)


def controlled_cases(config: dict[str, object]) -> list[dict[str, object]]:
    section = config["controlled_campaign"]
    cases = []
    for seed in section["rng_seeds"]:
        cases.append(make_growth_case(config, section, alpha=0.0, policy="smooth_average", seed=int(seed)))
        for policy in config["edge_policies"]:
            for alpha in (-float(config["alpha_magnitude"]), float(config["alpha_magnitude"])):
                cases.append(make_growth_case(config, section, alpha=alpha, policy=str(policy), seed=int(seed)))
    return cases


def initial_surface_summary(config: dict[str, object], section: dict[str, object]) -> dict[str, float]:
    rotation = phase2a3b1.rotation_from_angles(section["rotation_xyz_deg"])
    state = phase2a3b1.continuous_cube_state(
        shape=tuple(section["shape"]),
        substrate_layers=1,
        grain_size_m=float(config["grain_size_um"]) * 1.0e-6,
        voxel_size_m=float(config["voxel_size_um"]) * 1.0e-6,
        rotation=rotation,
    )
    coordinates = phase2a3b1.candidate_coordinates_from_material(state.material)
    records = crystallographic_candidate_diagnostics(
        state.material == Material.SALT,
        coordinates,
        voxel_size_m=float(config["voxel_size_um"]) * 1.0e-6,
        fitting_radius_m=float(config["fitting_radius_over_grain_size"]) * float(config["grain_size_um"]) * 1.0e-6,
        grain_to_lab=rotation,
        alpha=0.0,
        edge_policy="smooth_average",
        quality_minimum=float(config["quality_minimum"]),
        near_tolerance_deg=float(config["near_family_tolerance_deg"]),
    )
    return phase2a3b1.family_surface_summary(records)


def cumulative_family_counts(event_rows: Sequence[dict[str, object]], stage: float) -> dict[str, int]:
    selected = [row for row in event_rows if float(row["stage_fraction"]) <= stage]
    return {family: sum(row["orientation_family"] == family for row in selected) for family in FAMILY_NAMES}


def controlled_diagnostics(config, results) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    section = config["controlled_campaign"]
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    rotation = phase2a3b1.rotation_from_angles(section["rotation_xyz_deg"])
    initial_state = phase2a3b1.continuous_cube_state(
        shape=tuple(section["shape"]),
        substrate_layers=1,
        grain_size_m=float(config["grain_size_um"]) * 1.0e-6,
        voxel_size_m=voxel_m,
        rotation=rotation,
    )
    initial_surface = initial_surface_summary(config, section)
    directions = ideal_cubic_direction_families()
    diagnostic_rows = []
    case_summary_rows = []
    for result in results:
        stage_map = {float(row["stage_fraction"]): row for row in result["stage_rows"]}
        for stage, mask in zip(section["stage_fractions"], result["stage_masks"], strict=True):
            row = stage_map[float(stage)]
            counts = cumulative_family_counts(result["event_rows"], float(stage))
            for family in FAMILY_NAMES:
                column = FAMILY_COLUMNS[family]
                initial_area = float(initial_surface[f"surface_area_{column}_m2"])
                current_area = float(row[f"surface_area_{column}_m2"])
                mean_area = 0.5 * (initial_area + current_area)
                if family == "other":
                    advances = np.asarray([np.nan])
                else:
                    advances = support_advances_m(
                        initial_state.seed_mask,
                        mask,
                        voxel_size_m=voxel_m,
                        grain_to_lab=rotation,
                        directions_crystal=directions[family],
                    )
                diagnostic_rows.append(
                    {
                        "case_id": result["case_id"],
                        "edge_policy": row["edge_policy"],
                        "alpha": row["alpha"],
                        "rng_seed": row["rng_seed"],
                        "stage_fraction": stage,
                        "orientation_family": family,
                        "cumulative_events": counts[family],
                        "initial_surface_area_m2": initial_area,
                        "current_surface_area_m2": current_area,
                        "surface_area_fraction": row[f"surface_area_fraction_{column}"],
                        "surface_area_fraction_change": float(row[f"surface_area_fraction_{column}"]) - float(initial_surface[f"surface_area_fraction_{column}"]),
                        "equivalent_facet_displacement_um": equivalent_displacement_m(counts[family], voxel_m, mean_area) * 1.0e6 if mean_area > 0.0 else float("nan"),
                        "support_advance_mean_voxels": float(np.nanmean(advances) / voxel_m) if family != "other" else float("nan"),
                        "support_advance_median_voxels": float(np.nanmedian(advances) / voxel_m) if family != "other" else float("nan"),
                        "support_advance_max_voxels": float(np.nanmax(advances) / voxel_m) if family != "other" else float("nan"),
                    }
                )
            case_summary_rows.append(
                {
                    **row,
                    "domain_boundary_salt_voxels": boundary_voxel_count(mask),
                }
            )
    return diagnostic_rows, case_summary_rows


def result_lookup(results) -> dict[str, dict[str, object]]:
    return {str(result["case_id"]): result for result in results}


def controlled_effects(config, results, diagnostics) -> tuple[list[dict[str, object]], dict[str, object]]:
    section = config["controlled_campaign"]
    voxel_m = float(config["voxel_size_um"]) * 1.0e-6
    lookup = result_lookup(results)
    diag = {(row["case_id"], float(row["stage_fraction"]), row["orientation_family"]): row for row in diagnostics}
    rows = []
    for policy in config["edge_policies"]:
        for alpha in (-float(config["alpha_magnitude"]), float(config["alpha_magnitude"])):
            for seed in section["rng_seeds"]:
                candidate = lookup[case_id(alpha, policy, int(seed))]
                baseline = lookup[case_id(0.0, "smooth_average", int(seed))]
                for stage in section["stage_fractions"]:
                    candidate_points = phase2a3b1.event_points(candidate["event_rows"], candidate["case_id"], float(stage))
                    baseline_points = phase2a3b1.event_points(baseline["event_rows"], baseline["case_id"], float(stage))
                    candidate_family = phase2a3b1.family_fraction_vector(candidate["event_rows"], candidate["case_id"], float(stage))
                    baseline_family = phase2a3b1.family_fraction_vector(baseline["event_rows"], baseline["case_id"], float(stage))
                    candidate_equiv = np.asarray([float(diag[(candidate["case_id"], float(stage), family)]["equivalent_facet_displacement_um"]) for family in ("<100>", "<110>", "<111>")])
                    baseline_equiv = np.asarray([float(diag[(baseline["case_id"], float(stage), family)]["equivalent_facet_displacement_um"]) for family in ("<100>", "<110>", "<111>")])
                    candidate_support = np.asarray([float(diag[(candidate["case_id"], float(stage), family)]["support_advance_mean_voxels"]) for family in ("<100>", "<110>", "<111>")])
                    baseline_support = np.asarray([float(diag[(baseline["case_id"], float(stage), family)]["support_advance_mean_voxels"]) for family in ("<100>", "<110>", "<111>")])
                    rows.append(
                        {
                            "edge_policy": policy,
                            "alpha": alpha,
                            "rng_seed": seed,
                            "stage_fraction": stage,
                            "growth_only_chamfer_voxels": symmetric_chamfer_distance(candidate_points, baseline_points) / voxel_m,
                            "event_family_fraction_l1": float(np.sum(np.abs(candidate_family - baseline_family))),
                            "equivalent_displacement_l1_um": float(np.sum(np.abs(candidate_equiv - baseline_equiv))),
                            "support_advance_l1_voxels": float(np.sum(np.abs(candidate_support - baseline_support))),
                            "paired_f100_minus_f111_equivalent_contrast_um": float((candidate_equiv[0] - candidate_equiv[2]) - (baseline_equiv[0] - baseline_equiv[2])),
                            "paired_f100_minus_f111_support_contrast_voxels": float((candidate_support[0] - candidate_support[2]) - (baseline_support[0] - baseline_support[2])),
                        }
                    )
    final_stage = max(float(value) for value in section["stage_fractions"])
    stochastic = []
    seeds = [int(value) for value in section["rng_seeds"]]
    for index, first_seed in enumerate(seeds):
        for second_seed in seeds[index + 1 :]:
            first = lookup[case_id(0.0, "smooth_average", first_seed)]
            second = lookup[case_id(0.0, "smooth_average", second_seed)]
            first_points = phase2a3b1.event_points(first["event_rows"], first["case_id"], final_stage)
            second_points = phase2a3b1.event_points(second["event_rows"], second["case_id"], final_stage)
            stochastic.append(symmetric_chamfer_distance(first_points, second_points) / voxel_m)
    baseline = float(np.median(stochastic))
    conditions = {}
    all_signal_pass = True
    all_sign_pass = True
    for policy in config["edge_policies"]:
        for alpha in (-float(config["alpha_magnitude"]), float(config["alpha_magnitude"])):
            selected = [row for row in rows if row["edge_policy"] == policy and np.isclose(float(row["alpha"]), alpha) and np.isclose(float(row["stage_fraction"]), final_stage)]
            chamfer = np.asarray([float(row["growth_only_chamfer_voxels"]) for row in selected])
            contrast = np.asarray([float(row["paired_f100_minus_f111_equivalent_contrast_um"]) for row in selected])
            expected_sign = 1.0 if alpha > 0.0 else -1.0
            sign_fraction = float(np.mean(expected_sign * contrast > 0.0))
            ratio = float(np.median(chamfer) / baseline)
            key = f"{policy}|{alpha:+.12g}"
            conditions[key] = {
                "median_chamfer_voxels": float(np.median(chamfer)),
                "chamfer_range_voxels": [float(np.min(chamfer)), float(np.max(chamfer))],
                "alpha_zero_seed_baseline_median_voxels": baseline,
                "signal_to_baseline_ratio": ratio,
                "equivalent_displacement_contrast_um": [float(value) for value in contrast],
                "expected_sign_fraction": sign_fraction,
            }
            all_signal_pass &= ratio >= float(config["acceptance"]["retained_phase2a3b1_signal_ratio_min"])
            all_sign_pass &= sign_fraction >= float(config["acceptance"]["minimum_seed_sign_fraction"])
    return rows, {
        "alpha_zero_seed_chamfer_values_voxels": stochastic,
        "alpha_zero_seed_chamfer_median_voxels": baseline,
        "conditions": conditions,
        "retained_morphology_signal_gate": all_signal_pass,
        "facet_contrast_sign_consistency": all_sign_pass,
    }


def controlled_stage(config: dict[str, object], output: Path) -> dict[str, object]:
    cases = controlled_cases(config)
    workers = min(int(config["controlled_campaign"]["maximum_workers"]), len(cases))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(phase2a3b1.locked_growth_worker, cases))
    stages = [row for result in results for row in result["stage_rows"]]
    events = [row for result in results for row in result["event_rows"]]
    write_csv(output / "controlled_growth_stages.csv", stages)
    write_csv(output / "controlled_growth_events.csv", events)
    diagnostics, case_summaries = controlled_diagnostics(config, results)
    write_csv(output / "controlled_facet_diagnostics.csv", diagnostics)
    write_csv(output / "controlled_case_summaries.csv", case_summaries)
    effects, response = controlled_effects(config, results, diagnostics)
    write_csv(output / "controlled_paired_effects.csv", effects)
    write_json(output / "controlled_response_summary.json", response)
    np.savez_compressed(
        output / "controlled_final_masks.npz",
        **{result["case_id"]: result["final_mask"] for result in results},
    )
    voxel_mass = float(config["controlled_campaign"]["solid_density_kg_m3"]) * (float(config["voxel_size_um"]) * 1.0e-6) ** 3
    cache_pass = all(float(row["maximum_weight_error"]) <= 1.0e-12 and int(row["site_class_mismatch_count"]) == 0 and int(row["anisotropy_eligibility_mismatch_count"]) == 0 for row in stages)
    mass_pass = all(abs(float(row["mass_closure_error_kg"])) <= closure_tolerance_kg(float(row["cumulative_input_mass_kg"]), voxel_mass) for row in stages)
    boundary_pass = all(int(row["domain_boundary_salt_voxels"]) == 0 for row in case_summaries)
    gates = {
        "preserved_phase2a3b1_no_go": True,
        "periodic_full_cache_exact": cache_pass,
        "mass_closure": mass_pass,
        "domain_boundary_clear": boundary_pass,
        "retained_phase2a3b1_morphology_signal_gate": response["retained_morphology_signal_gate"],
        "facet_contrast_sign_consistency": response["facet_contrast_sign_consistency"],
        "site_weighted_law_has_unique_physical_velocity_mapping": False,
    }
    gates["current_site_weighted_model_go"] = all(
        gates[key]
        for key in (
            "periodic_full_cache_exact",
            "mass_closure",
            "domain_boundary_clear",
            "retained_phase2a3b1_morphology_signal_gate",
            "facet_contrast_sign_consistency",
            "site_weighted_law_has_unique_physical_velocity_mapping",
        )
    )
    gates["decision"] = "GO" if gates["current_site_weighted_model_go"] else "NO-GO: modified kinetic growth law required before porous-salt continuation"
    write_json(output / "gate_summary.json", gates)
    plot_controlled(config, output, diagnostics, effects, response)
    return gates


def plot_controlled(config, output, diagnostics, effects, response) -> None:
    mpl.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.3), layout="constrained")
    final_stage = max(config["controlled_campaign"]["stage_fractions"])
    colors = {"<100>": "#0072B2", "<110>": "#E69F00", "<111>": "#D55E00"}
    selected = [row for row in diagnostics if row["edge_policy"] == "smooth_average" and int(row["rng_seed"]) == 101 and row["orientation_family"] != "other"]
    for family in ("<100>", "<110>", "<111>"):
        rows = [row for row in selected if row["orientation_family"] == family and float(row["alpha"]) > 0.0]
        axes[0].plot([100.0 * float(row["stage_fraction"]) for row in rows], [float(row["equivalent_facet_displacement_um"]) for row in rows], marker="o", color=colors[family], label=family)
    axes[0].set(xlabel="Added solid volume (%)", ylabel="Equivalent facet displacement (um)")
    axes[0].set_title("smooth, positive alpha, seed 101")
    axes[0].legend(frameon=False)
    for policy, marker in (("smooth_average", "o"), ("isotropic_ambiguous", "s")):
        for alpha, linestyle in ((-float(config["alpha_magnitude"]), "--"), (float(config["alpha_magnitude"]), "-")):
            rows = [row for row in effects if row["edge_policy"] == policy and np.isclose(float(row["alpha"]), alpha)]
            medians = []
            for stage in config["controlled_campaign"]["stage_fractions"]:
                values = [float(row["growth_only_chamfer_voxels"]) for row in rows if np.isclose(float(row["stage_fraction"]), float(stage))]
                medians.append(float(np.median(values)))
            policy_label = "smooth" if policy == "smooth_average" else "isotropic"
            axes[1].plot(100.0 * np.asarray(config["controlled_campaign"]["stage_fractions"]), medians, marker=marker, linestyle=linestyle, label=f"{policy_label}, {'+' if alpha > 0 else '-'}")
    axes[1].set(xlabel="Added solid volume (%)", ylabel="Paired Chamfer vs alpha zero (voxels)")
    axes[1].legend(frameon=False, fontsize=6)
    keys = list(response["conditions"])
    ratios = [response["conditions"][key]["signal_to_baseline_ratio"] for key in keys]
    short_labels = [
        ("smooth" if key.startswith("smooth_average") else "isotropic")
        + (" +" if "|+" in key else " -")
        for key in keys
    ]
    axes[2].bar(np.arange(len(keys)), ratios)
    axes[2].axhline(float(config["acceptance"]["retained_phase2a3b1_signal_ratio_min"]), color="black", linestyle="--", linewidth=1)
    axes[2].set(xticks=np.arange(len(keys)), xticklabels=short_labels, ylabel="Final signal / stochastic baseline")
    for axis in axes:
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        path = output / f"controlled_deep_growth.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white", metadata={"Date": None} if suffix == "svg" else None)
        if suffix == "svg":
            normalize_svg(path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("exploratory", "controlled"), required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.mkdir(parents=True)
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if config["stage"] != args.stage:
        raise ValueError("config stage does not match requested stage")
    started = perf_counter()
    gates = exploratory_stage(config, args.output) if args.stage == "exploratory" else controlled_stage(config, args.output)
    manifest = {
        "stage": args.stage,
        "base_commit": config["base_commit"],
        "config": args.config.resolve().relative_to(REPOSITORY_ROOT.resolve()).as_posix(),
        "config_sha256": sha256(args.config),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "runtime_s": perf_counter() - started,
        "maximum_workers": 1 if args.stage == "exploratory" else config["controlled_campaign"]["maximum_workers"],
        "phase2a3b1_no_go_preserved": True,
        "physical_validation_claimed": False,
    }
    write_json(args.output / "manifest.json", manifest)
    write_checksums(args.output)
    print(json.dumps(gates, indent=2))


if __name__ == "__main__":
    main()
