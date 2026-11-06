"""Run the Phase 2A.2 morphology feasibility and grid-robustness audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import matplotlib
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca" / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase1 import EvaporationParameters, compute_mass_budget  # noqa: E402
from salt_phase2 import (  # noqa: E402
    AttachmentWeights,
    CaseConfiguration,
    Material,
    analyze_morphology,
    column_envelope,
    initialize_physical_state,
    run_case,
)
from salt_phase2.plotting_phase2a2 import (  # noqa: E402
    ClassifiedMorphologyView,
    plot_envelope_sensitivity,
    plot_grid_sensitivity,
    render_classified_morphologies,
)


BASE_COMMIT = "098f5bb8cec41ab28502fa3d6b1126057a14901a"
DOMAIN_M = (200.0e-6, 200.0e-6, 180.0e-6)
SUBSTRATE_THICKNESS_M = 20.0e-6
SEED_LOWER_M = (80.0e-6, 80.0e-6, 20.0e-6)
SEED_SIZE_M = 20.0e-6
DURATION_S = 600.0
TIMESTEP_S = 60.0
SOLID_DENSITY_KG_M3 = 2160.0
GRID_SEEDS = (101, 202, 303)
VOXEL_SIZES_M = (5.0e-6, 10.0e-6, 20.0e-6)
WEIGHT_SETS = {
    "uniform": AttachmentWeights(1.0, 1.0, 1.0, 1.0),
    "strong_high_coord": AttachmentWeights(1.0, 4.0, 16.0, 64.0),
}
PHASE2A1_SEEDS = ((6, 6, 1), (6, 13, 1), (13, 6, 1), (13, 13, 1))
TARGET_NEW_VOLUME_M3 = 160 * (10.0e-6) ** 3


def _integer_ratio(numerator: float, denominator: float, name: str) -> int:
    value = numerator / denominator
    rounded = int(round(value))
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-10):
        raise ValueError(f"{name} is not an integer number of voxels")
    return rounded


def _grid_geometry(
    voxel_size_m: float,
) -> tuple[tuple[int, int, int], int, tuple[tuple[int, int, int], ...]]:
    shape = tuple(
        _integer_ratio(length, voxel_size_m, "domain length")
        for length in DOMAIN_M
    )
    substrate_layers = _integer_ratio(
        SUBSTRATE_THICKNESS_M, voxel_size_m, "substrate thickness"
    )
    lower = tuple(
        _integer_ratio(value, voxel_size_m, "seed lower coordinate")
        for value in SEED_LOWER_M
    )
    side = _integer_ratio(SEED_SIZE_M, voxel_size_m, "seed side")
    coordinates = tuple(
        (x, y, z)
        for x in range(lower[0], lower[0] + side)
        for y in range(lower[1], lower[1] + side)
        for z in range(lower[2], lower[2] + side)
    )
    return shape, substrate_layers, coordinates


def _evaporation_area_scale() -> float:
    full_area = DOMAIN_M[0] * DOMAIN_M[1]
    parameters = EvaporationParameters(
        temperature_k=298.15,
        water_activity=0.753,
        relative_humidity=0.30,
        mass_transfer_coefficient_m_s=0.005,
        evaporation_area_m2=full_area,
        timestep_s=TIMESTEP_S,
        salt_mass_fraction=0.264,
        solid_density_kg_m3=SOLID_DENSITY_KG_M3,
    )
    unscaled_step_mass = compute_mass_budget(parameters).precipitated_salt_mass_kg
    target_mass = SOLID_DENSITY_KG_M3 * TARGET_NEW_VOLUME_M3
    return target_mass / ((DURATION_S / TIMESTEP_S) * unscaled_step_mass)


def build_grid_configurations() -> list[CaseConfiguration]:
    area_scale = _evaporation_area_scale()
    configurations = []
    for rule_name, weights in WEIGHT_SETS.items():
        for voxel_size_m in VOXEL_SIZES_M:
            shape, substrate_layers, seed_coordinates = _grid_geometry(voxel_size_m)
            for rng_seed in GRID_SEEDS:
                voxel_size_um = int(round(voxel_size_m * 1.0e6))
                configurations.append(
                    CaseConfiguration(
                        case_id=(
                            f"grid_{rule_name}_{voxel_size_um}um_seed{rng_seed}"
                        ),
                        rule_name=rule_name,
                        shape=shape,
                        substrate_layers=substrate_layers,
                        seed_coordinates=seed_coordinates,
                        wetted_seed_indices=None,
                        enforce_wetting=True,
                        growth_mode="stochastic",
                        update_mode="sequential",
                        attachment_weights=weights,
                        rng_seed=rng_seed,
                        timestep_s=TIMESTEP_S,
                        total_duration_s=DURATION_S,
                        evaporation_area_scale=area_scale,
                        voxel_size_m=voxel_size_m,
                        solid_density_kg_m3=SOLID_DENSITY_KG_M3,
                    )
                )
    return configurations


def _write_csv(path: Path, rows) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _state_from_mask(
    salt_mask: np.ndarray,
    *,
    substrate_layers: int,
    seed_coordinates,
):
    state = initialize_physical_state(
        shape=salt_mask.shape,
        substrate_layers=substrate_layers,
        seed_coordinates=seed_coordinates,
    )
    if np.any(salt_mask[:, :, :substrate_layers]):
        raise ValueError("salt mask overlaps the substrate")
    if not np.all(salt_mask[state.seed_mask]):
        raise ValueError("salt mask does not preserve every declared seed")
    state.material[:, :, substrate_layers:] = Material.EMPTY
    state.material[salt_mask] = Material.SALT
    state.validate()
    return state


def _physical_extent_um(mask: np.ndarray, voxel_size_m: float) -> tuple[float, float]:
    coordinates = np.argwhere(mask)
    if not coordinates.size:
        return 0.0, 0.0
    x_extent = (int(np.ptp(coordinates[:, 0])) + 1) * voxel_size_m * 1.0e6
    y_extent = (int(np.ptp(coordinates[:, 1])) + 1) * voxel_size_m * 1.0e6
    return x_extent, y_extent


def _window_fraction_rows(result) -> list[dict[str, object]]:
    salt = result.final_salt_mask
    dx = result.configuration.voxel_size_m
    rows = []
    for window_um in (80.0, 120.0, 160.0):
        lower_xy_m = 0.5 * (DOMAIN_M[0] - window_um * 1.0e-6)
        upper_xy_m = lower_xy_m + window_um * 1.0e-6
        lower = _integer_ratio(lower_xy_m, dx, "window lower")
        upper = _integer_ratio(upper_xy_m, dx, "window upper")
        lower_z = result.configuration.substrate_layers
        upper_z = min(
            salt.shape[2],
            lower_z + _integer_ratio(window_um * 1.0e-6, dx, "window height"),
        )
        window = salt[lower:upper, lower:upper, lower_z:upper_z]
        rows.append(
            {
                "case_id": result.configuration.case_id,
                "rule_name": result.configuration.rule_name,
                "rng_seed": result.configuration.rng_seed,
                "voxel_size_um": dx * 1.0e6,
                "window_side_um": window_um,
                "window_voxels": int(window.size),
                "salt_voxels": int(np.count_nonzero(window)),
                "solid_volume_fraction": float(np.mean(window)),
            }
        )
    return rows


def _grid_summary_row(result) -> dict[str, object]:
    summary = result.summary
    analysis = result.final_analysis
    dx = result.configuration.voxel_size_m
    x_extent_um, y_extent_um = _physical_extent_um(result.final_salt_mask, dx)
    target_mass = SOLID_DENSITY_KG_M3 * TARGET_NEW_VOLUME_M3
    return {
        **summary,
        "voxel_size_um": dx * 1.0e6,
        "domain_x_um": result.configuration.shape[0] * dx * 1.0e6,
        "domain_y_um": result.configuration.shape[1] * dx * 1.0e6,
        "domain_z_um": result.configuration.shape[2] * dx * 1.0e6,
        "substrate_thickness_um": result.configuration.substrate_layers * dx * 1.0e6,
        "initial_seed_volume_um3": int(np.count_nonzero(result.seed_mask))
        * (dx * 1.0e6) ** 3,
        "target_supplied_mass_kg": target_mass,
        "supplied_mass_relative_error": (
            float(summary["total_supplied_mass_kg"]) - target_mass
        )
        / target_mass,
        "maximum_height_um": float(summary["maximum_height_m"]) * 1.0e6,
        "footprint_area_um2": float(summary["footprint_area_m2"]) * 1.0e12,
        "x_extent_um": x_extent_um,
        "y_extent_um": y_extent_um,
        "surface_to_volume_ratio_per_mm": (
            float(analysis.surface_to_volume_ratio_m_inv) / 1000.0
        ),
        "characteristic_solid_size_um": (
            float(analysis.characteristic_solid_size_m) * 1.0e6
        ),
        "bbox_void_fraction": analysis.bounding_box.void_volume_fraction,
        "column_void_fraction": analysis.column_envelope.void_volume_fraction,
        "face_shell_void_fraction": analysis.face_shell_envelope.void_volume_fraction,
        "column_external_fraction_of_void": (
            analysis.column_envelope.external_air_connected_void_voxels
            / analysis.column_envelope.void_voxels
            if analysis.column_envelope.void_voxels
            else 0.0
        ),
    }


def _descriptive_statistics(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = (
        "maximum_height_um",
        "footprint_area_um2",
        "x_extent_um",
        "y_extent_um",
        "surface_to_volume_ratio_per_mm",
        "characteristic_solid_size_um",
        "bbox_void_fraction",
        "column_void_fraction",
        "runtime_s",
    )
    output = []
    for rule_name in WEIGHT_SETS:
        for voxel_size_m in VOXEL_SIZES_M:
            voxel_size_um = voxel_size_m * 1.0e6
            selected = [
                row
                for row in rows
                if row["rule_name"] == rule_name
                and math.isclose(float(row["voxel_size_um"]), voxel_size_um)
            ]
            for metric in metrics:
                values = np.array([float(row[metric]) for row in selected])
                output.append(
                    {
                        "rule_name": rule_name,
                        "voxel_size_um": voxel_size_um,
                        "metric": metric,
                        "n_simulation_seeds": len(values),
                        "mean": float(np.mean(values)),
                        "sample_sd": float(np.std(values, ddof=1)),
                        "minimum": float(np.min(values)),
                        "maximum": float(np.max(values)),
                    }
                )
    return output


def _grid_changes(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    metrics = (
        "maximum_height_um",
        "footprint_area_um2",
        "surface_to_volume_ratio_per_mm",
        "characteristic_solid_size_um",
        "bbox_void_fraction",
        "column_void_fraction",
    )
    indexed = {
        (row["rule_name"], int(row["rng_seed"]), float(row["voxel_size_um"])): row
        for row in rows
    }
    output = []
    for rule_name in WEIGHT_SETS:
        for rng_seed in GRID_SEEDS:
            reference = indexed[(rule_name, rng_seed, 5.0)]
            for voxel_size_um in (10.0, 20.0):
                candidate = indexed[(rule_name, rng_seed, voxel_size_um)]
                for metric in metrics:
                    reference_value = float(reference[metric])
                    candidate_value = float(candidate[metric])
                    output.append(
                        {
                            "rule_name": rule_name,
                            "rng_seed": rng_seed,
                            "comparison": f"{voxel_size_um:g}um-minus-5um",
                            "metric": metric,
                            "absolute_difference": candidate_value - reference_value,
                            "relative_difference": (
                                (candidate_value - reference_value) / reference_value
                                if reference_value != 0.0
                                else ""
                            ),
                        }
                    )
    return output


def _read_hash_inventory(path: Path) -> dict[str, str]:
    inventory = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", maxsplit=1)
        inventory[relative] = digest
    return inventory


def _verify_reused_mask(path: Path, inventory: dict[str, str], root: Path) -> str:
    relative = path.relative_to(root).as_posix()
    expected = inventory.get(relative)
    if expected is None:
        raise RuntimeError(f"missing source hash for {relative}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"source hash mismatch for {relative}")
    return actual


def reanalyze_phase2a1_masks() -> tuple[list[dict[str, object]], dict[str, object]]:
    source_root = REPOSITORY_ROOT / "phase2_ca" / "outputs" / "phase2a1"
    inventory = _read_hash_inventory(source_root / "SHA256SUMS.txt")
    rows = []
    views = {}
    source_hashes = {}
    for rule_name in WEIGHT_SETS:
        for rng_seed in (101, 202, 303, 404, 505):
            case_id = f"rule_{rule_name}_seed{rng_seed}"
            path = source_root / "masks" / f"{case_id}_final.npy"
            source_hashes[path.relative_to(source_root).as_posix()] = _verify_reused_mask(
                path, inventory, source_root
            )
            salt = np.load(path, allow_pickle=False).astype(bool, copy=False)
            state = _state_from_mask(
                salt,
                substrate_layers=1,
                seed_coordinates=PHASE2A1_SEEDS,
            )
            analysis = analyze_morphology(state, voxel_size_m=10.0e-6)
            envelopes = {
                "bounding_box": analysis.bounding_box,
                "column": analysis.column_envelope,
                "face_shell": analysis.face_shell_envelope,
            }
            row = {
                "case_id": case_id,
                "rule_name": rule_name,
                "rng_seed": rng_seed,
                "total_salt_voxels": int(np.count_nonzero(salt)),
                "salt_cluster_count": analysis.salt_cluster_count,
                "salt_cluster_sizes": json.dumps(analysis.salt_cluster_sizes),
                "computational_bottom_to_top_percolates": (
                    analysis.computational_box.bottom_to_top_void_percolates
                ),
                "computational_bottom_to_lateral_percolates": (
                    analysis.computational_box.bottom_to_lateral_void_percolates
                ),
                "column_continuous_candidate_path": (
                    analysis.column_envelope.bottom_to_top_void_percolates
                    or analysis.column_envelope.bottom_to_lateral_void_percolates
                ),
                "tested_envelope_pathway_agreement": (
                    analysis.tested_envelope_pathway_agreement
                ),
                "tested_envelope_void_fraction_range": (
                    analysis.tested_envelope_void_fraction_range
                ),
                "intrinsic_porosity_identifiable": (
                    analysis.intrinsic_porosity_identifiable
                ),
                "transport_suitability_established": False,
            }
            for name, envelope in envelopes.items():
                row.update(
                    {
                        f"{name}_void_fraction": envelope.void_volume_fraction,
                        f"{name}_void_components": envelope.void_component_count,
                        f"{name}_external_void_voxels": (
                            envelope.external_air_connected_void_voxels
                        ),
                        f"{name}_enclosed_void_voxels": envelope.enclosed_void_voxels,
                        f"{name}_bottom_to_top_percolates": (
                            envelope.bottom_to_top_void_percolates
                        ),
                        f"{name}_bottom_to_lateral_percolates": (
                            envelope.bottom_to_lateral_void_percolates
                        ),
                    }
                )
            row["column_external_fraction_of_void"] = (
                analysis.column_envelope.external_air_connected_void_voxels
                / analysis.column_envelope.void_voxels
                if analysis.column_envelope.void_voxels
                else 0.0
            )
            rows.append(row)
            if rng_seed == 101:
                views[rule_name] = ClassifiedMorphologyView(
                    label=(
                        "Open uniform, seed 101"
                        if rule_name == "uniform"
                        else "Compact high-coordination, seed 101"
                    ),
                    salt_mask=salt,
                    seed_mask=state.seed_mask.copy(),
                    candidate_pore_mask=column_envelope(state) & ~salt,
                    substrate_layers=1,
                    voxel_size_um=10.0,
                )
    return rows, {"source_hashes": source_hashes, "views": views}


def _hash_outputs(output_dir: Path) -> None:
    records = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            records.append(
                f"{hashlib.sha256(path.read_bytes()).hexdigest()}  "
                f"{path.relative_to(output_dir).as_posix()}"
            )
    (output_dir / "SHA256SUMS.txt").write_text(
        "\n".join(records) + "\n", encoding="utf-8"
    )


def run_campaign(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "masks").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    configurations = build_grid_configurations()
    results = []
    for index, configuration in enumerate(configurations, start=1):
        print(f"[{index}/{len(configurations)}] {configuration.case_id}", flush=True)
        result = run_case(configuration)
        results.append(result)
        np.save(
            output_dir / "masks" / f"{configuration.case_id}_final.npy",
            result.final_salt_mask,
            allow_pickle=False,
        )

    target_mass = SOLID_DENSITY_KG_M3 * TARGET_NEW_VOLUME_M3
    for result in results:
        if not math.isclose(
            float(result.summary["total_supplied_mass_kg"]),
            target_mass,
            rel_tol=2.0e-14,
            abs_tol=0.0,
        ):
            raise RuntimeError("grid cases do not share the target supplied salt mass")
        expected = round(TARGET_NEW_VOLUME_M3 / result.configuration.voxel_size_m**3)
        if int(result.summary["new_salt_voxels"]) != expected:
            raise RuntimeError("whole-voxel realization differs across nested grids")

    grid_rows = [_grid_summary_row(result) for result in results]
    window_rows = [
        row for result in results for row in _window_fraction_rows(result)
    ]
    _write_csv(output_dir / "grid_cases.csv", grid_rows)
    _write_csv(output_dir / "grid_window_fractions.csv", window_rows)
    _write_csv(
        output_dir / "grid_descriptive_statistics.csv",
        _descriptive_statistics(grid_rows),
    )
    _write_csv(output_dir / "grid_changes_relative_to_5um.csv", _grid_changes(grid_rows))
    _write_csv(
        output_dir / "grid_conservation_history.csv",
        [row for result in results for row in result.history],
    )

    phase2a1_rows, reused = reanalyze_phase2a1_masks()
    _write_csv(output_dir / "phase2a1_morphology_reanalysis.csv", phase2a1_rows)
    plot_grid_sensitivity(grid_rows, output_dir / "figures" / "grid_sensitivity.png")
    plot_envelope_sensitivity(
        phase2a1_rows, output_dir / "figures" / "envelope_sensitivity.png"
    )

    representative_results = [
        next(
            result
            for result in results
            if result.configuration.rule_name == "uniform"
            and result.configuration.rng_seed == 101
            and math.isclose(result.configuration.voxel_size_m, voxel_size_m)
        )
        for voxel_size_m in VOXEL_SIZES_M
    ]
    grid_views = []
    for result in representative_results:
        state = _state_from_mask(
            result.final_salt_mask,
            substrate_layers=result.configuration.substrate_layers,
            seed_coordinates=result.configuration.seed_coordinates,
        )
        grid_views.append(
            ClassifiedMorphologyView(
                label=f"Uniform, {result.configuration.voxel_size_m * 1e6:.0f} µm",
                salt_mask=result.final_salt_mask,
                seed_mask=result.seed_mask,
                candidate_pore_mask=column_envelope(state)
                & ~result.final_salt_mask,
                substrate_layers=result.configuration.substrate_layers,
                voxel_size_um=result.configuration.voxel_size_m * 1.0e6,
            )
        )
    slice_records = render_classified_morphologies(
        grid_views,
        output_dir / "figures" / "grid_morphologies.png",
        title="Grid sensitivity: salt, substrate, exterior air, and column-mask pore",
    )
    phase2a1_views = [
        reused["views"]["uniform"],
        reused["views"]["strong_high_coord"],
    ]
    slice_records.extend(
        render_classified_morphologies(
            phase2a1_views,
            output_dir / "figures" / "phase2a1_open_compact_pores.png",
            title="Matched inventory: open versus compact Phase 2A.1 morphology",
        )
    )
    _write_csv(output_dir / "figure_slice_indices.csv", slice_records)

    manifest = {
        "provenance": {
            "base_commit": BASE_COMMIT,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "physical_validation_claimed": False,
            "permeability_calculated": False,
            "phase2b_started": False,
        },
        "grid_design": {
            "physical_domain_m": DOMAIN_M,
            "substrate_thickness_m": SUBSTRATE_THICKNESS_M,
            "seed_lower_m": SEED_LOWER_M,
            "seed_size_m": SEED_SIZE_M,
            "duration_s": DURATION_S,
            "timestep_s": TIMESTEP_S,
            "voxel_sizes_m": VOXEL_SIZES_M,
            "rng_seeds": GRID_SEEDS,
            "target_supplied_salt_mass_kg": target_mass,
            "target_new_salt_volume_m3": TARGET_NEW_VOLUME_M3,
            "effective_evaporation_area_m2": (
                DOMAIN_M[0] * DOMAIN_M[1] * _evaporation_area_scale()
            ),
            "statistics": (
                "raw simulation seeds, mean, sample SD, range, and paired "
                "differences; no inferential p-values"
            ),
        },
        "analysis_masks": {
            "computational_box": "whole domain above substrate",
            "bounding_box": "tight axis-aligned box around salt",
            "column": "occupied x-y columns from substrate to local salt top",
            "face_shell": "salt plus one six-face-neighbor voxel shell",
            "candidate_pore_figure": "empty voxels inside the column mask",
        },
        "reused_phase2a1_source_hashes": reused["source_hashes"],
        "case_ids": [configuration.case_id for configuration in configurations],
        "case_count": len(configurations),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    _hash_outputs(output_dir)
    return {
        "grid_rows": grid_rows,
        "phase2a1_rows": phase2a1_rows,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "phase2a2",
    )
    args = parser.parse_args()
    results = run_campaign(args.output_dir)
    print(
        f"Wrote {len(results['grid_rows'])} grid cases and "
        f"{len(results['phase2a1_rows'])} reused-mask analyses to "
        f"{args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
