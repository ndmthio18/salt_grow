"""Run the controlled Phase 2A.1 stochastic growth validation campaign."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import sys

import matplotlib
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca" / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase2 import (  # noqa: E402
    AttachmentWeights,
    CaseConfiguration,
    mask_difference_fraction,
)
from salt_phase2 import mask_jaccard, run_case  # noqa: E402
from salt_phase2.plotting import plot_rule_metrics, render_case_rows  # noqa: E402


SHAPE = (20, 20, 18)
SUBSTRATE_LAYERS = 1
DURATION_S = 600.0
RULE_SEEDS = (101, 202, 303, 404, 505)
PAIRED_SEEDS = (101, 202, 303)
SINGLE_SEEDS = ((10, 10, 1),)
MULTIPLE_SEEDS = ((6, 6, 1), (6, 13, 1), (13, 6, 1), (13, 13, 1))
NINE_SEEDS = tuple((x, y, 1) for x in (4, 10, 16) for y in (4, 10, 16))

WEIGHT_SETS = {
    "uniform": AttachmentWeights(1.0, 1.0, 1.0, 1.0),
    "strong_high_coord": AttachmentWeights(1.0, 4.0, 16.0, 64.0),
    "mild_high_coord": AttachmentWeights(1.0, 1.5, 2.0, 3.0),
}


def _config(
    case_id: str,
    rule_name: str,
    *,
    seeds=MULTIPLE_SEEDS,
    wetted_seed_indices=None,
    enforce_wetting=True,
    growth_mode="stochastic",
    update_mode="sequential",
    weights=WEIGHT_SETS["uniform"],
    rng_seed=101,
    timestep_s=60.0,
    duration_s=DURATION_S,
    area_scale=1.0,
):
    return CaseConfiguration(
        case_id=case_id,
        rule_name=rule_name,
        shape=SHAPE,
        substrate_layers=SUBSTRATE_LAYERS,
        seed_coordinates=tuple(seeds),
        wetted_seed_indices=wetted_seed_indices,
        enforce_wetting=enforce_wetting,
        growth_mode=growth_mode,
        update_mode=update_mode,
        attachment_weights=weights,
        rng_seed=rng_seed,
        timestep_s=timestep_s,
        total_duration_s=duration_s,
        evaporation_area_scale=area_scale,
    )


def build_configurations():
    rule_configs = [
        _config(
            "rule_legacy_deterministic",
            "legacy_deterministic",
            growth_mode="deterministic",
            update_mode="synchronous",
            rng_seed=0,
        )
    ]
    for rule_name, weights in WEIGHT_SETS.items():
        for seed in RULE_SEEDS:
            rule_configs.append(
                _config(
                    f"rule_{rule_name}_seed{seed}",
                    rule_name,
                    weights=weights,
                    rng_seed=seed,
                )
            )

    timestep_configs = []
    for update_mode in ("synchronous", "sequential"):
        for timestep in (30.0, 60.0, 120.0):
            for seed in PAIRED_SEEDS:
                timestep_configs.append(
                    _config(
                        f"dt_{update_mode}_{int(timestep)}s_seed{seed}",
                        f"dt_{update_mode}",
                        seeds=NINE_SEEDS,
                        update_mode=update_mode,
                        rng_seed=seed,
                        timestep_s=timestep,
                        area_scale=0.5,
                    )
                )

    seed_configs = []
    for seed in PAIRED_SEEDS:
        seed_configs.extend(
            [
                _config(
                    f"nucleation_single_seed{seed}",
                    "single_seed",
                    seeds=SINGLE_SEEDS,
                    rng_seed=seed,
                ),
                _config(
                    f"nucleation_multiple_seed{seed}",
                    "multiple_seed",
                    seeds=MULTIPLE_SEEDS,
                    rng_seed=seed,
                ),
            ]
        )

    dry_seeds = ((4, 4, 1), (16, 16, 1))
    wetting_configs = [
        _config(
            "wetting_enforced_dry_cluster",
            "wetting_enforced",
            seeds=dry_seeds,
            wetted_seed_indices=(0,),
            enforce_wetting=True,
            rng_seed=303,
            duration_s=300.0,
        ),
        _config(
            "wetting_disabled_control",
            "wetting_disabled",
            seeds=dry_seeds,
            wetted_seed_indices=(0,),
            enforce_wetting=False,
            rng_seed=303,
            duration_s=300.0,
        ),
    ]

    capacity_configs = [
        _config(
            "capacity_stress_synchronous",
            "capacity_stress_sync",
            seeds=SINGLE_SEEDS,
            update_mode="synchronous",
            rng_seed=707,
            timestep_s=120.0,
            duration_s=120.0,
        ),
        _config(
            "capacity_stress_sequential",
            "capacity_stress_sequential",
            seeds=SINGLE_SEEDS,
            update_mode="sequential",
            rng_seed=707,
            timestep_s=120.0,
            duration_s=120.0,
        ),
    ]
    return {
        "rules": rule_configs,
        "timesteps": timestep_configs,
        "nucleation": seed_configs,
        "wetting": wetting_configs,
        "capacity": capacity_configs,
    }


def _write_csv(path: Path, rows):
    rows = list(rows)
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _descriptive_statistics(rule_results):
    metrics = (
        "new_salt_voxels",
        "maximum_height_voxels",
        "footprint_voxels",
        "salt_cluster_count",
        "bbox_salt_fraction",
        "column_salt_fraction",
        "bbox_enclosed_void_voxels",
        "column_enclosed_void_voxels",
    )
    rows = []
    for rule in ("legacy_deterministic", *WEIGHT_SETS):
        selected = [
            result for result in rule_results if result.configuration.rule_name == rule
        ]
        for metric in metrics:
            values = np.array([float(result.summary[metric]) for result in selected])
            rows.append(
                {
                    "rule_name": rule,
                    "metric": metric,
                    "n_simulation_seeds": len(values),
                    "mean": float(np.mean(values)),
                    "sample_sd": float(np.std(values, ddof=1)) if len(values) > 1 else "",
                    "minimum": float(np.min(values)),
                    "maximum": float(np.max(values)),
                }
            )
    return rows


def _paired_rule_differences(rule_results):
    indexed = {
        (result.configuration.rule_name, result.configuration.rng_seed): result
        for result in rule_results
    }
    rows = []
    for comparison in ("strong_high_coord", "mild_high_coord"):
        for seed in RULE_SEEDS:
            baseline = indexed[("uniform", seed)]
            candidate = indexed[(comparison, seed)]
            rows.append(
                {
                    "comparison": f"{comparison}-minus-uniform",
                    "rng_seed": seed,
                    "height_difference_voxels": (
                        int(candidate.summary["maximum_height_voxels"])
                        - int(baseline.summary["maximum_height_voxels"])
                    ),
                    "footprint_difference_voxels": (
                        int(candidate.summary["footprint_voxels"])
                        - int(baseline.summary["footprint_voxels"])
                    ),
                    "bbox_salt_fraction_difference": (
                        float(candidate.summary["bbox_salt_fraction"])
                        - float(baseline.summary["bbox_salt_fraction"])
                    ),
                    "column_salt_fraction_difference": (
                        float(candidate.summary["column_salt_fraction"])
                        - float(baseline.summary["column_salt_fraction"])
                    ),
                    "mask_jaccard": mask_jaccard(
                        baseline.final_salt_mask, candidate.final_salt_mask
                    ),
                    "mask_difference_fraction": mask_difference_fraction(
                        baseline.final_salt_mask, candidate.final_salt_mask
                    ),
                }
            )
    return rows


def _timestep_comparisons(results):
    indexed = {
        (
            result.configuration.update_mode,
            int(result.configuration.timestep_s),
            result.configuration.rng_seed,
        ): result
        for result in results
    }
    rows = []
    for update_mode in ("synchronous", "sequential"):
        for seed in PAIRED_SEEDS:
            baseline = indexed[(update_mode, 30, seed)]
            for timestep in (60, 120):
                candidate = indexed[(update_mode, timestep, seed)]
                matched = (
                    baseline.summary["new_salt_voxels"]
                    == candidate.summary["new_salt_voxels"]
                )
                rows.append(
                    {
                        "update_mode": update_mode,
                        "rng_seed": seed,
                        "baseline_timestep_s": 30,
                        "comparison_timestep_s": timestep,
                        "matched_new_salt_inventory": matched,
                        "baseline_new_voxels": baseline.summary["new_salt_voxels"],
                        "comparison_new_voxels": candidate.summary["new_salt_voxels"],
                        "mask_jaccard": (
                            mask_jaccard(
                                baseline.final_salt_mask, candidate.final_salt_mask
                            )
                            if matched
                            else ""
                        ),
                        "mask_difference_fraction": (
                            mask_difference_fraction(
                                baseline.final_salt_mask, candidate.final_salt_mask
                            )
                            if matched
                            else ""
                        ),
                        "height_difference_voxels": (
                            int(candidate.summary["maximum_height_voxels"])
                            - int(baseline.summary["maximum_height_voxels"])
                        ),
                        "footprint_difference_voxels": (
                            int(candidate.summary["footprint_voxels"])
                            - int(baseline.summary["footprint_voxels"])
                        ),
                        "baseline_retained_mass_kg": baseline.summary["retained_mass_kg"],
                        "comparison_retained_mass_kg": candidate.summary["retained_mass_kg"],
                    }
                )
    return rows


def _hash_outputs(output_dir: Path):
    records = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file() or path.name == "SHA256SUMS.txt":
            continue
        records.append(
            (
                path.relative_to(output_dir).as_posix(),
                hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        )
    with (output_dir / "SHA256SUMS.txt").open("w", encoding="utf-8") as handle:
        for name, digest in records:
            handle.write(f"{digest}  {name}\n")


def run_campaign(output_dir: Path):
    groups = build_configurations()
    results_by_group = {}
    for group_name, configurations in groups.items():
        print(f"Running {group_name}: {len(configurations)} cases")
        results_by_group[group_name] = [run_case(config) for config in configurations]
    all_results = [result for values in results_by_group.values() for result in values]

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "cases.csv", [result.summary for result in all_results])
    _write_csv(
        output_dir / "conservation_history.csv",
        [row for result in all_results for row in result.history],
    )
    _write_csv(
        output_dir / "rule_descriptive_statistics.csv",
        _descriptive_statistics(results_by_group["rules"]),
    )
    _write_csv(
        output_dir / "paired_rule_differences.csv",
        _paired_rule_differences(results_by_group["rules"]),
    )
    _write_csv(
        output_dir / "timestep_comparisons.csv",
        _timestep_comparisons(results_by_group["timesteps"]),
    )

    mask_dir = output_dir / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    for result in all_results:
        np.save(mask_dir / f"{result.configuration.case_id}_final.npy", result.final_salt_mask)
    rule_representatives = [
        results_by_group["rules"][0],
        *[
            next(
                result
                for result in results_by_group["rules"]
                if result.configuration.rule_name == rule
                and result.configuration.rng_seed == 101
            )
            for rule in WEIGHT_SETS
        ],
    ]
    for result in rule_representatives:
        for step, mask in result.snapshots.items():
            np.save(
                mask_dir / f"{result.configuration.case_id}_step{step:03d}.npy",
                mask,
            )

    figure_dir = output_dir / "figures"
    slice_records = []
    slice_records.extend(
        render_case_rows(
            rule_representatives,
            ["legacy deterministic", "uniform", "strong", "mild"],
            figure_dir / "rule_morphologies.png",
            title="Rule comparison at matched supplied mass",
        )
    )
    plot_rule_metrics(
        results_by_group["rules"],
        figure_dir / "rule_metric_distributions.png",
        rule_order=(
            "legacy_deterministic",
            "uniform",
            "strong_high_coord",
            "mild_high_coord",
        ),
    )
    timestep_index = {
        (result.configuration.update_mode, int(result.configuration.timestep_s)): result
        for result in results_by_group["timesteps"]
        if result.configuration.rng_seed == 101
    }
    for update_mode in ("synchronous", "sequential"):
        cases = [timestep_index[(update_mode, timestep)] for timestep in (30, 60, 120)]
        slice_records.extend(
            render_case_rows(
                cases,
                [f"{timestep} s" for timestep in (30, 60, 120)],
                figure_dir / f"timestep_{update_mode}.png",
                title=f"{update_mode.capitalize()} growth: paired RNG seed 101",
            )
        )
    nucleation_cases = [
        next(
            result
            for result in results_by_group["nucleation"]
            if result.configuration.rule_name == rule
            and result.configuration.rng_seed == 101
        )
        for rule in ("single_seed", "multiple_seed")
    ]
    slice_records.extend(
        render_case_rows(
            nucleation_cases,
            ["single seed", "four seeds"],
            figure_dir / "nucleation_morphologies.png",
            title="Nucleation topology at matched supplied mass",
        )
    )
    slice_records.extend(
        render_case_rows(
            results_by_group["wetting"],
            ["wetting enforced", "wetting disabled control"],
            figure_dir / "wetting_control.png",
            title="Dry disconnected seed control",
        )
    )
    slice_records.extend(
        render_case_rows(
            results_by_group["capacity"],
            ["synchronous frozen list", "sequential events"],
            figure_dir / "capacity_stress.png",
            title="One-step capacity stress (unmatched realized inventories)",
        )
    )
    _write_csv(output_dir / "figure_slice_indices.csv", slice_records)

    manifest = {
        "provenance": {
            "reviewed_commit": "96effabf9bca5c8825d2a158dc7c59feb0dd56c1",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "physical_validation_claimed": False,
        },
        "common_configuration": {
            "shape": SHAPE,
            "substrate_layers": SUBSTRATE_LAYERS,
            "duration_s": DURATION_S,
            "rule_rng_seeds": RULE_SEEDS,
            "paired_rng_seeds": PAIRED_SEEDS,
            "weight_sets": {name: asdict(weights) for name, weights in WEIGHT_SETS.items()},
            "statistics": (
                "descriptive mean, sample SD, range, paired mask Jaccard and difference; "
                "no inferential p-values because simulation seeds are not experimental replicates"
            ),
            "envelopes": (
                "axis-aligned salt bounding box and footprint columns to local salt top"
            ),
        },
        "case_ids": [result.configuration.case_id for result in all_results],
        "case_count": len(all_results),
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    _hash_outputs(output_dir)
    return results_by_group


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "phase2a1",
    )
    args = parser.parse_args()
    results = run_campaign(args.output_dir)
    print(f"Wrote Phase 2A.1 validation outputs to {args.output_dir.resolve()}")
    print(f"Completed {sum(len(group) for group in results.values())} cases")


if __name__ == "__main__":
    main()
