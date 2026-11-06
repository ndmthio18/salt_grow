"""Run corrected pore benchmarks and selected Phase 2A.2 mask reanalysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca" / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase2 import (  # noqa: E402
    Material,
    PhysicalState,
    analyze_morphology,
    analyze_void_space,
    build_geometry_benchmarks,
    column_envelope,
    evaluate_geometry_benchmark,
    external_air_connected_mask,
    rectangular_boundary_surfaces,
    rectangular_window_mask,
)


BASE_COMMIT = "3f6fc075e88fed9a175960554b8325c13ccb2045"
PHYSICAL_SHELL_THICKNESS_M = 20.0e-6
SELECTED_RULES = ("uniform", "strong_high_coord")
SELECTED_VOXEL_SIZES_UM = (5, 10, 20)
SELECTED_RNG_SEED = 101


def _write_csv(path: Path, rows) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_hash_inventory(path: Path) -> dict[str, str]:
    inventory = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", maxsplit=1)
        inventory[relative] = digest
    return inventory


def _verify_source(path: Path, inventory: dict[str, str], root: Path) -> str:
    relative = path.relative_to(root).as_posix()
    expected = inventory.get(relative)
    if expected is None:
        raise RuntimeError(f"source hash is absent for {relative}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"source hash mismatch for {relative}")
    return actual


def _integer_ratio(numerator: float, denominator: float) -> int:
    value = numerator / denominator
    rounded = int(round(value))
    if not math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-10):
        raise ValueError("physical coordinate is not aligned with this grid")
    return rounded


def _seed_coordinates(voxel_size_m: float) -> tuple[tuple[int, int, int], ...]:
    lower_x = _integer_ratio(80.0e-6, voxel_size_m)
    lower_y = _integer_ratio(80.0e-6, voxel_size_m)
    lower_z = _integer_ratio(20.0e-6, voxel_size_m)
    side = _integer_ratio(20.0e-6, voxel_size_m)
    return tuple(
        (x, y, z)
        for x in range(lower_x, lower_x + side)
        for y in range(lower_y, lower_y + side)
        for z in range(lower_z, lower_z + side)
    )


def _state_from_mask(salt: np.ndarray, voxel_size_m: float) -> PhysicalState:
    substrate_layers = _integer_ratio(20.0e-6, voxel_size_m)
    seed_coordinates = _seed_coordinates(voxel_size_m)
    material = np.full(salt.shape, Material.EMPTY, dtype=np.uint8)
    material[:, :, :substrate_layers] = Material.SUBSTRATE
    material[salt] = Material.SALT
    seed_mask = np.zeros(salt.shape, dtype=bool)
    liquid_accessible = np.zeros(salt.shape, dtype=bool)
    for coordinate in seed_coordinates:
        if not salt[coordinate]:
            raise RuntimeError("saved mask does not contain the fixed physical seed")
        seed_mask[coordinate] = True
        liquid_accessible[coordinate] = True
    state = PhysicalState(
        material=material,
        liquid_accessible=liquid_accessible,
        seed_mask=seed_mask,
        substrate_layers=substrate_layers,
    )
    state.validate()
    return state


def _fixed_window(
    shape: tuple[int, int, int], voxel_size_m: float
):
    lower = (
        _integer_ratio(60.0e-6, voxel_size_m),
        _integer_ratio(60.0e-6, voxel_size_m),
        _integer_ratio(20.0e-6, voxel_size_m),
    )
    upper = (
        _integer_ratio(140.0e-6, voxel_size_m),
        _integer_ratio(140.0e-6, voxel_size_m),
        _integer_ratio(100.0e-6, voxel_size_m),
    )
    window = rectangular_window_mask(
        shape, lower=lower, upper_exclusive=upper
    )
    boundaries = rectangular_boundary_surfaces(
        shape,
        lower=lower,
        upper_exclusive=upper,
        name="fixed_80um_window_bottom20um_top100um",
    )
    return window, boundaries


def run_benchmarks() -> list[dict[str, object]]:
    rows = []
    for case in build_geometry_benchmarks():
        result = evaluate_geometry_benchmark(case)
        checks = (
            result["window_bottom_to_top"]
            == result["expected_window_bottom_to_top"],
            result["material_bottom_to_top"]
            == result["expected_material_bottom_to_top"],
            result["material_external_void_voxels"]
            == result["expected_material_external_void_voxels"],
            result["material_enclosed_void_voxels"]
            == result["expected_material_enclosed_void_voxels"],
            result["accessible_intergranular_void_voxels"]
            == result["expected_accessible_intergranular_void_voxels"],
        )
        rows.append({**result, "all_reference_answers_match": all(checks)})
    if not all(bool(row["all_reference_answers_match"]) for row in rows):
        raise RuntimeError("a synthetic geometry benchmark failed")
    return rows


def reanalyze_selected_masks() -> tuple[list[dict[str, object]], dict[str, str]]:
    source_root = REPOSITORY_ROOT / "phase2_ca" / "outputs" / "phase2a2"
    inventory = _read_hash_inventory(source_root / "SHA256SUMS.txt")
    source_hashes: dict[str, str] = {}
    rows = []
    for rule_name in SELECTED_RULES:
        for voxel_size_um in SELECTED_VOXEL_SIZES_UM:
            case_id = (
                f"grid_{rule_name}_{voxel_size_um}um_seed{SELECTED_RNG_SEED}"
            )
            path = source_root / "masks" / f"{case_id}_final.npy"
            relative = path.relative_to(source_root).as_posix()
            source_hashes[relative] = _verify_source(path, inventory, source_root)
            salt = np.load(path, allow_pickle=False).astype(bool, copy=False)
            voxel_size_m = voxel_size_um * 1.0e-6
            state = _state_from_mask(salt, voxel_size_m)
            window, boundaries = _fixed_window(state.shape, voxel_size_m)
            external_air = external_air_connected_mask(state)
            fixed_window = analyze_void_space(
                name="fixed_80um_window",
                analysis_mask=window,
                salt_mask=salt,
                external_air_mask=external_air,
                boundaries=boundaries,
                voxel_size_m=voxel_size_m,
            )
            analysis = analyze_morphology(
                state,
                voxel_size_m=voxel_size_m,
                boundary_surfaces=boundaries,
                fixed_window_mask=window,
                physical_shell_thickness_m=PHYSICAL_SHELL_THICKNESS_M,
            )
            column = column_envelope(state)
            column_top = int(np.max(np.argwhere(column)[:, 2]))
            implicit_top_empty = int(
                np.count_nonzero(
                    column[:, :, column_top]
                    & (state.material[:, :, column_top] == Material.EMPTY)
                )
            )
            rows.append(
                {
                    "case_id": case_id,
                    "rule_name": rule_name,
                    "voxel_size_um": voxel_size_um,
                    "rng_seed": SELECTED_RNG_SEED,
                    "fixed_window_solid_fraction": (
                        analysis.fixed_window_solid_volume_fraction
                    ),
                    "fixed_window_void_fraction": fixed_window.void_volume_fraction,
                    "fixed_window_bottom_to_top": (
                        fixed_window.bottom_to_top_void_percolates
                    ),
                    "fixed_window_exterior_void_fraction": (
                        fixed_window.external_air_connected_void_voxels
                        / fixed_window.void_voxels
                        if fixed_window.void_voxels
                        else 0.0
                    ),
                    "fixed_window_exterior_void_volume_um3": (
                        fixed_window.external_air_connected_void_voxels
                        * voxel_size_um**3
                    ),
                    "bounding_box_void_fraction": (
                        analysis.bounding_box.void_volume_fraction
                    ),
                    "column_void_fraction": (
                        analysis.column_envelope.void_volume_fraction
                    ),
                    "physical_shell_void_fraction": (
                        analysis.face_shell_envelope.void_volume_fraction
                    ),
                    "physical_shell_thickness_um": (
                        analysis.physical_shell_thickness_m * 1.0e6
                    ),
                    "physical_shell_bottom_to_top": (
                        analysis.face_shell_envelope.bottom_to_top_void_percolates
                    ),
                    "envelope_void_fraction_range": (
                        analysis.tested_envelope_void_fraction_range
                    ),
                    "bounding_box_enclosed_void_voxels": (
                        analysis.bounding_box.enclosed_void_voxels
                    ),
                    "bounding_box_enclosed_void_volume_um3": (
                        analysis.bounding_box.enclosed_void_voxels
                        * voxel_size_um**3
                    ),
                    "ambiguous_exterior_void_voxels": (
                        analysis.ambiguous_exterior_void_voxels
                    ),
                    "ambiguous_exterior_void_volume_um3": (
                        analysis.ambiguous_exterior_void_voxels
                        * voxel_size_um**3
                    ),
                    "column_implicit_top_empty_voxels": implicit_top_empty,
                    "column_pathway_interpretation": (
                        "not_applicable_envelope_top_is_salt_by_construction"
                    ),
                    "accessible_intergranular_void_voxels": "",
                    "accessible_intergranular_void_volume_um3": "",
                    "accessible_intergranular_status": (
                        "unresolved_no_persistent_grain_ids"
                    ),
                    "intrinsic_porosity": "",
                    "intrinsic_porosity_status": (
                        analysis.intrinsic_porosity_flag
                    ),
                    "liquid_supply_pathway_established": False,
                    "liquid_supply_pathway_status": (
                        "not_evaluable_without_physical_liquid_film_mask"
                    ),
                    "boundary_definition": fixed_window.boundary_definition,
                }
            )
    return rows, source_hashes


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


def run_audit(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_rows = run_benchmarks()
    mask_rows, source_hashes = reanalyze_selected_masks()
    _write_csv(output_dir / "geometry_benchmarks.csv", benchmark_rows)
    _write_csv(output_dir / "selected_mask_reanalysis.csv", mask_rows)
    analytical_property = {
        "column_envelope_definition": (
            "for each occupied x-y column, include substrate surface through "
            "that column's highest salt voxel"
        ),
        "proof": (
            "at the global maximum envelope height, every included column has "
            "a highest salt voxel on that plane; therefore the envelope-empty "
            "intersection on its implicit upper plane is empty"
        ),
        "automated_test": (
            "test_column_envelope_implicit_global_top_has_no_empty_voxels_by_construction"
        ),
        "corrected_use": (
            "void paths are evaluated against caller-defined physical surfaces; "
            "column envelopes are retained only for volume sensitivity"
        ),
    }
    (output_dir / "column_envelope_analytical_property.json").write_text(
        json.dumps(analytical_property, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "provenance": {
            "base_commit": BASE_COMMIT,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "physical_validation_claimed": False,
            "permeability_calculated": False,
            "grain_growth_implemented": False,
            "phase2b_started": False,
        },
        "benchmark_count": len(benchmark_rows),
        "selected_mask_count": len(mask_rows),
        "fixed_window": {
            "x_um": [60.0, 140.0],
            "y_um": [60.0, 140.0],
            "z_um": [20.0, 100.0],
            "inlet": "z=20 um",
            "top_outlet": "z=100 um minus one voxel",
            "physical_shell_thickness_um": 20.0,
        },
        "source_hashes": source_hashes,
        "claim_boundary": (
            "geometric void connectivity is reported separately from exterior "
            "air and does not establish liquid-film connectivity"
        ),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    _hash_outputs(output_dir)
    return {
        "benchmarks": benchmark_rows,
        "selected_masks": mask_rows,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "phase2a2r",
    )
    args = parser.parse_args()
    results = run_audit(args.output_dir)
    print(
        f"Wrote {len(results['benchmarks'])} benchmark results and "
        f"{len(results['selected_masks'])} selected-mask analyses to "
        f"{args.output_dir.resolve()}"
    )


if __name__ == "__main__":
    main()
