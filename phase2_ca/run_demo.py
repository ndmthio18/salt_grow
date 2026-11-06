"""Run a deterministic small Phase 2A evaporation-driven growth demonstration."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "phase2_ca" / "src"))
sys.path.insert(0, str(REPOSITORY_ROOT / "phase1_mass_budget" / "src"))

from salt_phase1 import EvaporationParameters, compute_mass_budget  # noqa: E402
from salt_phase2 import SaltGrowthCA, initialize_physical_state  # noqa: E402


SHAPE = (16, 16, 12)
SUBSTRATE_LAYERS = 1
VOXEL_SIZE_M = 10.0e-6
SOLID_DENSITY_KG_M3 = 2160.0
TIMESTEP_S = 60.0
STEPS = 12
SEEDS = ((4, 4, 1), (4, 11, 1), (11, 4, 1), (11, 11, 1))


def generate_demo(output_path: Path) -> dict[str, object]:
    state = initialize_physical_state(
        shape=SHAPE,
        substrate_layers=SUBSTRATE_LAYERS,
        seed_coordinates=SEEDS,
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=SOLID_DENSITY_KG_M3,
        voxel_size_m=VOXEL_SIZE_M,
        enforce_wetting=True,
        tie_break_seed=20261106,
    )
    mass_budget = compute_mass_budget(
        EvaporationParameters(
            temperature_k=298.15,
            water_activity=0.753,
            relative_humidity=0.30,
            mass_transfer_coefficient_m_s=0.005,
            evaporation_area_m2=(SHAPE[0] * VOXEL_SIZE_M)
            * (SHAPE[1] * VOXEL_SIZE_M),
            timestep_s=TIMESTEP_S,
            salt_mass_fraction=0.264,
            solid_density_kg_m3=SOLID_DENSITY_KG_M3,
        )
    )

    steps: list[dict[str, object]] = []
    for step_index in range(1, STEPS + 1):
        result = ca.grow(mass_budget.precipitated_salt_mass_kg)
        steps.append(
            {
                "step": step_index,
                "time_s": step_index * TIMESTEP_S,
                "candidate_voxel_count": result.candidate_voxel_count,
                "selected_coordinates": result.selected_coordinates,
                "transaction": asdict(result.transaction),
                "morphology": asdict(result.morphology),
            }
        )

    ca.assert_invariants()
    final_salt_coordinates = [
        tuple(int(value) for value in coordinate)
        for coordinate in np.argwhere(ca.state.material == 2)
    ]
    summary: dict[str, object] = {
        "provenance": {
            "reviewed_phase1_1_commit": "ae4c7ae6bb641e2f2bbda5a2d01e434b473966d3",
            "ca_implementation": "independent CPU Phase 2A prototype",
            "physical_validation_claimed": False,
        },
        "configuration": {
            "shape": SHAPE,
            "substrate_layers": SUBSTRATE_LAYERS,
            "voxel_size_m": VOXEL_SIZE_M,
            "solid_density_kg_m3": SOLID_DENSITY_KG_M3,
            "timestep_s": TIMESTEP_S,
            "steps": STEPS,
            "seed_coordinates": SEEDS,
            "enforce_wetting": True,
            "neighbor_rule": "six-face",
            "candidate_priority": (
                "descending salt coordination then fixed-seed mixed-coordinate priority"
            ),
            "tie_break_seed": 20261106,
        },
        "evaporation_mass_budget": asdict(mass_budget),
        "steps": steps,
        "final_salt_coordinates": final_salt_coordinates,
        "final_morphology": asdict(ca.morphology_diagnostics()),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    with (output_path.parent / "SHA256SUMS.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"{digest}  {output_path.name}\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "demo_summary.json",
    )
    args = parser.parse_args()
    summary = generate_demo(args.output)
    final = summary["final_morphology"]
    print(f"Wrote {args.output.resolve()}")
    print(
        "salt voxels="
        f"{final['total_salt_voxels']} "
        f"new={final['new_salt_voxels']} "
        f"height={final['maximum_height_voxels_above_substrate']} voxels "
        f"closure={final['mass_closure_error_kg']:.3e} kg"
    )


if __name__ == "__main__":
    main()
