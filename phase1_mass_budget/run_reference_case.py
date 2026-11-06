"""Generate deterministic Phase 1 reference tables and figures."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt


SOURCE_ROOT = Path(__file__).resolve().parent / "src"
sys.path.insert(0, str(SOURCE_ROOT))

from salt_phase1 import (  # noqa: E402
    EvaporationParameters,
    MassToVoxelAllocator,
    compute_mass_budget,
    diagnose_domain_capacity,
    required_domain_height_m,
)


TEMPERATURE_K = 298.15
WATER_ACTIVITY = 0.753
SALT_MASS_FRACTION = 0.264
SOLID_DENSITY_KG_M3 = 2160.0
MASS_TRANSFER_COEFFICIENT_M_S = 0.005
GRID_SIZE = 64
VOXEL_SIZE_M = 10.0e-6
EVAPORATION_AREA_M2 = (GRID_SIZE * VOXEL_SIZE_M) ** 2
TIMESTEP_S = 60.0
RELATIVE_HUMIDITIES = (0.30, 0.50, 0.70, 0.753, 0.80)
TOTAL_DURATION_S = 24.0 * 3600.0


def parameters(relative_humidity: float, timestep_s: float = TIMESTEP_S):
    return EvaporationParameters(
        temperature_k=TEMPERATURE_K,
        water_activity=WATER_ACTIVITY,
        relative_humidity=relative_humidity,
        mass_transfer_coefficient_m_s=MASS_TRANSFER_COEFFICIENT_M_S,
        evaporation_area_m2=EVAPORATION_AREA_M2,
        timestep_s=timestep_s,
        salt_mass_fraction=SALT_MASS_FRACTION,
        solid_density_kg_m3=SOLID_DENSITY_KG_M3,
    )


def build_time_grid(total_duration_s: float, timestep_s: float) -> list[float]:
    """Return output times in seconds, including the configured final time."""

    if not math.isfinite(total_duration_s) or total_duration_s < 0.0:
        raise ValueError("total_duration_s must be finite and non-negative")
    if not math.isfinite(timestep_s) or timestep_s <= 0.0:
        raise ValueError("timestep_s must be finite and positive")
    full_steps = int(math.floor(total_duration_s / timestep_s))
    times = [index * timestep_s for index in range(full_steps + 1)]
    if not math.isclose(times[-1], total_duration_s, rel_tol=0.0, abs_tol=1.0e-12):
        times.append(total_duration_s)
    else:
        times[-1] = total_duration_s
    return times


def generate(
    output_dir: Path,
    *,
    total_duration_s: float = TOTAL_DURATION_S,
    timestep_s: float = TIMESTEP_S,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []

    for rh in RELATIVE_HUMIDITIES:
        result = compute_mass_budget(parameters(rh))
        allocator = MassToVoxelAllocator(
            solid_density_kg_m3=SOLID_DENSITY_KG_M3,
            voxel_size_m=VOXEL_SIZE_M,
            initial_seed_voxel_count=0,
        )
        allocation = allocator.allocate(result.precipitated_salt_mass_kg)
        rows.append(
            {
                "relative_humidity": rh,
                "regime": result.regime,
                "c_sat_kg_m3": result.saturation_vapor_concentration_kg_m3,
                "evaporation_flux_kg_m2_s": result.evaporation_flux_kg_m2_s,
                "water_rate_kg_s": result.water_evaporation_rate_kg_s,
                "water_mass_per_60s_kg": result.evaporated_water_mass_kg,
                "salt_rate_kg_s": result.salt_precipitation_rate_kg_s,
                "salt_mass_per_60s_kg": result.precipitated_salt_mass_kg,
                "salt_volume_per_60s_m3": result.precipitated_salt_volume_m3,
                "new_voxels_per_60s": allocation.new_solid_voxels,
                "residual_mass_after_60s_kg": allocation.residual_mass_kg,
            }
        )

    sweep_path = output_dir / "rh_sweep.csv"
    with sweep_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    time_seconds = build_time_grid(total_duration_s, timestep_s)
    time_hours = [value / 3600.0 for value in time_seconds]
    cumulative_rows: list[dict[str, float]] = []
    rates = {
        rh: compute_mass_budget(parameters(rh)).salt_precipitation_rate_kg_s
        for rh in RELATIVE_HUMIDITIES
    }
    for time_h in time_hours:
        row: dict[str, float] = {"time_h": time_h}
        for rh in RELATIVE_HUMIDITIES:
            row[f"rh_{rh:g}_salt_mass_kg"] = rates[rh] * time_h * 3600.0
        cumulative_rows.append(row)

    cumulative_path = output_dir / "cumulative_salt_mass.csv"
    with cumulative_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cumulative_rows[0]))
        writer.writeheader()
        writer.writerows(cumulative_rows)

    rhs = [float(row["relative_humidity"]) for row in rows]
    fluxes = [float(row["evaporation_flux_kg_m2_s"]) for row in rows]
    salt_rates = [float(row["salt_rate_kg_s"]) for row in rows]

    plt.rcParams.update(
        {"font.size": 9, "figure.dpi": 150, "svg.hashsalt": "salt-phase1-v1"}
    )
    figure, axes = plt.subplots(1, 3, figsize=(10.5, 3.3))
    axes[0].plot(rhs, fluxes, "o-", color="#2166ac")
    axes[0].axvline(WATER_ACTIVITY, color="0.45", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Environmental RH [-]", ylabel="Water flux [kg m$^{-2}$ s$^{-1}$]")
    axes[0].grid(alpha=0.25)

    axes[1].plot(rhs, salt_rates, "o-", color="#b2182b")
    axes[1].axvline(WATER_ACTIVITY, color="0.45", linestyle="--", linewidth=1)
    axes[1].set(xlabel="Environmental RH [-]", ylabel="NaCl rate [kg s$^{-1}$]")
    axes[1].grid(alpha=0.25)

    for rh in RELATIVE_HUMIDITIES:
        axes[2].plot(
            time_hours,
            [rates[rh] * value * 3600.0 for value in time_hours],
            label=f"RH={rh:g}",
        )
    axes[2].set(xlabel="Time [h]", ylabel="Cumulative NaCl mass [kg]")
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=7, frameon=False)
    figure.tight_layout()
    figure.savefig(
        output_dir / "phase1_reference_results.png",
        bbox_inches="tight",
        metadata={"Software": "salt_phase1 deterministic reference"},
    )
    svg_path = output_dir / "phase1_reference_results.svg"
    figure.savefig(
        svg_path,
        bbox_inches="tight",
        metadata={"Date": None, "Creator": "salt_phase1 deterministic reference"},
    )
    plt.close(figure)
    # Matplotlib emits path-data lines with trailing spaces. Normalize the text
    # artifact so repository whitespace checks stay meaningful and reproducible.
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )

    reference_rate = rates[0.30] / SOLID_DENSITY_KG_M3
    cumulative_reference_volume = reference_rate * total_duration_s
    capacity = diagnose_domain_capacity(
        grid_shape=(GRID_SIZE, GRID_SIZE, GRID_SIZE),
        voxel_size_m=VOXEL_SIZE_M,
        cumulative_precipitated_salt_volume_m3=cumulative_reference_volume,
        initial_seed_voxel_count=0,
        maximum_permitted_solid_volume_fraction=1.0,
        precipitation_volume_rate_m3_s=reference_rate,
    )
    required_heights = {
        f"solid_fraction_{fraction:g}": required_domain_height_m(
            total_solid_volume_m3=capacity.total_solid_volume_m3,
            base_grid_shape=(GRID_SIZE, GRID_SIZE),
            voxel_size_m=VOXEL_SIZE_M,
            target_solid_volume_fraction=fraction,
        )
        for fraction in (0.3, 0.5, 0.7)
    }

    summary = {
        "provenance": {
            "reference_ca_repository": "https://github.com/avredkov/3D-CA-Crystal-Life",
            "reference_ca_commit": "995f197192988483140e6f7e773c496aa16cf76b",
            "reference_ca_source_copied": False,
            "gpu_baseline_executed": False,
            "vapor_pressure_relation": "Buck (1981) over liquid water",
        },
        "assumptions": {
            "temperature_k": TEMPERATURE_K,
            "water_activity": WATER_ACTIVITY,
            "salt_mass_fraction": SALT_MASS_FRACTION,
            "solid_density_kg_m3": SOLID_DENSITY_KG_M3,
            "mass_transfer_coefficient_m_s": MASS_TRANSFER_COEFFICIENT_M_S,
            "grid": [GRID_SIZE, GRID_SIZE, GRID_SIZE],
            "voxel_size_m": VOXEL_SIZE_M,
            "evaporation_area_m2": EVAPORATION_AREA_M2,
            "timestep_s": timestep_s,
            "total_duration_s": total_duration_s,
        },
        "salt_to_water_mass_ratio": SALT_MASS_FRACTION / (1.0 - SALT_MASS_FRACTION),
        "voxel_mass_kg": SOLID_DENSITY_KG_M3 * VOXEL_SIZE_M**3,
        "time_grid_point_count": len(time_seconds),
        "capacity_diagnostic_rh_0.3": asdict(capacity),
        "required_domain_heights_m": required_heights,
        "rows": rows,
    }
    with (output_dir / "reference_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    artifact_names = (
        "rh_sweep.csv",
        "cumulative_salt_mass.csv",
        "reference_summary.json",
        "phase1_reference_results.png",
        "phase1_reference_results.svg",
    )
    with (output_dir / "SHA256SUMS.txt").open("w", encoding="utf-8") as handle:
        for name in artifact_names:
            digest = hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
            handle.write(f"{digest}  {name}\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "outputs" / "reference",
    )
    parser.add_argument("--total-duration-s", type=float, default=TOTAL_DURATION_S)
    parser.add_argument("--timestep-s", type=float, default=TIMESTEP_S)
    args = parser.parse_args()
    summary = generate(
        args.output_dir,
        total_duration_s=args.total_duration_s,
        timestep_s=args.timestep_s,
    )
    print(f"Wrote Phase 1 reference outputs to {args.output_dir.resolve()}")
    print(f"c_sat(298.15 K) = {summary['rows'][0]['c_sat_kg_m3']:.8e} kg/m3")
    for row in summary["rows"]:
        print(
            f"RH={row['relative_humidity']:.3f} "
            f"Jw={row['evaporation_flux_kg_m2_s']:.8e} kg/(m2 s) "
            f"Ms60={row['salt_mass_per_60s_kg']:.8e} kg "
            f"Nnew={row['new_voxels_per_60s']} "
            f"regime={row['regime']}"
        )


if __name__ == "__main__":
    main()
