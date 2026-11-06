"""Post-lock robustness audit for the Phase 2A.3a volume-resolution gate.

This diagnostic does not alter the locked protocol or its pass/fail result.  It
tests whether a favourable single alignment is stable to independent cube
rotation and sub-voxel translation.
"""

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


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT.parent / "phase1_mass_budget" / "src"))

from salt_phase2.orientation_verification import geometry_catalogue  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    catalogue = geometry_catalogue()
    width_m = 2.0 * catalogue["cube_rot_a"].half_width_m
    cells_across_values = (1, 2, 4, 8, 16, 32)
    offsets = ((0.0, 0.0, 0.0), (0.37, -0.23, 0.41), (-0.43, 0.31, -0.19))
    case_names = ("cube_axis", "cube_rot_a", "cube_rot_b")
    raw_rows = []
    summary_rows = []
    for cells_across in cells_across_values:
        voxel_m = width_m / cells_across
        errors = []
        for case_name in case_names:
            base = catalogue[case_name]
            for offset_index, fractions in enumerate(offsets):
                geometry = replace(
                    base,
                    domain_length_m=3.0 * width_m,
                    centre_m=tuple(fraction * voxel_m for fraction in fractions),
                )
                mask = geometry.voxelize(voxel_m)
                volume_m3 = np.count_nonzero(mask) * voxel_m**3
                error = abs(volume_m3 / geometry.exact_volume_m3 - 1.0)
                errors.append(error)
                raw_rows.append(
                    {
                        "cells_across_seed": cells_across,
                        "voxel_size_um": voxel_m * 1.0e6,
                        "geometry": case_name,
                        "subvoxel_offset_index": offset_index,
                        "offset_x_voxels": fractions[0],
                        "offset_y_voxels": fractions[1],
                        "offset_z_voxels": fractions[2],
                        "volume_relative_error": error,
                    }
                )
        summary_rows.append(
            {
                "cells_across_seed": cells_across,
                "voxel_size_um": voxel_m * 1.0e6,
                "median_volume_relative_error": float(np.median(errors)),
                "maximum_volume_relative_error": float(np.max(errors)),
                "alignment_count": len(errors),
            }
        )
    for filename, rows in (
        ("volume_alignment_cases.csv", raw_rows),
        ("volume_alignment_summary.csv", summary_rows),
    ):
        with (output / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    robust_volume_minimum = min(
        row["cells_across_seed"]
        for row in summary_rows
        if row["maximum_volume_relative_error"] <= 0.10
        and row["cells_across_seed"] > 2
    )
    requirements = {
        "status": "post-lock robustness diagnostic; locked thresholds and results unchanged",
        "minimum_cells_for_robust_volume_representation": robust_volume_minimum,
        "minimum_cells_for_surface_normal_estimation": 8,
        "minimum_cells_for_facet_orientation": 16,
        "minimum_cells_across_narrow_gap": 5,
        "one_or_two_voxel_seed": "NO-GO for orientation: the 120 um physical fitting support cannot be formed, and maximum volume error over independent alignment/rotation is 100% and 25%, respectively.",
        "scope": "Numerical geometry requirements only; not a transport, kinetics, or physical-calibration claim.",
    }
    (output / "robust_resolution_requirements.json").write_text(
        json.dumps(requirements, indent=2), encoding="utf-8"
    )
    mpl.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    fig, axis = plt.subplots(figsize=(4.7, 3.3), layout="constrained")
    x = np.array([row["cells_across_seed"] for row in summary_rows])
    median = 100 * np.array([row["median_volume_relative_error"] for row in summary_rows])
    maximum = 100 * np.array([row["maximum_volume_relative_error"] for row in summary_rows])
    axis.plot(x, maximum, "o-", color="#D55E00", label="maximum over 9 alignments")
    axis.plot(x, median, "s--", color="#0072B2", label="median over 9 alignments")
    axis.axhline(10.0, color="#333333", linestyle=":", label="predeclared 10% volume gate")
    axis.axvspan(0.8, 2.2, color="#999999", alpha=0.18, label="excluded from facet evidence")
    axis.set_xscale("log", base=2)
    axis.set(
        xlabel="Seed width (voxels)",
        ylabel="Absolute volume error (%)",
        xticks=x,
        xticklabels=[str(value) for value in x],
    )
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, fontsize=7)
    for suffix in ("png", "svg"):
        fig.savefig(output / f"volume_alignment_robustness.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)
    locked_protocol = ROOT / "configs" / "phase2a3a_locked_protocol.json"
    manifest = {
        "post_lock_diagnostic": True,
        "locked_protocol_sha256": sha256(locked_protocol),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": mpl.__version__,
        "geometry_count": len(case_names),
        "offset_count_per_geometry": len(offsets),
        "resolution_count": len(cells_across_values),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    checksum_lines = []
    for path in sorted(output.iterdir()):
        if path.name != "SHA256SUMS.txt":
            checksum_lines.append(f"{sha256(path)}  {path.name}")
    (output / "SHA256SUMS.txt").write_text(
        "\n".join(checksum_lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
