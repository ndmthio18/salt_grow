"""Create deterministic audit figures for Phase 2A.3b-4R."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


PALETTE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
}
FAMILIES = ("<100>", "<110>", "<111>")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig: mpl.figure.Figure, stem: Path) -> None:
    png = stem.with_suffix(".png")
    svg = stem.with_suffix(".svg")
    fig.savefig(png, dpi=301, facecolor="white", transparent=False)
    fig.savefig(svg, facecolor="white", transparent=False, metadata={"Date": None})
    plt.close(fig)
    with Image.open(png) as source:
        rgb = source.convert("RGB")
    rgb.save(png, dpi=(301, 301))
    lines = svg.read_text(encoding="utf-8").splitlines()
    svg.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def measurement_benchmark(rows: list[dict[str, str]], output: Path) -> None:
    methods = [
        ("current_area_weighted_advance_um", "Area-weighted centroid", PALETTE["blue"], "o"),
        ("fixed_region_realized_volume_advance_um", "Fixed-region realized volume", PALETTE["vermillion"], "s"),
        ("interior_height_median_advance_um", "Interior height median", PALETTE["green"], "^"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 112 / 25.4), layout="constrained", sharex="row", sharey="row")
    for row_index, dx in enumerate((10.0, 5.0)):
        for column, family in enumerate(FAMILIES):
            ax = axes[row_index, column]
            selected = sorted(
                (r for r in rows if float(r["voxel_size_um"]) == dx and r["orientation_family"] == family),
                key=lambda r: float(r["true_translation_um"]),
            )
            true = np.asarray([float(r["true_translation_um"]) for r in selected])
            ax.plot(true, true, color=PALETTE["black"], linestyle="--", linewidth=1.0, label="Analytical translation")
            for key, label, color, marker in methods:
                ax.plot(true, [float(r[key]) for r in selected], color=color, marker=marker, linewidth=1.2, markersize=3.5, label=label)
            ax.set_title(f"{family}; $\\Delta x={dx:g}$ µm")
            ax.grid(True, color="#D9D9D9", linewidth=0.5)
            if row_index == 1:
                ax.set_xlabel("Prescribed plane translation (µm)")
            if column == 0:
                ax.set_ylabel("Measured advance (µm)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=2, frameon=False)
    save(fig, output / "fig1_translated_plane_measurement")


def locked_attribution(rows: list[dict[str, str]], output: Path) -> None:
    components = [
        ("evolving_area_allocation_component_um", "Evolving-area allocation", PALETTE["sky"], "//"),
        ("binary_quantization_component_um", "Binary quantization", PALETTE["vermillion"], "xx"),
        ("measurement_definition_component_um", "Measurement definition", PALETTE["green"], ".."),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(178 / 25.4, 78 / 25.4), layout="constrained", sharey=True)
    for ax, dx in zip(axes, (10.0, 5.0), strict=True):
        selected = {r["orientation_family"]: r for r in rows if float(r["voxel_size_um"]) == dx}
        x = np.arange(len(FAMILIES))
        positive = np.zeros(len(FAMILIES))
        negative = np.zeros(len(FAMILIES))
        for key, label, color, hatch in components:
            values = np.asarray([float(selected[family][key]) for family in FAMILIES])
            bottoms = np.where(values >= 0.0, positive, negative)
            ax.bar(x, values, bottom=bottoms, color=color, edgecolor="black", linewidth=0.5, hatch=hatch, label=label)
            positive += np.where(values >= 0.0, values, 0.0)
            negative += np.where(values < 0.0, values, 0.0)
        totals = np.asarray([float(selected[family]["total_error_um"]) for family in FAMILIES])
        ax.scatter(x, totals, color=PALETTE["black"], marker="D", s=22, zorder=5, label="Total locked error")
        ax.axhline(0.0, color=PALETTE["black"], linewidth=0.8)
        ax.set_xticks(x, FAMILIES)
        ax.set_title(f"Locked $\\alpha=0$; $\\Delta x={dx:g}$ µm")
        ax.set_xlabel("Plane orientation")
        ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.5)
    axes[0].set_ylabel("Signed contribution to 20 µm advance error (µm)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=2, frameon=False)
    save(fig, output / "fig2_locked_error_attribution")


def controls_window_reservoirs(
    variants: list[dict[str, str]],
    finite: list[dict[str, str]],
    reservoirs: list[dict[str, str]],
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(178 / 25.4, 132 / 25.4), layout="constrained")
    ax = axes[0, 0]
    selected = {
        r["orientation_family"]: r
        for r in variants
        if float(r["voxel_size_um"]) == 5.0
        and float(r["growth_radius_um"]) == 70.0
        and float(r["stage_fraction"]) == 1.0
        and float(r["alpha"]) > 0.0
    }
    x = np.arange(3)
    width = 0.24
    for offset, (key, label, color, marker) in enumerate([
        ("variant_a_share", "A: analytical normal + area", PALETTE["blue"], "o"),
        ("variant_b_share", "B: reconstructed normal + ref. area", PALETTE["orange"], "s"),
        ("variant_c_share", "C: reconstructed normal + area", PALETTE["green"], "^"),
    ]):
        values = [float(selected[f][key]) for f in FAMILIES]
        ax.bar(x + (offset - 1) * width, values, width=width, color=color, edgecolor="black", linewidth=0.5, label=label)
        ax.scatter(x + (offset - 1) * width, values, color="black", marker=marker, s=10, zorder=3)
    ax.set_xticks(x, FAMILIES)
    ax.set_ylabel("Frozen allocation fraction")
    ax.set_title("a  Anisotropic probe ($\\alpha>0$), 5 µm grid")
    ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.5)
    ax.legend(frameon=False, fontsize=7)

    final = [r for r in finite if float(r["voxel_size_um"]) == 5.0]
    lookup = {r["orientation_family"]: r for r in final}
    for panel, (small_key, large_key, ylabel, title) in enumerate([
        ("small_interior_fractional_advance_um", "large_interior_fractional_advance_um", "Interior deposited-field advance (µm)", "b  Fixed interior: conserved deposited field"),
        ("small_interior_binary_advance_um", "large_interior_binary_advance_um", "Interior binary advance (µm)", "c  Fixed interior: binary geometry"),
    ]):
        ax = axes[0, 1] if panel == 0 else axes[1, 0]
        x = np.arange(3)
        small = [float(lookup[f][small_key]) for f in FAMILIES]
        large = [float(lookup[f][large_key]) for f in FAMILIES]
        ax.bar(x - 0.18, small, width=0.36, color=PALETTE["sky"], edgecolor="black", linewidth=0.5, hatch="//", label="40 µm growth radius")
        ax.bar(x + 0.18, large, width=0.36, color=PALETTE["orange"], edgecolor="black", linewidth=0.5, hatch="..", label="70 µm growth radius")
        ax.axhline(5.0, color=PALETTE["black"], linestyle="--", linewidth=1.0, label="Prescribed 5 µm")
        ax.set_xticks(x, FAMILIES)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, axis="y", color="#D9D9D9", linewidth=0.5)
        ax.legend(
            frameon=True,
            facecolor="white",
            framealpha=0.9,
            edgecolor="none",
            fontsize=7,
            loc="lower left" if panel == 0 else "upper right",
        )

    ax = axes[1, 1]
    bins = np.linspace(0.0, 1.0, 21)
    for family, front, color, linestyle in zip(
        FAMILIES, ("f100", "f110", "f111"),
        (PALETTE["blue"], PALETTE["vermillion"], PALETTE["green"]),
        ("-", "--", ":"), strict=True,
    ):
        values = [
            float(r["partial_fill_fraction"])
            for r in reservoirs
            if r["case"] == "dx5_r40" and r["front"] == front
        ]
        counts, edges = np.histogram(values, bins=bins)
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.step(centres, counts / max(1, np.sum(counts)), where="mid", color=color, linestyle=linestyle, linewidth=1.5, label=family)
    ax.set_xlabel("Residual local fill fraction")
    ax.set_ylabel("Fraction of candidate owners")
    ax.set_title("d  Unresolved fractional solid, 40 µm radius")
    ax.grid(True, color="#D9D9D9", linewidth=0.5)
    ax.legend(frameon=False)
    save(fig, output / "fig3_controls_window_fractional_solid")


def filtering_figure(rows: list[dict[str, str]], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(178 / 25.4, 76 / 25.4), layout="constrained", sharey=True)
    for ax, radius in zip(axes, (40.0, 70.0), strict=True):
        selected = [
            r for r in rows
            if float(r["voxel_size_um"]) == 5.0 and float(r["growth_radius_um"]) == radius
        ]
        for family, color, marker in zip(FAMILIES, (PALETTE["blue"], PALETTE["vermillion"], PALETTE["green"]), ("o", "s", "^"), strict=True):
            subset = sorted((r for r in selected if r["orientation_family"] == family), key=lambda r: float(r["stage_fraction"]))
            stages = np.asarray([float(r["stage_fraction"]) for r in subset])
            analytic_area = np.pi * (radius * 1e-6) ** 2
            area_ratio = np.asarray([float(r["active_area_m2"]) / analytic_area for r in subset])
            ax.plot(stages, area_ratio, color=color, marker=marker, linewidth=1.4, label=family)
        ax.axhline(1.0, color=PALETTE["black"], linestyle="--", linewidth=1.0, label="Analytical disk area")
        ax.set_xlabel("Deposited-volume stage fraction")
        ax.set_title(f"Growth radius {radius:g} µm")
        ax.grid(True, color="#D9D9D9", linewidth=0.5)
    axes[0].set_ylabel("Accepted corrected area / analytical disk area")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=4, frameon=False)
    save(fig, output / "fig4_active_interface_area")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()
    output = args.input.resolve()
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "svg.hashsalt": "phase2a3b4r",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    measurement_benchmark(read_rows(output / "translated_plane_benchmarks.csv"), output)
    locked_attribution(read_rows(output / "locked_error_attribution.csv"), output)
    controls_window_reservoirs(
        read_rows(output / "controlled_variants.csv"),
        read_rows(output / "finite_window_comparison.csv"),
        read_rows(output / "reservoir_fractions.csv"),
        output,
    )
    filtering_figure(read_rows(output / "interface_filtering.csv"), output)
    manifest = {
        "status": "general publication-style audit figures; journal requirements not specified",
        "background": "opaque white",
        "formats": ["opaque RGB PNG at at least 300 dpi metadata", "SVG with editable text"],
        "color_encoding": "Okabe-Ito colors plus markers, hatching, line styles, and panel separation",
        "uncertainty": "deterministic numerical diagnostics; no statistical interval is implied",
        "source_tables": {
            path.name: sha256(path)
            for path in sorted(output.glob("*.csv"))
        },
        "alt_text": {
            "fig1": "Six panels compare prescribed translated-plane motion with three front measurements at 10 and 5 micrometre grids; all methods show orientation-dependent voxel stair steps, while centroid and height median overlap.",
            "fig2": "Signed stacked bars decompose the locked alpha-zero 20 micrometre advance error into evolving-area allocation, binary quantization, and measurement-definition contributions at 10 and 5 micrometre grids.",
            "fig3": "Four panels compare analytical and reconstructed allocation controls, small and enlarged growth windows in a fixed interior, and residual local fill-fraction distributions.",
            "fig4": "Accepted corrected interface area relative to the analytical disk is tracked across deposition stages for three orientations and two growth radii.",
        },
    }
    (output / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    checksum_lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "SHA256SUMS.txt"
    ]
    (output / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
