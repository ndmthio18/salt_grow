"""Truthful static visualizations for Phase 2A.1 morphology comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from PIL import Image

from .validation import CaseResult


SALT_COLOR = "#0072B2"
SEED_COLOR = "#E69F00"
GRID_COLOR = "#666666"


def _save_opaque_png(figure: plt.Figure, output_path: Path) -> None:
    """Save a white-background PNG and remove the unnecessary alpha channel."""

    figure.savefig(
        output_path,
        facecolor="white",
        transparent=False,
        metadata={"Software": "salt_phase2a1 morphology validation"},
    )
    with Image.open(output_path) as image:
        dpi = image.info.get("dpi", (200.0, 200.0))
        image.convert("RGB").save(output_path, dpi=dpi)


def _slice_state(mask: np.ndarray, seeds: np.ndarray, axis: int, index: int) -> np.ndarray:
    state = np.zeros(mask.shape, dtype=np.uint8)
    state[mask] = 1
    state[seeds] = 2
    if axis == 2:
        return state[:, :, index].T
    if axis == 1:
        return state[:, index, :].T
    raise ValueError("axis must be 1 or 2")


def render_case_rows(
    cases: Sequence[CaseResult],
    labels: Sequence[str],
    output_path: Path,
    *,
    title: str,
) -> list[dict[str, object]]:
    """Render one 3D view and orthogonal cross-sections for every case row."""

    if len(cases) != len(labels) or not cases:
        raise ValueError("cases and labels must have the same non-zero length")
    mpl.rcParams.update(
        {
            "font.size": 8,
            "figure.dpi": 120,
            "savefig.dpi": 200,
            "svg.hashsalt": "salt-phase2a1",
        }
    )
    figure = plt.figure(
        figsize=(11.5, 3.1 * len(cases)),
        layout="constrained",
        facecolor="white",
    )
    subfigures = figure.subfigures(len(cases), 1, squeeze=False)
    slice_records: list[dict[str, object]] = []
    cmap = ListedColormap(["#FFFFFF", SALT_COLOR, SEED_COLOR])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    for row, (case, label) in enumerate(zip(cases, labels, strict=True)):
        subfigure = subfigures[row, 0]
        axes = subfigure.subplots(1, 3, subplot_kw={"projection": None})
        axes[0].remove()
        axis_3d = subfigure.add_subplot(1, 3, 1, projection="3d")
        horizontal_axis = axes[1]
        vertical_axis = axes[2]
        salt = case.final_salt_mask
        seeds = case.seed_mask
        colors = np.zeros(salt.shape + (4,), dtype=float)
        colors[salt] = mpl.colors.to_rgba(SALT_COLOR, alpha=0.72)
        colors[seeds] = mpl.colors.to_rgba(SEED_COLOR, alpha=0.95)
        axis_3d.voxels(
            salt,
            facecolors=colors,
            edgecolor=(1.0, 1.0, 1.0, 0.22),
            linewidth=0.15,
        )
        axis_3d.set(
            xlim=(0, salt.shape[0]),
            ylim=(0, salt.shape[1]),
            zlim=(0, salt.shape[2]),
            xlabel="x",
            ylabel="y",
            zlabel="z",
        )
        axis_3d.set_box_aspect(salt.shape)
        axis_3d.view_init(elev=24, azim=-58)
        axis_3d.set_title(f"{label}: 3D salt")

        salt_by_z = np.count_nonzero(salt, axis=(0, 1))
        horizontal_z = int(np.argmax(salt_by_z))
        salt_by_y = np.count_nonzero(salt, axis=(0, 2))
        vertical_y = int(np.argmax(salt_by_y))
        horizontal_axis.imshow(
            _slice_state(salt, seeds, 2, horizontal_z),
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
            aspect="equal",
        )
        horizontal_axis.set(
            title=f"horizontal z={horizontal_z}",
            xlabel="x voxel",
            ylabel="y voxel",
            xlim=(-0.5, salt.shape[0] - 0.5),
            ylim=(-0.5, salt.shape[1] - 0.5),
        )
        vertical_axis.imshow(
            _slice_state(salt, seeds, 1, vertical_y),
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
            aspect="equal",
        )
        vertical_axis.set(
            title=f"vertical y={vertical_y}",
            xlabel="x voxel",
            ylabel="z voxel",
            xlim=(-0.5, salt.shape[0] - 0.5),
            ylim=(-0.5, salt.shape[2] - 0.5),
        )
        for axis in (horizontal_axis, vertical_axis):
            axis.grid(False)
        subfigure.legend(
            handles=[
                Patch(facecolor=SALT_COLOR, label="new/existing salt"),
                Patch(facecolor=SEED_COLOR, label="initial seed"),
            ],
            loc="upper right",
            frameon=False,
        )
        slice_records.append(
            {
                "case_id": case.configuration.case_id,
                "horizontal_z_index": horizontal_z,
                "vertical_y_index": vertical_y,
            }
        )
    figure.suptitle(title, fontsize=11)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _save_opaque_png(figure, output_path)
    plt.close(figure)
    return slice_records


def plot_rule_metrics(
    cases: Sequence[CaseResult],
    output_path: Path,
    *,
    rule_order: Sequence[str],
) -> None:
    """Show raw simulation replicates and descriptive mean plus sample SD."""

    metrics = (
        ("maximum_height_voxels", "Maximum height [voxels]"),
        ("footprint_voxels", "Footprint [voxels]"),
        ("bbox_salt_fraction", "Bounding-box salt fraction [-]"),
        ("column_salt_fraction", "Column-envelope salt fraction [-]"),
    )
    palette = ["#000000", "#0072B2", "#D55E00", "#009E73"]
    markers = ["X", "o", "s", "^"]
    figure, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), layout="constrained")
    for axis, (metric, ylabel) in zip(axes.flat, metrics, strict=True):
        for group_index, rule in enumerate(rule_order):
            values = np.array(
                [
                    float(case.summary[metric])
                    for case in cases
                    if case.configuration.rule_name == rule
                ]
            )
            offsets = (
                np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else np.array([0.0])
            )
            axis.scatter(
                group_index + offsets,
                values,
                color=palette[group_index],
                marker=markers[group_index],
                s=28,
                alpha=0.82,
                label=rule if metric == metrics[0][0] else None,
                zorder=3,
            )
            if len(values) > 1:
                axis.errorbar(
                    group_index,
                    float(np.mean(values)),
                    yerr=float(np.std(values, ddof=1)),
                    color=palette[group_index],
                    marker="D",
                    markersize=4,
                    capsize=3,
                    linewidth=1.2,
                    zorder=4,
                )
        axis.set(
            ylabel=ylabel,
            xticks=range(len(rule_order)),
            xticklabels=rule_order,
        )
        axis.set_ylim(bottom=0.0)
        axis.tick_params(axis="x", rotation=20)
        axis.grid(axis="y", color=GRID_COLOR, alpha=0.18)
    axes.flat[0].legend(frameon=False, fontsize=7)
    figure.suptitle("Morphology across rules (points: seeds; diamonds/error bars: mean ± SD)")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _save_opaque_png(figure, output_path)
    plt.close(figure)
