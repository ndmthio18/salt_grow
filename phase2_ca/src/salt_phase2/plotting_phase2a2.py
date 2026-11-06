"""Phase 2A.2 figures for grid sensitivity and pore-space classification."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
import numpy as np
from numpy.typing import NDArray
from PIL import Image


EXTERIOR_COLOR = "#DDEFF2"
SUBSTRATE_COLOR = "#777777"
SALT_COLOR = "#0072B2"
PORE_COLOR = "#D55E00"
SEED_COLOR = "#E69F00"
GRID_COLOR = "#666666"


@dataclass(frozen=True)
class ClassifiedMorphologyView:
    label: str
    salt_mask: NDArray[np.bool_]
    seed_mask: NDArray[np.bool_]
    candidate_pore_mask: NDArray[np.bool_]
    substrate_layers: int
    voxel_size_um: float


def _save_rgb(figure: plt.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        output_path,
        facecolor="white",
        transparent=False,
        dpi=220,
        metadata={"Software": "salt_phase2a2 morphology feasibility audit"},
    )
    with Image.open(output_path) as image:
        dpi = image.info.get("dpi", (220.0, 220.0))
        image.convert("RGB").save(output_path, dpi=dpi)


def _classified_state(view: ClassifiedMorphologyView) -> NDArray[np.uint8]:
    state = np.zeros(view.salt_mask.shape, dtype=np.uint8)
    state[:, :, : view.substrate_layers] = 1
    state[view.candidate_pore_mask] = 3
    state[view.salt_mask] = 2
    state[view.seed_mask] = 4
    return state


def render_classified_morphologies(
    views: Sequence[ClassifiedMorphologyView],
    output_path: Path,
    *,
    title: str,
) -> list[dict[str, object]]:
    """Render true voxel geometry plus orthogonal categorical cross-sections."""

    if not views:
        raise ValueError("at least one morphology view is required")
    cmap = ListedColormap(
        [EXTERIOR_COLOR, SUBSTRATE_COLOR, SALT_COLOR, PORE_COLOR, SEED_COLOR]
    )
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1.0), cmap.N)
    figure = plt.figure(
        figsize=(11.5, 3.0 * len(views)),
        layout="constrained",
        facecolor="white",
    )
    subfigures = figure.subfigures(len(views), 1, squeeze=False)
    records: list[dict[str, object]] = []
    for row, view in enumerate(views):
        subfigure = subfigures[row, 0]
        axes = subfigure.subplots(1, 3)
        axes[0].remove()
        axis_3d = subfigure.add_subplot(1, 3, 1, projection="3d")
        horizontal_axis, vertical_axis = axes[1], axes[2]
        salt = view.salt_mask
        seed = view.seed_mask
        pore = view.candidate_pore_mask
        colors = np.zeros(salt.shape + (4,), dtype=float)
        colors[salt] = mpl.colors.to_rgba(SALT_COLOR, 0.78)
        colors[seed] = mpl.colors.to_rgba(SEED_COLOR, 0.98)
        colors[pore] = mpl.colors.to_rgba(PORE_COLOR, 0.30)
        visible = salt | pore
        axis_3d.voxels(
            visible,
            facecolors=colors,
            edgecolor=(1.0, 1.0, 1.0, 0.16),
            linewidth=0.10,
        )
        nx, ny, nz = salt.shape
        x_grid, y_grid = np.meshgrid(
            np.arange(nx + 1), np.arange(ny + 1), indexing="ij"
        )
        axis_3d.plot_surface(
            x_grid,
            y_grid,
            np.full_like(x_grid, view.substrate_layers, dtype=float),
            color=SUBSTRATE_COLOR,
            alpha=0.22,
            linewidth=0,
            shade=False,
        )
        axis_3d.set(
            xlim=(0, nx),
            ylim=(0, ny),
            zlim=(0, nz),
            xlabel="x [voxel]",
            ylabel="y [voxel]",
            zlabel="z [voxel]",
            title=f"{view.label}: 3D",
        )
        axis_3d.set_box_aspect(salt.shape)
        axis_3d.view_init(elev=24, azim=-58)

        categorical = _classified_state(view)
        above_substrate = np.arange(nz) >= view.substrate_layers
        salt_pore_by_z = np.count_nonzero(
            (salt | pore)[:, :, above_substrate], axis=(0, 1)
        )
        horizontal_z = view.substrate_layers + int(np.argmax(salt_pore_by_z))
        salt_pore_by_y = np.count_nonzero(salt | pore, axis=(0, 2))
        vertical_y = int(np.argmax(salt_pore_by_y))
        extent_xy = [
            0.0,
            nx * view.voxel_size_um,
            0.0,
            ny * view.voxel_size_um,
        ]
        horizontal_axis.imshow(
            categorical[:, :, horizontal_z].T,
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
            extent=extent_xy,
            aspect="equal",
        )
        horizontal_axis.set(
            title=f"horizontal z={horizontal_z * view.voxel_size_um:.0f} µm",
            xlabel="x [µm]",
            ylabel="y [µm]",
        )
        extent_xz = [
            0.0,
            nx * view.voxel_size_um,
            0.0,
            nz * view.voxel_size_um,
        ]
        vertical_axis.imshow(
            categorical[:, vertical_y, :].T,
            origin="lower",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
            extent=extent_xz,
            aspect="equal",
        )
        vertical_axis.set(
            title=f"vertical y={vertical_y * view.voxel_size_um:.0f} µm",
            xlabel="x [µm]",
            ylabel="z [µm]",
        )
        for axis in (horizontal_axis, vertical_axis):
            axis.grid(False)
        records.append(
            {
                "label": view.label,
                "horizontal_z_index": horizontal_z,
                "vertical_y_index": vertical_y,
                "voxel_size_um": view.voxel_size_um,
            }
        )
    figure.suptitle(title, fontsize=11)
    figure.legend(
        handles=[
            Patch(facecolor=SALT_COLOR, label="salt"),
            Patch(facecolor=SEED_COLOR, label="initial seed"),
            Patch(facecolor=SUBSTRATE_COLOR, label="substrate"),
            Patch(facecolor=PORE_COLOR, label="candidate pore (column mask)"),
            Patch(facecolor=EXTERIOR_COLOR, label="exterior air (transparent in 3D)"),
        ],
        loc="outside lower center",
        ncols=5,
        frameon=False,
        fontsize=7,
    )
    _save_rgb(figure, output_path)
    plt.close(figure)
    return records


def plot_grid_sensitivity(
    rows: Sequence[Mapping[str, object]],
    output_path: Path,
) -> None:
    """Show every simulation seed plus descriptive mean and sample SD."""

    metrics = (
        ("maximum_height_um", "Maximum height [µm]"),
        ("footprint_area_um2", "Footprint area [µm²]"),
        ("surface_to_volume_ratio_per_mm", "Surface/volume [mm⁻¹]"),
        ("characteristic_solid_size_um", "Characteristic size 6V/A [µm]"),
        ("bbox_void_fraction", "Bounding-box void fraction [-]"),
        ("column_void_fraction", "Column-envelope void fraction [-]"),
    )
    rules = ("uniform", "strong_high_coord")
    styles = {
        "uniform": ("#0072B2", "o", "-"),
        "strong_high_coord": ("#D55E00", "s", "--"),
    }
    voxel_sizes = (5.0, 10.0, 20.0)
    figure, axes = plt.subplots(2, 3, figsize=(10.8, 6.6), layout="constrained")
    for axis, (metric, ylabel) in zip(axes.flat, metrics, strict=True):
        for rule in rules:
            color, marker, line_style = styles[rule]
            means = []
            standard_deviations = []
            for voxel_size in voxel_sizes:
                values = np.array(
                    [
                        float(row[metric])
                        for row in rows
                        if row["rule_name"] == rule
                        and float(row["voxel_size_um"]) == voxel_size
                    ]
                )
                offsets = np.linspace(-0.35, 0.35, len(values))
                axis.scatter(
                    voxel_size + offsets,
                    values,
                    color=color,
                    marker=marker,
                    s=24,
                    alpha=0.65,
                    zorder=3,
                )
                means.append(float(np.mean(values)))
                standard_deviations.append(
                    float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
                )
            axis.errorbar(
                voxel_sizes,
                means,
                yerr=standard_deviations,
                color=color,
                marker=marker,
                linestyle=line_style,
                capsize=3,
                linewidth=1.3,
                label=rule,
                zorder=4,
            )
        axis.set(
            xlabel="Voxel size [µm]",
            ylabel=ylabel,
            xticks=voxel_sizes,
        )
        axis.grid(axis="y", color=GRID_COLOR, alpha=0.18)
        axis.set_ylim(bottom=0.0)
    axes[0, 0].legend(
        frameon=False,
        fontsize=7,
        title="points: seeds; line: mean ± sample SD",
    )
    figure.suptitle("Grid sensitivity at fixed physical domain, seed, duration, and supply")
    _save_rgb(figure, output_path)
    plt.close(figure)


def plot_envelope_sensitivity(
    rows: Sequence[Mapping[str, object]],
    output_path: Path,
) -> None:
    """Compare void fractions and pathway flags across explicit envelopes."""

    rules = ("uniform", "strong_high_coord")
    envelopes = ("bounding_box", "column", "face_shell")
    labels = ("Bounding box", "Column", "Face shell")
    colors = {"uniform": "#0072B2", "strong_high_coord": "#D55E00"}
    markers = {"uniform": "o", "strong_high_coord": "s"}
    figure, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), layout="constrained")
    for rule_index, rule in enumerate(rules):
        selected = [row for row in rows if row["rule_name"] == rule]
        for envelope_index, envelope in enumerate(envelopes):
            values = np.array(
                [float(row[f"{envelope}_void_fraction"]) for row in selected]
            )
            offsets = np.linspace(-0.08, 0.08, len(values))
            x = envelope_index + (rule_index - 0.5) * 0.22
            axes[0].scatter(
                x + offsets,
                values,
                color=colors[rule],
                marker=markers[rule],
                alpha=0.7,
                label=rule if envelope_index == 0 else None,
            )
            axes[0].errorbar(
                x,
                float(np.mean(values)),
                yerr=float(np.std(values, ddof=1)),
                color=colors[rule],
                marker="D",
                capsize=3,
            )
        pathway = np.array(
            [
                int(bool(row["column_continuous_candidate_path"]))
                for row in selected
            ]
        )
        external_fraction = np.array(
            [float(row["column_external_fraction_of_void"]) for row in selected]
        )
        axes[1].scatter(
            external_fraction,
            pathway + (rule_index - 0.5) * 0.04,
            color=colors[rule],
            marker=markers[rule],
            alpha=0.75,
            label=rule,
        )
    axes[0].set(
        xticks=range(len(labels)),
        xticklabels=labels,
        ylabel="Void fraction within analysis mask [-]",
        ylim=(0.0, 1.0),
    )
    axes[0].grid(axis="y", color=GRID_COLOR, alpha=0.18)
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].set(
        xlabel="Column-mask void connected to exterior air [-]",
        ylabel="Continuous bottom-to-top/lateral path [0/1]",
        yticks=(0, 1),
        ylim=(-0.15, 1.15),
        xlim=(0.0, 1.02),
    )
    axes[1].grid(color=GRID_COLOR, alpha=0.18)
    axes[1].legend(frameon=False, fontsize=7)
    figure.suptitle("Envelope choice and exterior-air access (Phase 2A.1 masks)")
    _save_rgb(figure, output_path)
    plt.close(figure)
