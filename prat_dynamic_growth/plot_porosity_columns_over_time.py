"""Plot physical-height porosity columns for the archived k_r study.

The two completed 6 h cases contain the requested quarter-time snapshots in
their archived NPZ files.  The two porosity-limited cases contain only their
initial and terminal states, so this script re-runs those cases in memory with
the same archived B2 ALE model/configuration and samples at 0, 25, 50, 75 and
100 percent of the archived terminal time.  No archived result is modified.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from math import sqrt
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


PROJECT_DIR = Path(__file__).resolve().parent
SRC_DIR = PROJECT_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from prat_growth.dynamic import DynamicALEModel, DynamicResult  # noqa: E402
from prat_growth.parameters import NumericalParameters, PhysicalParameters  # noqa: E402


FRACTIONS = np.array([0.0, 0.25, 0.50, 0.75, 1.0])
CASE_DEFINITIONS = (
    ("S3_kr2p3e6_power2", 2.3e-6),
    ("S3_kr2p3e5_power2", 2.3e-5),
    ("S3_kr2p3e4_power2", 2.3e-4),
    ("baseline_B2", 2.3e-3),
)


@dataclass(frozen=True)
class SnapshotSet:
    case_id: str
    k_reaction: float
    time: np.ndarray
    height: np.ndarray
    epsilon: np.ndarray
    faces: np.ndarray
    status: int
    message: str
    source: str


class ArchivedM1DynamicALEModel(DynamicALEModel):
    """Exact M1 flux closure used to create ``outputs/runs``.

    The archived sensitivity files were committed with the original solver at
    git revision 9d923ea.  M3 later corrected the common ALE face porosity and
    the top-boundary quadratic.  Only the two early-stopped cases need missing
    intermediate states reconstructed, so this local subclass preserves those
    two historical M1 methods without changing the current production solver.
    """

    SOURCE_REVISION = "9d923ea"

    def top_state(
        self,
        concentration_last: float,
        epsilon_last: float,
        d_eff_last: float,
        height: float,
    ) -> tuple[float, float]:
        p = self.p
        distance = height * (self.faces[-1] - self.centers[-1])
        conductance = p.rho_solution * epsilon_last * d_eff_last / distance
        kinetic_mass = p.beta * p.rho_solution * (1.0 - p.epsilon_birth) * p.k_reaction
        qa = kinetic_mass
        qb = p.evaporation_mass_flux - kinetic_mass * (1.0 + p.C_eq) - conductance
        qc = kinetic_mass * p.C_eq + conductance * concentration_last
        discriminant = qb * qb - 4.0 * qa * qc
        if discriminant < 0.0:
            raise RuntimeError(f"negative top-boundary discriminant: {discriminant}")
        root = sqrt(discriminant)
        candidates = np.array([(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)])
        physical = candidates[(candidates >= 0.0) & (candidates <= 1.0)]
        if physical.size:
            c_top = float(physical[np.argmin(np.abs(physical - p.C_eq))])
        else:
            c_top = float(candidates[np.argmin(np.abs(candidates - p.C_eq))])
        velocity = p.beta * p.rho_solution / p.rho_crystal * p.k_reaction * (c_top - p.C_eq)
        return c_top, velocity

    def flux_state(self, epsilon: np.ndarray, concentration: np.ndarray, height: float):
        p = self.p
        n = self.n.n_cells
        volumes = height * self.widths
        a_v, d_eff = self.coefficients(epsilon)
        reaction = p.rho_solution * a_v * p.k_reaction * (concentration - p.C_eq)
        c_top, velocity = self.top_state(concentration[-1], epsilon[-1], d_eff[-1], height)

        liquid_flux = np.empty(n + 1)
        structural_density_birth = (
            p.rho_solution * p.epsilon_birth
            + p.rho_crystal * (1.0 - p.epsilon_birth)
        )
        liquid_flux[-1] = structural_density_birth * velocity + p.evaporation_mass_flux
        density_factor = 1.0 - p.rho_solution / p.rho_crystal
        for index in range(n - 1, -1, -1):
            liquid_flux[index] = (
                liquid_flux[index + 1]
                + density_factor * reaction[index] * volumes[index]
            )

        mesh_velocity = self.faces * velocity
        epsilon_face = np.empty(n + 1)
        epsilon_face[0] = epsilon[0]
        epsilon_face[-1] = p.epsilon_birth
        epsilon_face[1:-1] = 0.5 * (epsilon[:-1] + epsilon[1:])

        salt_flux = np.empty(n + 1)
        bottom_distance = height * (self.centers[0] - self.faces[0])
        bottom_diffusive = (
            -p.rho_solution
            * epsilon[0]
            * d_eff[0]
            * (concentration[0] - p.C_eq)
            / bottom_distance
        )
        salt_flux[0] = p.C_eq * liquid_flux[0] + bottom_diffusive

        for face in range(1, n):
            left = face - 1
            right = face
            distance = height * (self.centers[right] - self.centers[left])
            k_left = epsilon[left] * d_eff[left]
            k_right = epsilon[right] * d_eff[right]
            half_left = height * (self.faces[face] - self.centers[left])
            half_right = height * (self.centers[right] - self.faces[face])
            k_face = distance / (half_left / k_left + half_right / k_right)
            diffusive = (
                -p.rho_solution
                * k_face
                * (concentration[right] - concentration[left])
                / distance
            )
            relative_liquid = (
                liquid_flux[face]
                - p.rho_solution * epsilon_face[face] * mesh_velocity[face]
            )
            advected_c = concentration[left] if relative_liquid >= 0.0 else concentration[right]
            salt_flux[face] = relative_liquid * advected_c + diffusive

        salt_flux[-1] = p.rho_crystal * (1.0 - p.epsilon_birth) * velocity
        epsilon_mesh_flux = np.empty(n + 1)
        epsilon_mesh_flux[0] = 0.0
        epsilon_mesh_flux[1:-1] = -mesh_velocity[1:-1] * epsilon[1:]
        epsilon_mesh_flux[-1] = -velocity * p.epsilon_birth
        return {
            "volumes": volumes,
            "a_v": a_v,
            "d_eff": d_eff,
            "reaction": reaction,
            "c_top": c_top,
            "velocity": velocity,
            "liquid_flux": liquid_flux,
            "salt_flux": salt_flux,
            "epsilon_mesh_flux": epsilon_mesh_flux,
        }


def _load_result(runs_dir: Path, case_id: str) -> DynamicResult:
    metadata = json.loads((runs_dir / f"{case_id}.json").read_text(encoding="utf-8"))
    with np.load(runs_dir / f"{case_id}.npz") as arrays:
        return DynamicResult(
            closure=metadata["closure"],
            physical=PhysicalParameters(**metadata["physical"]),
            numerical=NumericalParameters(**metadata["numerical"]),
            faces=arrays["faces"].copy(),
            centers=arrays["centers"].copy(),
            widths=arrays["widths"].copy(),
            time=arrays["time"].copy(),
            epsilon=arrays["epsilon"].copy(),
            concentration=arrays["concentration"].copy(),
            height=arrays["height"].copy(),
            top_solid=arrays["top_solid"].copy(),
            internal_solid=arrays["internal_solid"].copy(),
            bottom_salt=arrays["bottom_salt"].copy(),
            bottom_solution=arrays["bottom_solution"].copy(),
            evaporated=arrays["evaporated"].copy(),
            status=int(metadata["solver_status"]),
            message=str(metadata["solver_message"]),
            initial_profile=metadata.get("initial_profile", "equilibrium"),
            solver_warnings=metadata.get("solver_warnings", []),
        )


def _matching_indices(time: np.ndarray, targets: np.ndarray) -> list[int] | None:
    indices: list[int] = []
    # Event localization can differ by sub-microseconds across SciPy builds;
    # this tolerance is still at least eight orders below the simulated span.
    tolerance = max(1e-8, 1e-8 * max(1.0, float(targets[-1])))
    for target in targets:
        index = int(np.argmin(np.abs(time - target)))
        if abs(float(time[index] - target)) > tolerance:
            return None
        indices.append(index)
    return indices


def _sample_case(runs_dir: Path, case_id: str, expected_k: float) -> SnapshotSet:
    archived = _load_result(runs_dir, case_id)
    if not np.isclose(archived.physical.k_reaction, expected_k, rtol=1e-12, atol=0.0):
        raise ValueError(
            f"{case_id}: archived k_r={archived.physical.k_reaction:g}, "
            f"expected {expected_k:g}"
        )

    archived_final = float(archived.time[-1])
    targets = FRACTIONS * archived_final
    indices = _matching_indices(archived.time, targets)
    source = "archived"
    sampled = archived

    if indices is None:
        # t_eval does not affect solve_ivp's adaptive integration steps.  Using
        # the archived terminal time therefore reproduces the stopped case while
        # exposing the requested quarter-time states without writing a new run.
        numerical = archived.numerical.with_updates(
            t_end=archived_final,
            output_interval=archived_final / 4.0,
        )
        sampled = ArchivedM1DynamicALEModel(
            archived.physical, numerical, archived.closure
        ).run(archived.initial_profile)
        indices = _matching_indices(sampled.time, targets)
        if indices is None:
            raise RuntimeError(
                f"{case_id}: ALE reconstruction did not return the requested times; "
                f"got {sampled.time}, expected {targets}"
            )
        if not np.isclose(sampled.height[-1], archived.height[-1], rtol=2e-8, atol=1e-12):
            raise RuntimeError(f"{case_id}: reconstructed terminal height differs from archive")
        if not np.allclose(
            sampled.epsilon[:, -1], archived.epsilon[:, -1], rtol=2e-7, atol=2e-9
        ):
            raise RuntimeError(f"{case_id}: reconstructed terminal porosity differs from archive")
        source = (
            "ALE reconstruction (archived M1 equations/configuration, "
            f"revision {ArchivedM1DynamicALEModel.SOURCE_REVISION})"
        )

    assert indices is not None
    selected_time = sampled.time[indices]
    selected_height = sampled.height[indices]
    selected_epsilon = sampled.epsilon[:, indices]
    _validate_snapshots(
        case_id,
        targets,
        selected_time,
        selected_height,
        selected_epsilon,
        sampled.faces,
    )
    return SnapshotSet(
        case_id=case_id,
        k_reaction=expected_k,
        time=selected_time,
        height=selected_height,
        epsilon=selected_epsilon,
        faces=sampled.faces,
        status=archived.status,
        message=archived.message,
        source=source,
    )


def _validate_snapshots(
    case_id: str,
    targets: np.ndarray,
    time: np.ndarray,
    height: np.ndarray,
    epsilon: np.ndarray,
    faces: np.ndarray,
) -> None:
    if time.shape != (5,) or height.shape != (5,) or epsilon.shape[1] != 5:
        raise ValueError(f"{case_id}: expected exactly five snapshots")
    if epsilon.shape[0] + 1 != faces.size:
        raise ValueError(f"{case_id}: porosity cells and ALE faces are inconsistent")
    if not np.allclose(time, targets, rtol=0.0, atol=max(1e-8, targets[-1] * 1e-8)):
        raise ValueError(f"{case_id}: selected times do not match target fractions")
    if not np.all(np.diff(time) > 0.0):
        raise ValueError(f"{case_id}: selected times are not strictly increasing")
    if not np.isclose(faces[0], 0.0) or not np.isclose(faces[-1], 1.0):
        raise ValueError(f"{case_id}: ALE faces do not span [0, 1]")
    if np.any(height <= 0.0) or np.any(~np.isfinite(epsilon)):
        raise ValueError(f"{case_id}: non-physical or non-finite snapshot data")
    if np.any((epsilon < 0.0) | (epsilon > 1.0)):
        raise ValueError(f"{case_id}: porosity lies outside [0, 1]")
    # The upper physical face must be exactly the plotted layer height.
    for column_height in height:
        if not np.isclose(column_height * faces[-1], column_height, rtol=0.0, atol=1e-15):
            raise ValueError(f"{case_id}: physical face mapping does not end at h(t)")


def _scientific_k(value: float) -> str:
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / 10.0**exponent
    return rf"{coefficient:.1f}$\times$10$^{{{exponent}}}$"


def _time_label(value: float) -> str:
    if np.isclose(value, round(value), atol=0.05):
        return f"{round(value):d}"
    return f"{value:.1f}"


def _terminal_title(case: SnapshotSet) -> str:
    final_time = float(case.time[-1])
    if case.status == 1:
        epsilon_min = float(np.min(case.epsilon[:, -1]))
        return f"stopped at {_time_label(final_time)} s ($\\epsilon_{{min}}$ = {epsilon_min:.2f})"
    return f"final time = {_time_label(final_time)} s"


def _draw_panel(
    ax: plt.Axes,
    case: SnapshotSet,
    norm: mpl.colors.Normalize,
    cmap: mpl.colors.Colormap,
    y_limit_mm: float,
    panel_label: str | None = None,
) -> mpl.collections.QuadMesh:
    final_time = float(case.time[-1])
    column_width = 0.075 * final_time
    mesh: mpl.collections.QuadMesh | None = None

    for index, (time, height, epsilon) in enumerate(
        zip(case.time, case.height, case.epsilon.T, strict=True)
    ):
        x_faces = np.array([time - column_width / 2.0, time + column_width / 2.0])
        z_faces_mm = 1e3 * height * case.faces
        x_grid, z_grid = np.meshgrid(x_faces, z_faces_mm)
        mesh = ax.pcolormesh(
            x_grid,
            z_grid,
            epsilon[:, None],
            cmap=cmap,
            norm=norm,
            shading="flat",
            edgecolors="none",
            rasterized=True,
            zorder=2,
        )
        ax.add_patch(
            Rectangle(
                (x_faces[0], 0.0),
                column_width,
                1e3 * height,
                fill=False,
                edgecolor="0.45",
                linewidth=0.65,
                zorder=3,
            )
        )
        # Verify that the artist's frame encodes the physical h(t), not a
        # normalized height.  This is intentionally checked during plotting.
        frame_height = ax.patches[-1].get_height()
        if not np.isclose(frame_height, 1e3 * height, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"{case.case_id}: plotted frame height mismatch at index {index}")

    if mesh is None:
        raise RuntimeError(f"{case.case_id}: no columns were drawn")

    margin = 0.09 * final_time
    ax.set_xlim(-margin, final_time + margin)
    ax.set_ylim(0.0, y_limit_mm)
    ax.set_xticks(case.time)
    ax.set_xticklabels([_time_label(value) for value in case.time], rotation=25, ha="right")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Salt layer height $z$ (mm)")
    ax.set_title(
        rf"$k_r$ = {_scientific_k(case.k_reaction)} m s$^{{-1}}$" + "\n" + _terminal_title(case),
        pad=7,
    )
    ax.set_facecolor("white")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.90", linewidth=0.55)
    ax.xaxis.grid(False)
    for spine in ax.spines.values():
        spine.set_color("0.35")
        spine.set_linewidth(0.7)
    if panel_label is not None:
        ax.text(
            -0.12,
            1.08,
            panel_label,
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            ha="left",
            va="top",
        )
    return mesh


def _configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _save_figure(fig: plt.Figure, output_path: Path, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, facecolor="white", bbox_inches=None)
    plt.close(fig)


def _label_colorbar(
    colorbar: mpl.colorbar.Colorbar, norm: mpl.colors.Normalize
) -> None:
    ticks = np.linspace(float(norm.vmin), float(norm.vmax), 6)
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([f"{tick:.2f}" for tick in ticks])
    colorbar.set_label("Porosity $\\epsilon$ (-)")


def _make_combined_figure(
    cases: list[SnapshotSet],
    output_path: Path,
    norm: mpl.colors.Normalize,
    cmap: mpl.colors.Colormap,
    y_limit_mm: float,
    dpi: int,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), layout="constrained", sharey=True)
    meshes = []
    for ax, case, panel_label in zip(axes.flat, cases, ("(a)", "(b)", "(c)", "(d)"), strict=True):
        meshes.append(_draw_panel(ax, case, norm, cmap, y_limit_mm, panel_label))
    colorbar = fig.colorbar(
        meshes[-1],
        ax=axes,
        orientation="vertical",
        shrink=0.92,
        pad=0.025,
        aspect=30,
    )
    _label_colorbar(colorbar, norm)
    _save_figure(fig, output_path, dpi)


def _make_single_figures(
    cases: list[SnapshotSet],
    output_dir: Path,
    norm: mpl.colors.Normalize,
    cmap: mpl.colors.Colormap,
    y_limit_mm: float,
    dpi: int,
) -> list[Path]:
    paths: list[Path] = []
    for case in cases:
        exponent = int(round(-np.log10(case.k_reaction / 2.3)))
        path = output_dir / f"saltgrow_porosity_columns_kr_2p3em{exponent}.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.4), layout="constrained")
        mesh = _draw_panel(ax, case, norm, cmap, y_limit_mm)
        colorbar = fig.colorbar(mesh, ax=ax, pad=0.025, aspect=28)
        _label_colorbar(colorbar, norm)
        _save_figure(fig, path, dpi)
        paths.append(path)
    return paths


def _write_audit(
    cases: list[SnapshotSet],
    output_dir: Path,
    y_limit_mm: float,
    color_min: float,
    color_max: float,
    combined_path: Path,
    single_paths: list[Path],
) -> tuple[Path, Path]:
    csv_path = output_dir / "saltgrow_porosity_columns_snapshot_audit.csv"
    rows: list[dict[str, object]] = []
    for case in cases:
        for fraction, time, height, epsilon in zip(
            FRACTIONS, case.time, case.height, case.epsilon.T, strict=True
        ):
            rows.append(
                {
                    "case_id": case.case_id,
                    "k_reaction_m_s": case.k_reaction,
                    "source": case.source,
                    "fraction_of_final_time": fraction,
                    "time_s": time,
                    "height_m": height,
                    "height_mm": 1e3 * height,
                    "epsilon_min": float(np.min(epsilon)),
                    "epsilon_max": float(np.max(epsilon)),
                    "n_cells": epsilon.size,
                    "column_top_matches_height": bool(
                        np.isclose(1e3 * height * case.faces[-1], 1e3 * height)
                    ),
                }
            )
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    manifest_path = output_dir / "saltgrow_porosity_columns_manifest.json"
    manifest = {
        "figure_type": "physical-height column snapshots with internal porosity cloud map",
        "fractions_of_final_time": FRACTIONS.tolist(),
        "common_y_limit_mm": y_limit_mm,
        "shared_color_limits": {"epsilon_min": color_min, "epsilon_max": color_max},
        "colormap": "viridis",
        "combined_figure": combined_path.name,
        "single_figures": [path.name for path in single_paths],
        "cases": [
            {
                "case_id": case.case_id,
                "k_reaction_m_s": case.k_reaction,
                "source": case.source,
                "archived_status": case.status,
                "archived_message": case.message,
                "selected_times_s": case.time.tolist(),
                "selected_heights_mm": (1e3 * case.height).tolist(),
            }
            for case in cases
        ],
        "quality_checks": {
            "five_times_per_case": True,
            "time_positions_match_requested_fractions": True,
            "physical_column_tops_match_h_t": True,
            "porosity_arrays_match_selected_result_fields": True,
            "blank_region_above_each_column": True,
            "no_porosity_line_plots": True,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return csv_path, manifest_path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=PROJECT_DIR / "outputs" / "runs",
        help="Directory containing archived NPZ/JSON case pairs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "outputs" / "figures",
        help="Directory for the new figures and audit files.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    _configure_style()
    cases = [_sample_case(args.runs_dir, case_id, k_value) for case_id, k_value in CASE_DEFINITIONS]

    all_epsilon = np.concatenate([case.epsilon.ravel() for case in cases])
    color_min = float(np.min(all_epsilon))
    color_max = float(np.max(all_epsilon))
    norm = mpl.colors.Normalize(vmin=color_min, vmax=color_max)
    cmap = mpl.colormaps["viridis"].with_extremes(bad="#b3b3b3")

    maximum_height_mm = max(float(np.max(case.height)) for case in cases) * 1e3
    y_limit_mm = 1.05 * maximum_height_mm
    if y_limit_mm <= maximum_height_mm:
        raise RuntimeError("common y-limit does not leave blank space above the tallest column")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_path = args.output_dir / "saltgrow_porosity_columns_2x2.png"
    _make_combined_figure(cases, combined_path, norm, cmap, y_limit_mm, args.dpi)
    single_paths = _make_single_figures(cases, args.output_dir, norm, cmap, y_limit_mm, args.dpi)
    csv_path, manifest_path = _write_audit(
        cases,
        args.output_dir,
        y_limit_mm,
        color_min,
        color_max,
        combined_path,
        single_paths,
    )

    print(f"Wrote {combined_path}")
    for path in single_paths:
        print(f"Wrote {path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
