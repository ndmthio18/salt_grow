from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from .analytic import general_growth_velocity, large_da_growth_velocity
from .dynamic import DynamicResult, snapshot_diagnostics


def _save(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def _model_a_height(result: DynamicResult) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = result.physical
    t = result.time
    sol = solve_ivp(
        lambda _t, h: [general_growth_velocity(p, float(h[0]))],
        (float(t[0]), float(t[-1])),
        [p.h0],
        t_eval=t,
        rtol=1e-11,
        atol=1e-14,
    )
    height = sol.y[0]
    velocity = np.array([general_growth_velocity(p, h) for h in height])
    return t, height, velocity


def make_required_figures(results: dict[str, DynamicResult], output: Path) -> None:
    if "B2" not in results:
        raise ValueError("B2 baseline is required for the figure package")
    b2 = results["B2"]
    t_a, h_a, v_a = _model_a_height(b2)
    time_h = b2.time / 3600.0

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(t_a / 3600.0, 1e3 * h_a, "k--", label="Model A (general analytical)")
    for label, result in results.items():
        ax.plot(result.time / 3600.0, 1e3 * result.height, label=f"Model {label}")
    ax.set(xlabel="Time (h)", ylabel="Height (mm)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    _save(fig, output, "figure_1_height_evolution")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(t_a / 3600.0, 86400.0 * v_a * 1e3, "k--", label="Model A")
    for label, result in results.items():
        velocity = [snapshot_diagnostics(result, i)["velocity"] for i in range(result.time.size)]
        ax.plot(result.time / 3600.0, np.asarray(velocity) * 86400.0 * 1e3,
                label=f"Model {label}")
    ax.set(xlabel="Time (h)", ylabel="Growth velocity (mm day$^{-1}$)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    _save(fig, output, "figure_2_growth_velocity")

    indices = sorted(set([0, max(1, b2.time.size // 6), b2.time.size // 2, b2.time.size - 1]))
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for i in indices:
        diag = snapshot_diagnostics(b2, i)
        distance = b2.height[i] * (1.0 - b2.centers)
        reaction = np.maximum(np.asarray(diag["reaction"]), 0.0)
        ax.semilogy(1e6 * distance, np.maximum(reaction, 1e-12),
                    label=f"{b2.time[i] / 3600:.2g} h")
    final = snapshot_diagnostics(b2, -1)
    x_limit = max(10.0e-6, 6.0 * float(final["delta_C"]))
    distance_ref = np.linspace(0.0, x_limit, 300)
    r0 = max(float(np.asarray(final["reaction"])[-1]), 1e-12)
    ax.semilogy(1e6 * distance_ref, r0 * np.exp(-distance_ref / float(final["delta_C"])),
                "k--", label="exp(-d/delta_C), final")
    ax.set(xlabel="Distance below moving interface (micrometres)",
           ylabel="Positive internal precipitation (kg m$^{-3}$ s$^{-1}$)")
    ax.set_xlim(left=0.0, right=1e6 * x_limit)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    _save(fig, output, "figure_3_internal_precipitation")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for i in indices:
        z = b2.height[i] * b2.centers
        ax.plot(1e3 * z, b2.epsilon[:, i], label=f"{b2.time[i] / 3600:.2g} h")
    ax.axvline(1e3 * b2.physical.h0, color="0.4", linestyle=":",
               label="initial top (old/new layer boundary)")
    ax.set(xlabel="z (mm)", ylabel="Porosity", ylim=(0.0, 0.5))
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    _save(fig, output, "figure_4_porosity_evolution")

    dissolved = []
    residual = []
    for i in range(b2.time.size):
        diag = snapshot_diagnostics(b2, i)
        dissolved.append(diag["dissolved_salt"])
        residual.append(diag["salt_residual"])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(time_h, b2.top_solid, label="new-layer solid (cumulative)")
    ax.plot(time_h, b2.internal_solid, label="internal solid (cumulative)")
    ax.plot(time_h, dissolved, label="dissolved salt inventory")
    ax.set(xlabel="Time (h)", ylabel="Salt mass per area (kg m$^{-2}$)")
    ax2 = ax.twinx()
    ax2.plot(time_h, residual, "k:", label="salt conservation residual")
    ax2.set_ylabel("Residual (kg m$^{-2}$)")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], frameon=False,
              loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2)
    ax.grid(alpha=0.2)
    _save(fig, output, "figure_5_salt_mass_partition")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    v_ref = 0.01 * large_da_growth_velocity(b2.physical)
    for label, result in results.items():
        vb = np.array([snapshot_diagnostics(result, i)["velocity"]
                       for i in range(result.time.size)])
        va = np.array([general_growth_velocity(result.physical, h) for h in result.height])
        error = np.abs(vb - va) / np.maximum(np.abs(vb), v_ref)
        ax.plot(result.time / 3600.0, 100.0 * error, label=f"Model {label}")
    ax.set(xlabel="Time (h)", ylabel="Analytical velocity error (%)")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    ax.text(0.02, 0.98, r"$v_{ref}=0.01v_{A,large-Da}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=9)
    _save(fig, output, "figure_6_analytical_growth_law_error")
