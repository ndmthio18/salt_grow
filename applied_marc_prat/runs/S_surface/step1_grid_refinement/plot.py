"""Step 1B figures from results/step1_refinement.json (run.py first).

    $PY runs/S_surface/step1_grid_refinement/plot.py

fig1_integrated_flux_convergence.{png,svg}   20 vs 10 um: E_total, E_top, E_side, j_top_mean, j_side_mean (A1, A2)
fig2_side_effect_comparison.{png,svg}        the small difference dE = E_A2 - E_A1, 100 (R_sides - 1), f_side
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; FIG = HERE / "figures"


def save(fig, name):
    FIG.mkdir(exist_ok=True); fig.savefig(FIG / f"{name}.png", dpi=160); fig.savefig(FIG / f"{name}.svg"); plt.close(fig)


def fig1(r):
    panels = (("E_total", "E_total [kg/s]", 1.0), ("E_top", "E_top [kg/s]", 1.0), ("E_side", "E_side [kg/s] (A2)", 1.0),
              ("j_top_mean", r"j_top_mean = E_top / A_top [kg m$^{-2}$ s$^{-1}$]", 1.0),
              ("j_side_mean", r"j_side_mean = E_side / A_side [kg m$^{-2}$ s$^{-1}$] (A2)", 1.0))
    fig, ax = plt.subplots(1, 5, figsize=(17, 4.2))
    for k, (key, lab, sc) in enumerate(panels):
        a = ax[k]; x = 0
        for cond, c in (("A1", "tab:blue"), ("A2", "tab:red")):
            vals = [r["coarse"][cond][key], r["fine"][cond][key]]
            if vals[0] is None or vals[0] == 0.0:
                continue
            bars = a.bar([x, x + 1], vals, color=[c, c], alpha=[0.55, 1.0][0], width=0.8, edgecolor="k", lw=0.5)
            bars[1].set_alpha(1.0); bars[0].set_alpha(0.5)
            rel = (vals[1] - vals[0]) / vals[0]
            for b, v in zip(bars, vals):
                a.text(b.get_x() + b.get_width() / 2, b.get_height() * 0.5, f"{v:.6e}", ha="center", va="center", fontsize=6, rotation=90, color="w")
            a.text(x + 0.5, max(vals) * 1.32, f"{cond}: {rel:+.2e}", ha="center", fontsize=8, color=c)
            x += 2.6
        a.set_xticks([0, 1] + ([2.6, 3.6] if x > 3 else [])); a.set_xticklabels(["20 um", "10 um"] * (2 if x > 3 else 1), fontsize=8)
        a.set_ylabel(lab, fontsize=8); a.set_title(f"({'abcde'[k]}) {key}", fontsize=9)
        a.set_ylim(0, max(v for cond in ("A1", "A2") for v in (r["coarse"][cond][key] or 0, r["fine"][cond][key] or 0)) * 1.5)
    fig.suptitle("Fig. 1  Integrated fluxes at h = 0.3 mm, 20 um (light) vs 10 um (dark), A1 blue / A2 red; the annotated number is the "
                 "relative change (10 - 20) / 20.  Same physical geometry (540 um opening, 300 um height) on both grids.", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92)); save(fig, "fig1_integrated_flux_convergence")


def fig2(r):
    pc, pf = r["paired"]["coarse"], r["paired"]["fine"]
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 4.0))
    for k, (key, lab, scale, fmt) in enumerate((("dE", r"$\Delta E = E_{A2} - E_{A1}$ [kg/s]", 1.0, "%.4e"),
                                                ("R_sides_minus_1_pct", r"$100\,(R_{sides} - 1)$ [%]", 1.0, "%.4f"),
                                                ("f_side", r"$f_{side} = E_{side,A2} / E_{A2}$ [-]", 1.0, "%.5f"))):
        vals = [pc[key] * scale, pf[key] * scale]
        bars = ax[k].bar([0, 1], vals, color=["0.6", "tab:red"], width=0.6, edgecolor="k", lw=0.5)
        for b, v in zip(bars, vals):
            ax[k].text(b.get_x() + b.get_width() / 2, b.get_height(), fmt % v, ha="center", va="bottom", fontsize=8)
        ax[k].set_xticks([0, 1]); ax[k].set_xticklabels(["20 um", "10 um"]); ax[k].set_ylabel(lab)
        ax[k].set_ylim(0, max(vals) * 1.25)
        ax[k].set_title(f"({'abc'[k]})  change 20 -> 10 um: {(vals[1] - vals[0]) / vals[0]:+.2e} (rel)", fontsize=9)
    ax[0].text(0.5, max(pc["dE"], pf["dE"]) * 1.12, f"E_total ~ {r['fine']['A2']['E_total']:.3e} kg/s:\ndE is {pf['dE'] / r['fine']['A2']['E_total'] * 100:.3f} % of it",
               ha="center", fontsize=7.5)
    fig.suptitle("Fig. 2  The small side-face gain in TOTAL evaporation (a, b) against the large REDISTRIBUTION to the sides (c), h = 0.3 mm.\n"
                 "Axes start at zero.  The sides take ~27 % of the flux away from the top; the far-field resistance (20.6 mm) fixes the total, so the net gain is ~0.1 %.",
                 fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90)); save(fig, "fig2_side_effect_comparison")


def main():
    r = json.loads((RES / "step1_refinement.json").read_text()); fig1(r); fig2(r); print("figures written to", FIG)


if __name__ == "__main__":
    main()
