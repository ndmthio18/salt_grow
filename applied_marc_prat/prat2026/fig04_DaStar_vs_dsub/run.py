"""Figure 4 -- the height-independent Damkoehler number Da* vs grain size.

Da* is the one number that decides what fraction of the evaporation-supplied
salt becomes front advance, so this figure is the paper's parameter map. Eq (34)
makes the linearity explicit: Da* is proportional to d_sub.

Panel layout follows the paper: a linear main panel over [0, 1000] um, a left
inset showing F = sqrt(Da*)/(1+sqrt(Da*)) on a log abscissa, and a right inset
zooming into [0, 100] um with a LOG ordinate spanning 0.1 to 100. That log
ordinate matters -- on the main panel's linear axis everything below ~50 um is
indistinguishable from zero, so the right inset is the only place the small-grain
behaviour can actually be read, and it is the sharpest check on the published
curve.

KNOWN DISCREPANCY: our Da* sits ~24% above the curve the paper draws. The
published main panel and its right inset agree with each other to 2%, so the
figure is internally consistent; it is our formula -- eq (27)/(34) evaluated
with the paper's stated parameters -- that differs from it. See investigate.py.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import DioufHeave, NorouziRad

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"


def load(name):
    path = DATA / name
    if not path.exists():
        return None
    d = np.loadtxt(path, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def main():
    pub_main = load("fig04_model_curve.csv")
    pub_F = load("fig04_inset_F.csv")
    pub_zoom = load("fig04_inset_zoom.csv")

    fig, ax = plt.subplots(figsize=(6.8, 5.4))

    # --- main panel --------------------------------------------------------
    d_um = np.linspace(0.1, 1000.0, 800)
    if pub_main is not None:
        ax.plot(*pub_main, color="0.72", lw=4, label="published curve (digitised)")
    ax.plot(d_um, core.Da_star(d_sub=d_um * 1e-6), color="k", lw=1.5,
            label="eq (27)/(34), ours")

    ours_sand = core.Da_star(d_sub=NorouziRad.D_SUB)
    ax.plot([NorouziRad.D_SUB * 1e6], [NorouziRad.DA_STAR_PAPER], "v",
            color="crimson", ms=8, label="value quoted in the paper's text")
    ax.plot([NorouziRad.D_SUB * 1e6], [ours_sand], "o", mfc="none", color="k", ms=8)

    ax.set_xlim(0, 1000)
    ax.set_ylim(0, 300)
    ax.set_xlabel(r"$d_{sub}$ ($\mu$m)")
    ax.set_ylabel(r"$Da^*$")
    ax.set_title(r"Prat (2026) Figure 4 -- $Da^*$ from eq (27)/(34)", fontsize=11)
    # tucked under the left inset: the lower right is taken by the zoom inset
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(0.02, 0.545))

    if pub_main is not None:
        pub_at_ref = np.polyval(np.polyfit(*pub_main, 1), NorouziRad.D_SUB * 1e6)
        ax.annotate(
            f"at 0.58 mm\nours {ours_sand:.0f}\ntext {NorouziRad.DA_STAR_PAPER:.0f}\n"
            f"curve {pub_at_ref:.0f}",
            xy=(NorouziRad.D_SUB * 1e6, ours_sand), xytext=(330.0, 18.0),
            fontsize=8, arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
        )

    # --- left inset: kinetic efficiency F, log x ---------------------------
    ax_f = ax.inset_axes([0.11, 0.67, 0.31, 0.27])
    d_log = np.logspace(0, 3, 400)
    if pub_F is not None:
        ax_f.semilogx(*pub_F, color="0.72", lw=3)
    ax_f.semilogx(d_log, core.F(core.Da_star(d_sub=d_log * 1e-6)), color="k", lw=1.2)
    ax_f.set_xlim(1, 1000)
    ax_f.set_ylim(0, 1)
    ax_f.set_ylabel(r"$\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}$", fontsize=9)
    ax_f.set_xlabel(r"$d_{sub}$ ($\mu$m)", fontsize=7)
    ax_f.tick_params(labelsize=6)

    # --- right inset: zoom on [1, 100] um, LOG ordinate --------------------
    ax_zoom = ax.inset_axes([0.60, 0.13, 0.36, 0.31])
    d_zoom = np.linspace(0.5, 100.0, 400)
    if pub_zoom is not None:
        ax_zoom.semilogy(*pub_zoom, color="0.72", lw=3)
    ax_zoom.semilogy(d_zoom, core.Da_star(d_sub=d_zoom * 1e-6), color="k", lw=1.2)
    ax_zoom.set_xlim(0, 100)
    ax_zoom.set_ylim(0.1, 100)
    ax_zoom.set_xlabel(r"$d_{sub}$ ($\mu$m)", fontsize=7)
    ax_zoom.set_ylabel(r"$Da^*$", fontsize=8)
    ax_zoom.tick_params(labelsize=6)

    fig.tight_layout()
    fig.savefig(HERE / "fig04.png", dpi=180)
    plt.close(fig)

    # =====================================================================
    print("Figure 4 reproduction")
    print(f"  Da* is linear in d_sub: ours Da*/d_sub = "
          f"{core.Da_star(d_sub=1e-3) / 1e-3:.4g} 1/m")
    for d in (1, 10, 100, 580, 1000):
        ds = core.Da_star(d_sub=d * 1e-6)
        print(f"    d_sub = {d:5d} um -> Da* = {ds:8.2f}, F = {core.F(ds):.4f}")

    if pub_main is not None and pub_zoom is not None:
        s_main = np.polyfit(*pub_main, 1)[0] * 1e6
        zx, zy = pub_zoom
        keep = zx > 5.0
        s_zoom = float(np.sum(zx[keep] * zy[keep]) / np.sum(zx[keep] ** 2)) * 1e6
        s_ours = core.Da_star(d_sub=1e-3) / 1e-3
        print("\n  Da*/d_sub [1/m], from every source in the figure:")
        print(f"    published main panel   {s_main:.4g}")
        print(f"    published right inset  {s_zoom:.4g}   "
              f"(agrees with the main panel to {abs(s_zoom / s_main - 1) * 100:.1f}%)")
        print(f"    paper text, 156 @ 0.58 mm  "
              f"{NorouziRad.DA_STAR_PAPER / (NorouziRad.D_SUB):.4g}")
        print(f"    ours                   {s_ours:.4g}   "
              f"({s_ours / s_zoom:.3f}x the published figure)")

    print("\n  Da* the paper quotes elsewhere, against ours:")
    kr45 = core.kr_at_temperature(318.15)
    for label, ours, theirs in (
        ("sand, eps=0.43, 0.58 mm", core.Da_star(d_sub=NorouziRad.D_SUB), 156.0),
        ("heave, eps=0.43, 90 um, 45C",
         core.Da_star(eps_s=0.43, d_sub=DioufHeave.D_SUB, kr=kr45, Ds=DioufHeave.DS), 57.0),
        ("heave, eps=0.34, 90 um, 45C",
         core.Da_star(eps_s=0.34, d_sub=DioufHeave.D_SUB, kr=kr45, Ds=DioufHeave.DS), 92.0),
    ):
        print(f"    {label:30s} ours {ours:7.1f}  paper {theirs:5.0f}  "
              f"ratio {ours / theirs:.3f}  F differs by "
              f"{abs(core.F(ours) - core.F(theirs)) / core.F(theirs) * 100:.2f}%")


if __name__ == "__main__":
    main()
