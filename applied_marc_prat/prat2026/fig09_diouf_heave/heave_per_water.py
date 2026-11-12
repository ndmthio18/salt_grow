"""The U_evap-free test of the growth law on the Diouf surface-heave experiment.

Figure 9b compares heave against TIME, which needs U_evap -- and re-deriving
U_evap from Figure 9a gives 5.38 mm/day against the 5.1 the paper quotes, so a
5% ambiguity sits under that comparison.  Plotting heave against EVAPORATED
WATER removes it: eq (32) says

    dh/dt = G_e * j       and      dM_w/dt = j * A
    =>  dh/dM_w = G_e / A          with  G_e = C_eq F(Da*) / (rho_cr (1-eps_s)(1-C_eq))

which contains neither U_evap nor the back-calibrated brine density rho = 1200.
It is the same U-invariant construction used on Norouzi Rad's Figure 4(a) in
new_model/validation, applied to the second, independent experiment.

Why this case is the cleaner of the two:  Diouf's salt grows INSIDE the bead
packing and lifts the surface, so there is no structure protruding into the
vapour field, no flux concentration and no bare surface competing for the
evaporation.  The measured h/M_w is therefore constant from the first point
(5.57 / 5.37 / 5.61 / 5.63 mm/g, a 4.6% spread) instead of drifting upward the
way Norouzi Rad's ratio does.  So the coefficient can be read off directly.

Result: with the paper's BASE porosity the law over-predicts by ~20%; with the
porosity the paper fits for this figure (eps_s = 0.34) it over-predicts by 6%.

NOTHING HERE IS A SIMULATION, and the labels say so.  The two straight lines are
eq (32) evaluated as algebra -- one coefficient, G_e/A, times M_w.  There is no
3D run because there is no geometry to run: Diouf's salt precipitates INSIDE the
bead packing, so the moving vapour boundary that the efflorescence3d model exists
to resolve simply is not present.  Running it anyway would return these same
lines: the 1D anchor in new_model/tests/test_efflorescence3d.py verifies that a
flat salt layer advances at core.growth_rate(j) to 1e-12.  That degeneracy is
exactly what makes this case a clean test of the CLOSURE COEFFICIENT alone.

Run from prat2026/fig09_diouf_heave/:
    ../../.venv/bin/python heave_per_water.py
"""

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import core                                    # noqa: E402
import params as pp                            # noqa: E402
from params import DioufHeave as D             # noqa: E402

KR45 = core.kr_at_temperature(pp.T45_K)
AREA = np.pi * (D.CYLINDER_DIAMETER / 2.0) ** 2
DATA = HERE.parent / "data"


def load(name):
    d = np.loadtxt(DATA / name, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def predicted(eps_s):
    """dh/dM_w [m/kg] from eq (32) -- no U_evap, no rho."""
    ds = core.Da_star(eps_s=eps_s, d_sub=D.D_SUB, kr=KR45, Ds=D.DS)
    F = core.F(ds)
    G_e = pp.C_EQ * F / (pp.RHO_CR * (1.0 - eps_s) * (1.0 - pp.C_EQ))
    return G_e / AREA, ds, F


def main():
    t_w, M_w = load("fig09a_exp.csv")
    t_h, h = load("fig09b_exp.csv")

    # Fig 9b's last marker sits just past the end of the Fig 9a trace, so the
    # mass curve is linearly extended over its last third to cover it.
    tail = t_w > 25.0
    a1, b1 = np.polyfit(t_w[tail], M_w[tail], 1)
    Mw_at_h = np.where(t_h <= t_w.max(), np.interp(t_h, t_w, M_w), a1 * t_h + b1)

    (s, b), C = np.polyfit(Mw_at_h, h, 1, cov=True)
    se = np.sqrt(C[0, 0])
    s0 = float(np.sum(Mw_at_h * h) / np.sum(Mw_at_h ** 2))   # through origin

    print("Diouf surface heave -- U-invariant test of eq (32)")
    print(f"  vessel 7 mm -> A = {AREA*1e6:.2f} mm^2, "
          f"T = 45 C -> k_r = {KR45:.3e} m/s, D_s = {D.DS:.1e} m^2/s")
    print(f"\n  {'point':>5} {'t [hr]':>8} {'M_w [g]':>9} {'h [mm]':>8} "
          f"{'h/M_w':>8}")
    for i in range(len(t_h)):
        print(f"  {i:>5} {t_h[i]:>8.2f} {Mw_at_h[i]:>9.4f} {h[i]:>8.4f} "
              f"{h[i]/Mw_at_h[i]:>8.3f}")
    rat = h / Mw_at_h
    print(f"  ratio spread {(rat.max()-rat.min())/rat.mean()*100:.1f}% "
          f"-> h is proportional to M_w, the law's FORM is confirmed")

    print(f"\n  measured  free intercept  {s:.3f} +/- {se:.3f} mm/g "
          f"(intercept {b:+.4f} mm, i.e. zero)")
    print(f"  measured  through origin   {s0:.3f} mm/g   "
          f"(the paper's model lines also pass through the origin)")

    rows = []
    for eps_s, tag in ((0.43, "base value"), (0.34, "fitted for this figure")):
        v, ds, F = predicted(eps_s)
        rows.append((eps_s, tag, v, ds, F))
        print(f"\n  predicted eps_s = {eps_s} ({tag})")
        print(f"    Da* = {ds:.1f} (paper text "
              f"{D.DA_STAR_PAPER_BASE if eps_s == 0.43 else D.DA_STAR_PAPER_FITTED}), "
              f"F = {F:.4f}")
        print(f"    dh/dM_w = {v:.3f} mm/g   ->  {v/s0-1:+.1%} vs measured")
    v1, _, _ = predicted(0.34)
    ceil = pp.C_EQ * 1.0 / (pp.RHO_CR * (1 - 0.34) * (1 - pp.C_EQ)) / AREA
    print(f"\n  F = 1 ceiling at eps_s = 0.34: {ceil:.3f} mm/g "
          f"({ceil/s0-1:+.1%}) -- so F cannot absorb the gap")

    # ------------------------------------------------------------------ #
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    xx = np.linspace(0, 0.40, 20)
    ax.plot(Mw_at_h, h, "o", color="0.15", ms=9,
            label="DATA  Diouf heave, Fig. 9b (4 markers)")
    ax.plot(xx, s0 * xx, ":", color="0.15", lw=2,
            label=f"DATA  fit through origin: {s0:.2f} mm/g")
    for (eps_s, tag, v, ds, F), c in zip(rows, ("#d62728", "#1f77b4")):
        ax.plot(xx, v * xx, "-", color=c, lw=2.2,
                label=f"ALGEBRA  eq (32), $\\varepsilon_s$={eps_s} ({tag})\n"
                      f"    {v:.2f} mm/g, {v/s0-1:+.0%}")
    ax.set_xlim(0, 0.40)
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("cumulative evaporated water  $M_w$  [g]")
    ax.set_ylabel("surface heave  $h$  [mm]")
    ax.set_title("(a) heave per evaporated water — no $U_{evap}$, no $\\rho$")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    ax2.axhline(s0, color="0.15", ls=":", lw=2, label=f"measured {s0:.2f} mm/g")
    ax2.axhspan(s0 - se, s0 + se, color="0.6", alpha=0.25)
    ee = np.linspace(0.28, 0.48, 60)
    ax2.plot(ee, [predicted(e)[0] for e in ee], "-", color="#9467bd", lw=2.5,
             label="ALGEBRA  eq (32) vs assumed $\\varepsilon_s$")
    for (eps_s, tag, v, ds, F), c in zip(rows, ("#d62728", "#1f77b4")):
        ax2.plot([eps_s], [v], "o", color=c, ms=10,
                 label=f"$\\varepsilon_s$={eps_s} ({tag})")
    # where does eq (32) meet the measurement?
    from scipy.optimize import brentq
    eps_match = brentq(lambda e: predicted(e)[0] - s0, 0.28, 0.48)
    ax2.axvline(eps_match, color="#2ca02c", ls="--", lw=1.8,
                label=f"matches measurement at $\\varepsilon_s$={eps_match:.3f}")
    ax2.set_xlabel(r"assumed structure porosity  $\varepsilon_s$")
    ax2.set_ylabel("$dh/dM_w$  [mm/g]")
    ax2.set_title("(b) the gap is absorbed by the fitted porosity")
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(alpha=0.3)
    print(f"\n  eq (32) reproduces the measured slope at eps_s = {eps_match:.3f}")

    fig.suptitle("Diouf (2018) surface heave — no simulation here: eq (32) is "
                 "ALGEBRA (one coefficient), and it is ~20% high at base "
                 "porosity", fontsize=10.5)
    fig.tight_layout(rect=(0, 0.01, 1, 0.94))
    out = HERE / "heave_per_water.png"
    fig.savefig(out, dpi=150)
    print(f"\nfigure: {out}")


if __name__ == "__main__":
    sys.exit(main())
