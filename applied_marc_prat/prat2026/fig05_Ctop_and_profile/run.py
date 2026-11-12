"""Figure 5 -- surface supersaturation and the concentration profile.

(a) C_top/C_eq vs substrate grain size, from eq (31). The gap from equilibrium
    is of order 1e-5: NaCl's k_r is so large that a minute supersaturation
    already sustains the observed growth.
(b) C/C_eq through the structure, from eq (19). Everything sits at equilibrium
    except a boundary layer at the moving surface -- which is why the paper does
    not need to discretise the interior, and why anyone who does needs an
    extremely fine mesh there.

THE UNSTATED HEIGHT. The paper gives U_evap, d_sub, k_r and eps_s for panel (b)
but never says what structure height h it used -- and on a z/h axis the apparent
layer thickness is 1/sqrt(Da), which scales as 1/h. A first attempt at h = 1 mm
therefore produced a layer ~8x too thick compared with the published panel, even
though every formula was right. Fitting the digitised published curve recovers
sqrt(Da) ~ 2.6e4, i.e. h = 8.4 +- 0.4 mm; that is used below and the published
curve is overlaid so the agreement is checkable rather than asserted.

The physics worth keeping straight: delta_C in METRES is set by reaction and
diffusion alone and does not depend on h at all (325 nm here). Only its value
RELATIVE to h changes, which is purely an artefact of plotting against z/h.

NUMERICAL NOTE: this panel is why core.conc_profile does not use eq (19)
literally. Even at h = 1 mm, sqrt(Da) ~ 3078 and sinh(3078) overflows to inf;
at the inferred h it is ~2.6e4 and the overflow is far worse. The exponential
rewrite in core keeps every exponent non-positive.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import C_EQ, DAY

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"

U_EVAP = 0.01 / DAY  # 1 cm/day, paper text
EPS_S = 0.43
D_SUB_PROFILE = 10e-6  # paper's worked example for panel (b)

H_INFERRED = 8.4e-3
"""Structure height for panel (b) [m] -- NOT stated in the paper.

Recovered by fitting an exponential to the digitised published curve. Four
independent fitting windows give sqrt(Da) = 24600-26700, i.e. h = 8.0-8.7 mm,
and a pointwise inversion of the same curve gives 8.2-8.6 mm. Adopted value
8.4 +- 0.4 mm. See the panel README for the diagnosis.
"""

H_NAIVE = 1e-3
"""Kept as a second curve purely to show how strongly the z/h picture depends
on an h the paper never states."""


def load_published():
    path = DATA / "fig05b_profile.csv"
    if not path.exists():
        return None
    d = np.loadtxt(path, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def profile(h, z_star):
    """C/C_eq at z_star for a structure of height h."""
    Da = core.Da(h, eps_s=EPS_S, d_sub=D_SUB_PROFILE)
    c_top = core.C_top_star(U_EVAP, eps_s=EPS_S, d_sub=D_SUB_PROFILE)
    return 1.0 + c_top / C_EQ * core.conc_profile(z_star, Da), Da


def main():
    fig, (ax_top, ax_prof) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # --- (a) surface supersaturation vs grain size -------------------------
    d_um = np.linspace(1.0, 100.0, 400)
    c_top_star = core.C_top_star(U_EVAP, eps_s=EPS_S, d_sub=d_um * 1e-6)
    ax_top.plot(d_um, 1.0 + c_top_star / C_EQ, color="k")
    ax_top.set_xlim(0, 100)
    ax_top.set_ylim(1.0, 1.00006)
    ax_top.set_xlabel(r"$d_{sub}$ ($\mu$m)")
    ax_top.set_ylabel(r"$C_{top}\,/\,C_{eq}$")
    ax_top.set_title("(a) supersaturation at the moving surface, eq (31)", fontsize=10)
    ax_top.ticklabel_format(axis="y", useOffset=False)

    # --- (b) concentration profile -----------------------------------------
    z_star = np.linspace(0.999, 1.0, 20000)
    pub = load_published()
    if pub is not None:
        ax_prof.plot(pub[0], pub[1], color="0.72", lw=4,
                     label="published curve (digitised)")

    prof_naive, Da_naive = profile(H_NAIVE, z_star)
    ax_prof.plot(z_star, prof_naive, "b--", lw=1.2,
                 label=rf"eq (19), $h$ = 1 mm ($\sqrt{{Da}}$ = {np.sqrt(Da_naive):.0f})")

    prof_inf, Da_inf = profile(H_INFERRED, z_star)
    ax_prof.plot(z_star, prof_inf, "k-", lw=1.5,
                 label=rf"eq (19), $h$ = {H_INFERRED * 1e3:.1f} mm inferred "
                       rf"($\sqrt{{Da}}$ = {np.sqrt(Da_inf):.0f})")

    ax_prof.set_xlim(0.999, 1.0)
    ax_prof.set_ylim(0.99999, 1.00004)
    ax_prof.set_xlabel(r"$z\,/\,h$")
    ax_prof.set_ylabel(r"$C\,/\,C_{eq}$")
    ax_prof.set_title("(b) profile, $d_{sub}$ = 10 $\\mu$m", fontsize=10)
    ax_prof.ticklabel_format(axis="y", useOffset=False)
    ax_prof.legend(frameon=False, fontsize=7.5, loc="upper left")

    delta_C = core.boundary_layer_thickness(eps_s=EPS_S, d_sub=D_SUB_PROFILE)
    ax_prof.text(
        0.99912, 1.0000185,
        rf"$\delta_C$ = {delta_C * 1e9:.0f} nm regardless of $h$." "\n"
        rf"Only $\delta_C/h$ shifts, so the paper's" "\n"
        rf"unstated $h$ sets the apparent width.",
        fontsize=7.5, color="0.35",
    )

    fig.suptitle("Prat (2026) Figure 5 -- reproduced from eqs (31) and (19)", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "fig05.png", dpi=180)
    plt.close(fig)

    # =====================================================================
    print("Figure 5 reproduction")
    print("  (a) C_top/C_eq - 1  (eq 31, independent of h):")
    for d in (1, 10, 25, 50, 100):
        cs = core.C_top_star(U_EVAP, eps_s=EPS_S, d_sub=d * 1e-6)
        print(f"      d_sub = {d:3d} um -> {cs:.3e}   C_top/C_eq = {1 + cs / C_EQ:.7f}")

    print("\n  (b) the height the paper does not state")
    print(f"      delta_C = {delta_C * 1e9:.1f} nm -- fixed by reaction/diffusion,")
    print("        the same whatever h is (delta_C = h/sqrt(Da) and Da ~ h^2)")
    for h in (H_NAIVE, H_INFERRED):
        Da = core.Da(h, eps_s=EPS_S, d_sub=D_SUB_PROFILE)
        print(f"      h = {h * 1e3:4.1f} mm -> Da = {Da:.3e}, sqrt(Da) = {np.sqrt(Da):7.0f}, "
              f"delta_C/h = {1 / np.sqrt(Da):.2e}")

    if pub is not None:
        zp, Cp = pub
        v = Cp - 1.0
        c_top = core.C_top_star(U_EVAP, eps_s=EPS_S, d_sub=D_SUB_PROFILE) / C_EQ
        coef = np.sqrt(core.Da(1.0, eps_s=EPS_S, d_sub=D_SUB_PROFILE))
        print("\n      fitting the published curve, over several windows:")
        for lo, hi in ((0.02, 0.30), (0.03, 0.35), (0.05, 0.35), (0.02, 0.35)):
            sel = (v > lo * c_top) & (v < hi * c_top)
            s = -np.polyfit(1.0 - zp[sel], np.log(v[sel] / c_top), 1)[0]
            print(f"        v/C_top in [{lo:.2f}, {hi:.2f}]: {sel.sum():3d} pts -> "
                  f"sqrt(Da) = {s:6.0f}, h = {s / coef * 1e3:.2f} mm")

        # Compare HORIZONTALLY, not vertically. The curve is near-vertical over
        # the rise, so a vertical residual is dominated by a sub-pixel shift in
        # z and says nothing useful; the meaningful question is whether the two
        # curves reach the same concentration at the same depth.
        # np.interp needs a strictly increasing abscissa; the flat part of the
        # published curve is a long run of identical zeros, so restrict to the
        # rising portion and enforce monotonicity before inverting.
        rise = v > 0.01 * c_top
        v_r, z_r = v[rise], zp[rise]
        keep = np.concatenate([[True], np.diff(v_r) > 0])
        v_r, z_r = v_r[keep], z_r[keep]

        # only up to 0.3 C_top: the extraction stops short of the tick zone
        # near z/h = 1, so the very top of the rise is not available
        levels = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
        print(f"\n      horizontal agreement of the adopted h = {H_INFERRED * 1e3:.1f} mm curve"
              " (depth at which each level is reached):")
        dev = []
        for lev in levels:
            zp_l = np.interp(lev * c_top, v_r, z_r)
            zm_l = np.interp(lev * c_top, prof_inf - 1.0, z_star)
            dev.append(abs(zp_l - zm_l) / (1.0 / np.sqrt(Da_inf)))
            print(f"        C-C_eq = {lev:.2f} C_top:  published 1-z/h = {1 - zp_l:.3e},"
                  f"  ours {1 - zm_l:.3e}   ({dev[-1]:.2f} delta_C)")
        dev = np.array(dev)
        print(f"        worst offset {dev.max():.2f} boundary-layer thicknesses")
        assert dev.max() < 0.6, "adopted h does not reproduce the published profile"

    assert np.all(np.isfinite(prof_inf)), "profile must not contain nan/inf"


if __name__ == "__main__":
    main()
