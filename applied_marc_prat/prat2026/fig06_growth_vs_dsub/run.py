"""Figure 6 -- growth rate and induced Darcy velocity vs substrate grain size.

Main panel: (dh/dt)/U_evap from eq (32). The ratio lands in 0.1-0.3, i.e. the
salt front advances at a sizeable fraction of the evaporation velocity, and
grows with d_sub because Da* does.

Inset: U/U_evap from eq (33), always above 1. The liquid drawn up through the
structure must not only replace evaporated water but also fill the pore space
of the newly created structure and carry the salt that becomes crystal.

This figure carries the calibration of params.RHO: the paper states dh/dt =
0.25 cm/day at d_sub = 25 um and U_evap = 1 cm/day, which we reproduce as
0.254, and the curve passes through ~0.219 at 10 um matching the published
figure. Those two points are the only evidence for rho = 1200 kg/m^3.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core

HERE = pathlib.Path(__file__).resolve().parent

EPS_S = 0.43


def main():
    fig, ax = plt.subplots(figsize=(6.4, 4.8))

    d_um = np.linspace(1.0, 100.0, 500)
    ratio = core.growth_ratio(eps_s=EPS_S, d_sub=d_um * 1e-6)
    ax.plot(d_um, ratio, color="k")
    ax.set_xlim(0, 100)
    ax.set_ylim(0.1, 0.35)
    ax.set_xlabel(r"$d_{sub}$ ($\mu$m)")
    ax.set_ylabel(r"$(dh/dt)\,/\,U_{evap}$")
    ax.set_title("Prat (2026) Figure 6 -- growth law, eq (32)", fontsize=11)

    # The paper's own quoted anchor -- also the sole support for params.RHO.
    anchor = core.growth_ratio(eps_s=EPS_S, d_sub=25e-6)
    ax.plot([25.0], [anchor], "o", mfc="none", color="crimson", ms=8)
    ax.annotate(
        f"paper text: 0.25 cm/day at $d_{{sub}}$ = 25 $\\mu$m\nours: {anchor:.3f}"
        "\n(the calibration anchor for $\\rho$)",
        xy=(25.0, anchor),
        xytext=(31.0, 0.315),
        fontsize=8,
        color="crimson",
        arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8),
    )

    # --- inset: induced Darcy velocity ------------------------------------
    ax_u = ax.inset_axes([0.52, 0.17, 0.42, 0.34])
    ax_u.plot(d_um, core.U_over_Uevap(eps_s=EPS_S, d_sub=d_um * 1e-6), color="k")
    ax_u.set_xlim(0, 100)
    ax_u.set_xlabel(r"$d_{sub}$ ($\mu$m)", fontsize=7)
    ax_u.set_ylabel(r"$U\,/\,U_{evap}$", fontsize=8)
    ax_u.tick_params(labelsize=6)

    fig.tight_layout()
    fig.savefig(HERE / "fig06.png", dpi=180)
    plt.close(fig)

    print("Figure 6 reproduction")
    for d in (1, 10, 25, 50, 100):
        r = core.growth_ratio(eps_s=EPS_S, d_sub=d * 1e-6)
        u = core.U_over_Uevap(eps_s=EPS_S, d_sub=d * 1e-6)
        print(f"    d_sub = {d:3d} um -> (dh/dt)/Uevap = {r:.4f},  U/Uevap = {u:.4f}")
    print("\n  RHO calibration check (the only support for params.RHO = 1200):")
    print(f"    25 um -> {core.growth_ratio(eps_s=EPS_S, d_sub=25e-6):.4f}  (paper text: 0.25)")
    print(f"    10 um -> {core.growth_ratio(eps_s=EPS_S, d_sub=10e-6):.4f}  (paper figure: ~0.21-0.22)")
    assert 0.24 < core.growth_ratio(eps_s=EPS_S, d_sub=25e-6) < 0.26


if __name__ == "__main__":
    main()
