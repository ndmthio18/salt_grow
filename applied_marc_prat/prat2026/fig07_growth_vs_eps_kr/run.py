"""Figure 7 -- sensitivity of the growth rate to porosity and reaction coefficient.

(a) vs salt structure porosity eps_s. The dependence is weak, and it is mostly
    the explicit 1/(1-eps_s) prefactor of eq (32) rather than Da*. The practical
    reading: if precipitation clogs the structure and lowers eps_s, growth slows
    down, but only mildly.

(b) vs precipitation coefficient k_r over three decades. Here the dependence is
    strong -- an order of magnitude in dh/dt. This panel is the paper's argument
    for taking k_r = 2.3e-3 m/s (the microfluidic measurement of ref 23) rather
    than the ~1e-6 m/s used in an earlier numerical study, which the paper
    attributes to that study being unable to resolve the thin boundary layer.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import KR_T21

HERE = pathlib.Path(__file__).resolve().parent

D_SUB = 20e-6  # paper text, for both panels
EPS_S = 0.43  # paper text, for panel (b)


def main():
    fig, (ax_eps, ax_kr) = plt.subplots(1, 2, figsize=(10.5, 4.2))

    # --- (a) vs structure porosity ----------------------------------------
    eps = np.linspace(0.1, 0.5, 400)
    ax_eps.plot(eps, core.growth_ratio(eps_s=eps, d_sub=D_SUB), color="k")
    ax_eps.set_xlim(0.1, 0.5)
    ax_eps.set_xlabel(r"$\varepsilon_s$")
    ax_eps.set_ylabel(r"$(dh/dt)\,/\,U_{evap}$")
    ax_eps.set_title(rf"(a) vs porosity ($d_{{sub}}$ = 20 $\mu$m)", fontsize=10)

    # --- (b) vs reaction coefficient, log x --------------------------------
    kr = np.logspace(np.log10(2.3e-6), np.log10(2.3e-3), 400)
    ax_kr.semilogx(kr, core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=kr), color="k")
    ax_kr.set_xlim(1e-6, 1e-2)
    ax_kr.set_ylim(0, 0.27)
    ax_kr.set_xlabel(r"$k_r$ (m/s)")
    ax_kr.set_ylabel(r"$(dh/dt)\,/\,U_{evap}$")
    ax_kr.set_title(
        rf"(b) vs reaction coefficient ($\varepsilon_s$ = 0.43, $d_{{sub}}$ = 20 $\mu$m)",
        fontsize=10,
    )
    ax_kr.plot(
        [KR_T21],
        [core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=KR_T21)],
        "o",
        mfc="none",
        color="crimson",
        ms=7,
    )
    ax_kr.annotate(
        "measured $k_r$\n(ref 23)",
        xy=(KR_T21, core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=KR_T21)),
        xytext=(6e-5, 0.215),
        fontsize=8,
        color="crimson",
        arrowprops=dict(arrowstyle="->", color="crimson", lw=0.8),
    )

    fig.suptitle("Prat (2026) Figure 7 -- parameter sensitivity of eq (32)", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "fig07.png", dpi=180)
    plt.close(fig)

    print("Figure 7 reproduction")
    print("  (a) vs porosity, d_sub = 20 um:")
    for e in (0.1, 0.2, 0.3, 0.43, 0.5):
        print(f"      eps_s = {e:.2f} -> (dh/dt)/Uevap = {core.growth_ratio(eps_s=e, d_sub=D_SUB):.4f}")
    print("  (b) vs reaction coefficient, eps_s = 0.43:")
    for k in (2.3e-6, 2.3e-5, 2.3e-4, 2.3e-3):
        print(f"      kr = {k:.1e} -> (dh/dt)/Uevap = {core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=k):.4f}")
    lo = core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=2.3e-6)
    hi = core.growth_ratio(eps_s=EPS_S, d_sub=D_SUB, kr=2.3e-3)
    print(f"\n  k_r spans 3 decades -> growth spans a factor {hi / lo:.1f}")
    print(
        "  porosity spans 0.1-0.5 -> growth spans a factor "
        f"{core.growth_ratio(eps_s=0.5, d_sub=D_SUB) / core.growth_ratio(eps_s=0.1, d_sub=D_SUB):.2f}"
        "  (the 'weak dependence' the paper notes)"
    )


if __name__ == "__main__":
    main()
