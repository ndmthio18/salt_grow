"""Figure 2 -- how Da and Pe vary with salt structure thickness.

(a) Da(h) for three substrate grain sizes, log y. Establishes that Da >> 1 for
    any structure thicker than a few microns, which licenses the Da >> 1 step
    that turns eq (28) into eq (29).
(b) Pe(h) at U_evap = 1 cm/day, linear. Establishes Pe = O(1), so Pe << Da and
    the convective term of eq (13) can be dropped to reach eq (18).

Together these two panels are the entire justification for the paper's two
strongest approximations, which is why they come first.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import DAY

HERE = pathlib.Path(__file__).resolve().parent

D_SUB_VALUES = (1e-6, 10e-6, 100e-6)  # paper: 1, 10, 100 um
U_EVAP = 0.01 / DAY  # paper: 1 cm/day


def main():
    fig, (ax_da, ax_pe) = plt.subplots(1, 2, figsize=(10.5, 4.2))

    # --- (a) Damkoehler number, log y, h in [0, 1] mm ----------------------
    h_mm = np.linspace(1e-3, 1.0, 600)
    h = h_mm * 1e-3
    for d_sub in D_SUB_VALUES:
        ax_da.semilogy(
            h_mm, core.Da(h, d_sub=d_sub), label=rf"$d_{{sub}} = {d_sub * 1e6:.0f}\ \mu$m"
        )
    # Da = 10 (the eq-29 threshold) sits below this window entirely -- it is
    # reached within a few microns, which is the point. Reported in stdout.
    ax_da.set_xlim(0, 1)
    ax_da.set_ylim(1e2, 1e8)
    ax_da.set_xlabel("$h$ (mm)")
    ax_da.set_ylabel("$Da$")
    ax_da.set_title("(a) reaction vs diffusion", fontsize=10)
    ax_da.legend(frameon=False, fontsize=9)

    # --- (b) Peclet number, linear, h in [0, 10] mm ------------------------
    h_mm_b = np.linspace(0.0, 10.0, 400)
    ax_pe.plot(h_mm_b, core.Pe(U_EVAP, h_mm_b * 1e-3), color="k")
    ax_pe.set_xlim(0, 10)
    ax_pe.set_ylim(0, 4)
    ax_pe.set_xlabel("$h$ (mm)")
    ax_pe.set_ylabel("$Pe$")
    ax_pe.set_title(r"(b) convection vs diffusion, $U_{evap}$ = 1 cm/day", fontsize=10)

    fig.suptitle(
        "Prat (2026) Figure 2 -- reproduced from eq (13) definitions", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(HERE / "fig02.png", dpi=180)
    plt.close(fig)

    # --- numbers worth printing, so the run self-reports -------------------
    print("Figure 2 reproduction")
    print("  Da at h = 1 mm:")
    for d_sub in D_SUB_VALUES:
        print(f"    d_sub = {d_sub * 1e6:5.0f} um -> Da = {core.Da(1e-3, d_sub=d_sub):.3e}")
    print(f"  Pe at h = 10 mm, U_evap = 1 cm/day -> {core.Pe(U_EVAP, 10e-3):.3f}")
    print(
        "  Da exceeds 10 above h = "
        f"{_h_where_Da_equals(10.0, 100e-6) * 1e6:.2f} um (coarsest substrate)"
    )


def _h_where_Da_equals(target, d_sub):
    """Da grows as h^2, so invert analytically rather than searching."""
    return np.sqrt(target / core.Da(1.0, d_sub=d_sub))


if __name__ == "__main__":
    main()
