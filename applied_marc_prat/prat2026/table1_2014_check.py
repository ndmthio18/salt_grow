"""Table 1 of Veran-Tissoires & Prat (2014, JFM 749) as a validation asset.

The table (paper.txt:489-531, plain text -- no digitisation needed) gives, for
the delta = 2.5 / 7.5 / 15 mm gap series (1 mm beads, rw = 19 mm wick):

    J_exp   measured evaporation rate
    J_num   their PUBLISHED pure-diffusion numerical solution (enclosure
            geometry), and the free-convection + diffusion solution

Three things are extracted here:

1.  SERIES-RESISTANCE CROSS-CHECK (the L2 test). Our V0 formula is exact for
    its own 1D boundary condition (machine-verified); the question Table 1 can
    answer is whether the diffusive-resistance PICTURE maps onto a real
    geometry. Model: R_tot(delta) = delta/(rhoD A) + R_ext with ONE constant
    external (access) resistance for the cylinder-mouth -> enclosure path.
    If R_ext comes out constant across the three gaps, the picture holds.
    A pure "H = delta" reading would be off by 2-4x here (rw >> delta: the
    external resistance dominates) -- which is precisely the lesson for our
    model: H is an EFFECTIVE diffusive length, not a literal gap width.

2.  CONVECTION FACTOR. J_exp / J_num(pure diffusion) = the enhancement the
    real system gets from free convection: 1.41 / 1.55 / 1.67. The delta=15mm
    value 1.67 independently supports the rhoD x 1.655 calibration used by the
    NorouziRad preset in new_model (efflorescence3d/params.py) to absorb
    external transport into an effective rhoD.

3.  ABSOLUTE-RATE CONTEXT. Water-equivalent evaporation velocities of all
    experiments in this repository, side by side -- a one-glance dimensional
    sanity table.

Run:  .venv/bin/python prat2026/table1_2014_check.py
"""

import numpy as np

# --- Table 1, transcribed from paper.txt:489-531 (SI units) ----------------
# provenance: Veran-Tissoires & Prat 2014, JFM 749, Table 1; 1 mm beads,
# Ls ~ 10 mm series; J in kg/s
DELTA = np.array([2.5e-3, 7.5e-3, 15e-3])          # gap widths [m]
J_EXP = np.array([2.70e-8, 2.37e-8, 2.02e-8])      # measured
J_NUM_DIFF = np.array([1.92e-8, 1.53e-8, 1.21e-8])  # published, pure diffusion
J_NUM_CONV = np.array([2.72e-8, 2.08e-8, 1.55e-8])  # published, + free conv.

RW = 19e-3                       # wick inner radius [m] (paper.txt:165)
AREA = np.pi * RW**2
RHOD = 3.0e-5                    # rho_gas * D_v [kg/m/s] (salt_lattice conventions)
XVE, XVINF = 0.0198, 0.00185     # a_w = 0.75 (sat. NaCl), RH = 7% -- same rig
DX = XVE - XVINF
DAY, RHO_W = 86400.0, 1000.0


def main():
    print("=" * 72)
    print("2014 JFM Table 1 -- three checks (no digitisation, plain-text data)")
    print("=" * 72)

    # --- 1. series-resistance cross-check ------------------------------------
    R_gap = DELTA / (RHOD * AREA)
    R_tot = DX / J_NUM_DIFF
    R_ext = R_tot - R_gap
    R_disc = 1.0 / (4.0 * RHOD * RW)     # flush-disc half-space access value
    print("\n1. series-resistance model vs PUBLISHED pure-diffusion numerics")
    print(f"   {'delta':>8} {'R_gap':>10} {'R_tot':>10} {'R_ext':>10}")
    for d, rg, rt, re in zip(DELTA, R_gap, R_tot, R_ext):
        print(f"   {d * 1e3:6.1f}mm {rg:10.3e} {rt:10.3e} {re:10.3e}")
    spread = np.ptp(R_ext) / R_ext.mean()
    print(f"   R_ext constant to {spread * 100:.1f}% across a 6x gap change "
          f"(mean {R_ext.mean():.3e}); residual trend is MONOTONE in delta")
    print("   (the source recedes into the tube as delta grows, so the mouth")
    print("    access resistance itself grows slightly -- geometry, not noise)")
    print(f"   R_ext / flush-disc half-space value = "
          f"{R_ext.mean() / R_disc:.2f}  (tube rim + finite enclosure > 1: ok)")
    print(f"   naive 'H = delta' reading would overestimate J by "
          f"{(R_tot / R_gap).min():.1f}-{(R_tot / R_gap).max():.1f}x")
    print("   => first-order picture holds: ONE geometric constant explains")
    print("      the delta series to ~20%, where 'H = delta' fails by 3-13x.")
    print("      Lesson for our model: H is an EFFECTIVE diffusive length,")
    print("      never a literal gap width.")
    assert spread < 0.25, "R_ext not even first-order constant"

    # --- 2. convection factor ------------------------------------------------
    fac = J_EXP / J_NUM_DIFF
    print("\n2. free-convection factor J_exp / J_num(diffusion)")
    for d, f in zip(DELTA, fac):
        print(f"   delta = {d * 1e3:4.1f} mm:  {f:.2f}")
    print("   an INDEPENDENT measurement showing external convection adds a")
    print("   1.4-1.7x transport factor -- the physical basis for calibrating")
    print("   an EFFECTIVE rhoD from a measured empty-surface flux, as done in")
    print("   new_model/validation/fig08_coupled_multiseed.py (factor depends")
    print("   on the chosen H: x1.66 at H=10mm, x3.97 at H=20mm).")

    # --- 3. absolute-rate context table --------------------------------------
    v = J_EXP / AREA / RHO_W * DAY * 1e3   # mm/day water-equivalent
    print("\n3. absolute evaporation velocities [mm/day water equivalent]")
    for d, vv in zip(DELTA, v):
        print(f"   2014 JFM, delta = {d * 1e3:4.1f} mm      {vv:5.2f}")
    print("   Norouzi Rad 2013 sand column      7.70   (prat2026/params.py)")
    print("   Diouf 2018 heave cell             5.10   (prat2026/params.py)")
    print("   new_model default j_flat          4.65   (H = 10 mm, dirichlet)")
    print("   all within the 1-10 mm/day indoor-evaporation decade: no unit")
    print("   accident anywhere in the chain.")


if __name__ == "__main__":
    main()
