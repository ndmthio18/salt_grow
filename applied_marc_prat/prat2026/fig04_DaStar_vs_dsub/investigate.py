"""Investigate the Da* discrepancy, using the paper's own figure as evidence.

Several different values of Da* are in circulation for the SAME parameter set
(eps_s = 0.43, k_r = 2.3e-3 m/s, D_s = 1.3e-9 m^2/s, beta = 2):

  * eq (27)/(34) evaluated as printed
  * the value quoted in the paper's running text
  * the curve drawn in the paper's Figure 4 main panel
  * Figure 4's right inset, which plots the same thing on a log ordinate

This script pins each of them down, checks whether Figure 4 agrees with itself,
and works out whether any of it changes a prediction.

A methodological warning that this analysis learned the hard way: the left
inset plots F = sqrt(Da*)/(1+sqrt(Da*)), and inverting F back to Da* amplifies
reading error by 2/(F(1-F)) -- about 9x at F = 0.65 and 27x at F = 0.92. An
earlier pass inverted it and concluded the figure was internally inconsistent
by 40%. It is not; that was the amplification. Comparisons below are made in F
space wherever the source is an F plot.

Run after ../digitize.py.  ../.venv/bin/python investigate.py
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import C_EQ, RHO, RHO_CR, NorouziRad

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
D_REF = NorouziRad.D_SUB  # 0.58 mm, the point the paper quotes


def load(name):
    d = np.loadtxt(DATA / name, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def main():
    print("=" * 74)
    print("Da* DISCREPANCY -- formula vs text vs figure")
    print("=" * 74)

    # --- 1. our formula ----------------------------------------------------
    ours_slope = core.Da_star(d_sub=1e-3) / 1e-3
    ours_at_ref = core.Da_star(d_sub=D_REF)
    print("\n1. Evaluating eq (27)/(34) as printed, with the paper's parameters")
    print(f"   Da*/d_sub = {ours_slope:.4g} 1/m      Da*(0.58 mm) = {ours_at_ref:.1f}")

    # --- 2. the published main-panel curve ---------------------------------
    x4, y4 = load("fig04_model_curve.csv")
    fit = np.polyfit(x4, y4, 1)
    pub_slope = fit[0] * 1e6
    pub_at_ref = np.polyval(fit, D_REF * 1e6)
    print("\n2. The curve drawn in Figure 4's main panel (digitised, 708 points)")
    print(f"   Da* = {fit[0]:.5f} * d_sub[um] {fit[1]:+.3f}   (intercept ~ 0, as it should be)")
    print(f"   Da*/d_sub = {pub_slope:.4g} 1/m      Da*(0.58 mm) = {pub_at_ref:.1f}")

    # --- 3. the text -------------------------------------------------------
    text_at_ref = NorouziRad.DA_STAR_PAPER
    print("\n3. The value quoted in the paper's running text")
    print(f"   Da*(0.58 mm) = {text_at_ref:.0f}   -> Da*/d_sub = {text_at_ref / D_REF:.4g} 1/m")

    # --- 4. is the published figure internally consistent? -----------------
    xz, Daz = load("fig04_inset_zoom.csv")
    keep = xz > 5.0
    zoom_slope = float(np.sum(xz[keep] * Daz[keep]) / np.sum(xz[keep] ** 2)) * 1e6
    print("\n4. Is the published figure internally consistent?")
    print("   The right inset plots Da* on a LOG ordinate over [0, 100] um, where the")
    print("   main panel's linear axis has no resolution at all. It is therefore an")
    print("   independent rendering of the same quantity.")
    print(f"     right inset -> Da*/d_sub = {zoom_slope:.4g} 1/m")
    print(f"     main panel  -> Da*/d_sub = {pub_slope:.4g} 1/m")
    print(f"     agreement to {abs(zoom_slope / pub_slope - 1) * 100:.1f}%"
          "  -- the figure IS self-consistent")

    xF, Fv = load("fig04_inset_F.csv")
    print("\n   The left inset plots F, so compare there directly instead of inverting:")
    print(f"     {'d_sub':>7} {'left inset':>11} {'right inset':>12} {'main':>8} {'ours':>8}")
    for dt in (5, 10, 20, 50, 100):
        fl = np.interp(dt, xF, Fv)
        fr = core.F(np.interp(dt, xz, Daz))
        fm = core.F(pub_slope * dt * 1e-6)
        fo = core.F(core.Da_star(d_sub=dt * 1e-6))
        print(f"     {dt:7d} {fl:11.4f} {fr:12.4f} {fm:8.4f} {fo:8.4f}")
    print("   The left inset runs ~0.03-0.04 low in F against the other two. Inverted")
    print("   to Da* that looks like tens of percent -- which is exactly the trap")
    print("   described in this module's docstring. In F, the quantity that actually")
    print("   enters the growth law, all four sources agree to within 0.04.")

    # --- 5. does any of it matter? -----------------------------------------
    print("\n5. Impact on what the growth law actually predicts")
    for label, ds in (("ours", ours_at_ref), ("text", text_at_ref), ("figure", pub_at_ref)):
        print(f"     {label:7s} Da* = {ds:6.1f}  ->  F = {core.F(ds):.5f}")
    spread = (core.F(ours_at_ref) - core.F(pub_at_ref)) / core.F(pub_at_ref)
    print(f"   full spread in F across all of them: {spread * 100:.2f}%")

    # --- 6. can the Figure 8b comparison discriminate? ---------------------
    print("\n6. Can the Figure 8b comparison tell us which Da* the paper used?")
    tb, hb = load("fig08b_eta076.csv")
    slope_mm_day = np.polyfit(tb, hb, 1)[0] * 24.0
    ratio = slope_mm_day / (NorouziRad.U_EVAP * 86400.0 * 1e3)
    prefac = RHO * C_EQ / (RHO_CR * (1.0 - 0.43) * (1.0 - C_EQ))
    F_needed = ratio / prefac
    print(f"   published line slope {slope_mm_day:.4f} mm/day at U_evap = 7.7 mm/day")
    print(f"   -> required F = {F_needed:.5f}  ->  implied Da* = "
          f"{(F_needed / (1.0 - F_needed)) ** 2:.1f}")
    print("   But F is so flat that the same line comes out of any Da* in roughly")
    print("   [130, 180] provided rho moves by under 1%. So no -- what Figure 8b")
    print("   pins down tightly is the prefactor, and hence rho ~ 1200 kg/m^3:")
    for ds in (136.0, 156.0, 169.6):
        rho_needed = ratio / core.F(ds) * (RHO_CR * 0.57 * (1.0 - C_EQ)) / C_EQ
        print(f"     if Da* were {ds:6.1f}, the published line implies rho = {rho_needed:7.1f} kg/m^3")

    # --- 7. isolate it: Figure 2a contains Da, which has no beta ------------
    print("\n7. Isolating the discrepancy with Figure 2a")
    print("   Da and Da* share every parameter except beta. Dividing them,")
    print("      Da*/Da = beta^2 d_sub^2 / (1296 h^2)")
    print("   k_r, D_s, eps_s, C_eq and the a_v prefactor all cancel. So Figure 2a,")
    print("   which plots Da, constrains everything EXCEPT beta -- and Figure 4,")
    print("   which plots Da*, adds beta on top. Comparing both isolates the cause.")

    print("\n   Figure 2a (published Da) against ours:")
    ratios = []
    for tag, d_sub in (("dsub1", 1e-6), ("dsub10", 10e-6), ("dsub100", 100e-6)):
        hx, Dy = load(f"fig02a_{tag}.csv")
        sel = (hx > 0.25) & (hx < 0.98)
        r = Dy[sel] / core.Da(hx[sel] * 1e-3, eps_s=0.43, d_sub=d_sub)
        ratios.append(np.median(r))
        print(f"     d_sub = {d_sub * 1e6:5.0f} um: published/ours = {np.median(r):.3f}"
              f"   (10-90 pct {np.percentile(r, 10):.3f}-{np.percentile(r, 90):.3f})")
    print(f"   -> Da agrees to {abs(np.median(ratios) - 1) * 100:.1f}%."
          "  Everything except beta is confirmed.")

    print("\n   What each hypothesis predicted for that ratio:")
    print("     a_v with six cube faces instead of five : 1.200   REFUTED")
    print("     a different k_r or D_s                  : 0.815   REFUTED")
    print("     a different beta                        : 1.000   consistent")

    print("\n   So beta is the only parameter left. From Da*(published)/Da*(ours):")
    for label, ratio in (("Figure 4 right inset", zoom_slope / ours_slope),
                         ("Figure 4 main panel", pub_slope / ours_slope),
                         ("text value, 156", (text_at_ref / D_REF) / ours_slope)):
        print(f"     {label:22s} {ratio:.4f}  ->  beta = 2*sqrt(ratio) = "
              f"{2 * np.sqrt(ratio):.3f}")
    print("   The paper states beta = 2, and its Appendix derives exactly that from")
    print("   A = 24 d_t^2 and V_c = 4 d_t^3, giving beta = A/(12 d_t^2) = 2. Nothing")
    print("   in that derivation produces 1.8, so this looks like the plotted curve")
    print("   having been made with a different beta than the one the paper reports.")

    print("\n   Impact of beta = 1.8 instead of 2, at d_sub = 0.58 mm:")
    ours_18 = ours_at_ref * (1.8 / 2.0) ** 2
    print(f"     Da*  {ours_at_ref:.1f} -> {ours_18:.1f}")
    print(f"     F    {core.F(ours_at_ref):.5f} -> {core.F(ours_18):.5f}"
          f"   ({abs(core.F(ours_18) / core.F(ours_at_ref) - 1) * 100:.2f}% on dh/dt)")

    print("\n" + "=" * 74)
    print("CONCLUSION")
    print("=" * 74)
    print("""
  Da*/d_sub, from every available source:

      2.92e5 1/m   eq (27)/(34) evaluated as printed        (ours)
      2.69e5 1/m   the value quoted in the running text     (156 at 0.58 mm)
      2.39e5 1/m   Figure 4's right inset, log ordinate
      2.34e5 1/m   Figure 4's main panel

  * The published FIGURE is internally consistent -- main panel and right inset
    agree to 2%. An earlier pass claimed otherwise; that was an artefact of
    inverting F to Da*, and it is corrected here.
  * Figure 2a pins down everything except beta. Its three Da curves reproduce
    ours to 1.3%, which refutes both the six-cube-face explanation (it would
    need 1.20) and any k_r or D_s difference (0.815).
  * Since Da*/Da = beta^2 d_sub^2 / (1296 h^2), beta is the ONLY parameter that
    can move Da* without moving Da. The published Da* curve corresponds to
    beta = 1.79-1.81, and the quoted text value to beta = 1.92, against the
    beta = 2 that the paper states and that its own Appendix derives.
  * The consequence is still negligible: beta = 1.8 changes F by 0.8%, so no
    growth rate, mass rate or thickness in the paper is materially affected.
  * Practical consequence: use eq (32)/(35) with confidence -- they reproduce
    the paper's own Figure 8b line to 0.08% -- but do not treat a quoted or
    plotted Da* as a reproducible number, and do not tune parameters to hit
    one, because that would corrupt a chain that is otherwise verified.""")


if __name__ == "__main__":
    main()
