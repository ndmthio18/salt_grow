"""Figure 8 -- comparison with the sand-column experiment of Norouzi Rad et al.

(a) Cumulative precipitated salt mass, eq (35).
(b) The same data as an equivalent efflorescence thickness, eq (32).

THE MISSING AREA. Eq (35) needs A, the cross-section of the experimental
cylinder, which Prat never states (it lives in ref 18, behind a paywall). But
panels (a) and (b) show the SAME experiment in different units, related by
Ms = rho_cr(1-eps_s) h eta A. Digitising both therefore recovers the product
eta*A without any external data -- and the scatter of that product across the
13 time points is itself a check on both digitisations.

WHAT IS AND IS NOT AN INDEPENDENT TEST. Because eta*A comes from the ratio of
the two panels, panel (a) is a unit conversion of panel (b), not a second
validation. The genuinely independent test is the SLOPE of panel (b), which is
dh/dt from eq (32) and involves no fitted quantity at all.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
from params import DAY, RHO, RHO_CR, NorouziRad

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"

EPS_S = NorouziRad.EPS_S
D_SUB = NorouziRad.D_SUB
U_EVAP = NorouziRad.U_EVAP
J = RHO * U_EVAP


def load(name):
    d = np.loadtxt(DATA / name, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def main():
    t_a, Ms = load("fig08a_exp.csv")  # g
    t_b, h = load("fig08b_exp.csv")  # mm
    pub_a1 = load("fig08a_eta1.csv")
    pub_a076 = load("fig08a_eta076.csv")
    pub_b076 = load("fig08b_eta076.csv")

    # --- recover eta*A from the two panels ---------------------------------
    etaA = (Ms * 1e-3) / (RHO_CR * (1.0 - EPS_S) * (h * 1e-3))  # m^2
    etaA_mean, etaA_rsd = etaA.mean(), etaA.std() / etaA.mean()
    area = etaA_mean / NorouziRad.ETA

    # --- our predictions ----------------------------------------------------
    dh_dt = core.growth_rate(J, eps_s=EPS_S, d_sub=D_SUB)  # m/s
    dMs_dt_1 = core.dMs_dt(J, area, 1.0, eps_s=EPS_S, d_sub=D_SUB)  # kg/s
    dMs_dt_076 = core.dMs_dt(J, area, NorouziRad.ETA, eps_s=EPS_S, d_sub=D_SUB)

    # The paper anchors its lines on the first experimental point rather than
    # on the origin (efflorescence only starts once the surface saturates), so
    # we do the same in order to compare like with like.
    t0, Ms0, h0 = t_a[0], Ms[0], h[0]
    tt = np.linspace(t0, t_a[-1], 200)

    fig, (ax_m, ax_h) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # --- (a) precipitated mass ---------------------------------------------
    ax_m.plot(t_a, Ms, "ko", ms=7, label="experiment (ref 18, digitised)")
    ax_m.plot(*pub_a1, color="0.75", lw=3, label="published line, $\\eta$ = 1")
    ax_m.plot(*pub_a076, color="0.75", lw=3)
    ax_m.plot(
        tt, (Ms0 + dMs_dt_1 * (tt - t0) * 3600.0 * 1e3), "r-", lw=1.4,
        label=r"eq (35), $\eta$ = 1 (ours)",
    )
    ax_m.plot(
        tt, (Ms0 + dMs_dt_076 * (tt - t0) * 3600.0 * 1e3), "b-", lw=1.4,
        label=rf"eq (35), $\eta$ = {NorouziRad.ETA} (ours)",
    )
    ax_m.set_xlim(0, 25)
    ax_m.set_ylim(0, 0.30)
    ax_m.set_xlabel("Time (hr)")
    ax_m.set_ylabel("Precipitated salt (g)")
    ax_m.set_title("(a) precipitated mass, eq (35)", fontsize=10)
    ax_m.legend(frameon=False, fontsize=7.5, loc="upper left")

    # --- (b) equivalent thickness ------------------------------------------
    ax_h.plot(t_b, h, "ko", ms=7, label="from experiment (ref 18, digitised)")
    ax_h.plot(*pub_b076, color="0.75", lw=3, label="published line")
    ax_h.plot(
        tt, h0 + dh_dt * (tt - t0) * 3600.0 * 1e3, "b-", lw=1.4,
        label=r"eq (32) (ours)",
    )
    ax_h.set_xlim(0, 25)
    ax_h.set_ylim(0, 3)
    ax_h.set_xlabel("Time (hr)")
    ax_h.set_ylabel(r"$h_{equ}$ (mm)")
    ax_h.set_title("(b) equivalent thickness, eq (32)", fontsize=10)
    ax_h.legend(frameon=False, fontsize=7.5, loc="upper left")

    fig.suptitle(
        "Prat (2026) Figure 8 -- sand column of Norouzi Rad et al. (ref 18)", fontsize=11
    )
    fig.tight_layout()
    fig.savefig(HERE / "fig08.png", dpi=180)
    plt.close(fig)

    # =====================================================================
    # Report
    # =====================================================================
    print("Figure 8 reproduction")
    print(f"  U_evap = {U_EVAP * DAY * 1e3:.1f} mm/day, d_sub = {D_SUB * 1e3:.2f} mm, "
          f"eps_s = {EPS_S}")
    print(f"  Da* (ours) = {core.Da_star(eps_s=EPS_S, d_sub=D_SUB):.1f}   "
          f"(paper text quotes {NorouziRad.DA_STAR_PAPER:.0f})")

    print("\n  Cross-check of the two digitisations")
    print(f"    max |t_a - t_b| over the 13 markers = {np.abs(t_a - t_b).max():.3f} hr")
    print(f"    recovered eta*A = {etaA_mean * 1e4:.4f} cm^2, "
          f"relative scatter {etaA_rsd * 100:.2f}% (accept < 5%)")
    print(f"    with eta = {NorouziRad.ETA} -> A = {area * 1e4:.3f} cm^2 "
          f"-> equivalent diameter {2 * np.sqrt(area / np.pi) * 100:.2f} cm")
    print("      (a ~1 cm column is plausible for an X-ray microtomography sample,")
    print("       so the recovered area is physically sensible -- ref 18 is paywalled")
    print("       and could not be consulted to confirm it directly.)")

    print("\n  Growth rate -- the genuinely independent test")
    pub_slope = np.polyfit(*pub_b076, 1)[0] * 24.0  # mm/day
    ours = dh_dt * DAY * 1e3
    print(f"    published line  {pub_slope:.4f} mm/day")
    print(f"    ours, eq (32)   {ours:.4f} mm/day")
    print(f"    ratio           {ours / pub_slope:.4f}   <-- agreement to "
          f"{abs(ours / pub_slope - 1) * 100:.2f}%  (model vs the PAPER'S LINE)")

    # --- the experimental points themselves --------------------------------
    # The 0.08% above certifies "our implementation = Prat's implementation".
    # Fitting the 13 digitised EXPERIMENTAL markers directly is the stronger
    # statement: model vs nature (through one digitisation).
    exp_fit = np.polyfit(t_b, h, 1)
    exp_slope = exp_fit[0] * 24.0                                  # mm/day
    exp_resid = h - np.polyval(exp_fit, t_b)
    print("\n  Growth rate vs the 13 EXPERIMENTAL points (self-fit)")
    print(f"    experimental slope  {exp_slope:.4f} mm/day "
          f"(fit residual rms {exp_resid.std():.3f} mm)")
    print(f"    ours, eq (32)       {ours:.4f} mm/day")
    print(f"    ratio               {ours / exp_slope:.4f}   <-- agreement to "
          f"{abs(ours / exp_slope - 1) * 100:.2f}%")

    print("\n  Mass rates, eq (35)")
    pa1 = np.polyfit(*pub_a1, 1)[0] * 1e3  # mg/hr
    pa076 = np.polyfit(*pub_a076, 1)[0] * 1e3
    print(f"    eta = 1     published {pa1:7.4f} mg/hr   ours {dMs_dt_1 * 3600 * 1e6:7.4f} mg/hr")
    print(f"    eta = 0.76  published {pa076:7.4f} mg/hr   ours {dMs_dt_076 * 3600 * 1e6:7.4f} mg/hr")
    print(f"    published slope ratio {pa076 / pa1:.4f}  (should equal eta = {NorouziRad.ETA})")

    assert etaA_rsd < 0.05, "eta*A is not constant -- a digitisation is biased"
    assert abs(ours / pub_slope - 1) < 0.05, "growth rate disagrees with the published line"


if __name__ == "__main__":
    main()
