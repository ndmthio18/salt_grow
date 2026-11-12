"""Figure 9 -- comparison with the surface-heave experiment of Diouf (ref 41).

(a) Measured evaporated mass. The paper reads U_evap = 5.1 mm/day off this
    curve by a linear fit; since the vessel diameter (7 mm) IS given, we can
    redo that step ourselves and check it -- an independent test of both the
    digitisation and params.RHO.
(b) Subflorescence thickness against time, with eq (32) for two porosities.
    Unlike Figure 8, these model lines pass through the origin.

At 45 C both k_r and D_s differ from their room-temperature values: k_r is
rescaled by Arrhenius (49.1 kJ/mol) to 10.5e-3 m/s and D_s is taken as
2.5e-9 m^2/s from ref 42.
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import core
import params
from params import DAY, RHO, DioufHeave

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"

KR45 = core.kr_at_temperature(params.T45_K)
AREA = np.pi * (DioufHeave.CYLINDER_DIAMETER / 2.0) ** 2


def load(name):
    d = np.loadtxt(DATA / name, delimiter=",", skiprows=6)
    return d[:, 0], d[:, 1]


def main():
    t_ev, m_ev = load("fig09a_exp.csv")
    t_ex, h_ex = load("fig09b_exp.csv")
    pub43 = load("fig09b_eps043.csv")
    pub34 = load("fig09b_eps034.csv")

    # --- re-derive U_evap from panel (a) ------------------------------------
    # The curve bends after ~28 hr (evaporation slows), and the paper fits the
    # linear part, so restrict the fit to match what it did.
    lin = t_ev < 28.0
    slope_g_per_hr = np.polyfit(t_ev[lin], m_ev[lin], 1)[0]
    dm_dt = slope_g_per_hr * 1e-3 / 3600.0  # kg/s
    U_from_curve = dm_dt / (RHO * AREA)  # m/s

    # --- our model lines ----------------------------------------------------
    j = RHO * DioufHeave.U_EVAP
    rate43 = core.growth_rate(
        j, eps_s=DioufHeave.EPS_S_BASE, d_sub=DioufHeave.D_SUB, kr=KR45, Ds=DioufHeave.DS
    )
    rate34 = core.growth_rate(
        j, eps_s=DioufHeave.EPS_S_FITTED, d_sub=DioufHeave.D_SUB, kr=KR45, Ds=DioufHeave.DS
    )
    tt = np.linspace(0.0, 36.0, 100)

    fig, (ax_m, ax_h) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # --- (a) evaporated mass -------------------------------------------------
    ax_m.plot(t_ev, m_ev, "k-", lw=1.2, label="experiment (ref 41, digitised)")
    ax_m.plot(
        t_ev[lin], np.polyval(np.polyfit(t_ev[lin], m_ev[lin], 1), t_ev[lin]),
        "r--", lw=1.2,
        label=rf"linear fit $\to U_{{evap}}$ = {U_from_curve * DAY * 1e3:.2f} mm/day",
    )
    ax_m.set_xlim(0, 40)
    ax_m.set_ylim(0, 0.4)
    ax_m.set_xlabel("Time (hr)")
    ax_m.set_ylabel("Evaporated mass (g)")
    ax_m.set_title("(a) evaporation, and the $U_{evap}$ it implies", fontsize=10)
    ax_m.legend(frameon=False, fontsize=8, loc="upper left")

    # --- (b) surface heave ---------------------------------------------------
    ax_h.plot(*pub43, color="0.75", lw=3, label="published lines")
    ax_h.plot(*pub34, color="0.75", lw=3)
    ax_h.plot(t_ex, h_ex, "ko", ms=8, label="experiment (ref 41, digitised)")
    ax_h.plot(tt, rate43 * tt * 3600.0 * 1e3, "r-", lw=1.4,
              label=rf"eq (32), $\varepsilon_s$ = {DioufHeave.EPS_S_BASE} (ours)")
    ax_h.plot(tt, rate34 * tt * 3600.0 * 1e3, "b-", lw=1.4,
              label=rf"eq (32), $\varepsilon_s$ = {DioufHeave.EPS_S_FITTED} (ours)")
    ax_h.set_xlim(0, 40)
    ax_h.set_ylim(0, 2.2)
    ax_h.set_xlabel("Time (hr)")
    ax_h.set_ylabel(r"$h$ or $h_{top}$ (mm)")
    ax_h.set_title("(b) subflorescence growth, eq (32)", fontsize=10)
    ax_h.legend(frameon=False, fontsize=8, loc="upper left")

    fig.suptitle("Prat (2026) Figure 9 -- surface heave, Diouf thesis (ref 41)", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "fig09.png", dpi=180)
    plt.close(fig)

    # =====================================================================
    print("Figure 9 reproduction")
    print(f"  T = 45 C:  k_r = {KR45:.3e} m/s (paper: 10.5e-3), D_s = {DioufHeave.DS:.1e} m^2/s")
    print(f"  Da* (ours):  eps_s=0.43 -> "
          f"{core.Da_star(eps_s=0.43, d_sub=DioufHeave.D_SUB, kr=KR45, Ds=DioufHeave.DS):.1f} "
          f"(paper text 57);  eps_s=0.34 -> "
          f"{core.Da_star(eps_s=0.34, d_sub=DioufHeave.D_SUB, kr=KR45, Ds=DioufHeave.DS):.1f} "
          f"(paper text 92)")

    print("\n  (a) U_evap re-derived from the evaporated mass curve")
    print(f"      vessel diameter {DioufHeave.CYLINDER_DIAMETER * 1e3:.0f} mm "
          f"-> A = {AREA * 1e6:.2f} mm^2")
    print(f"      linear fit over t < 28 hr: dm/dt = {slope_g_per_hr * 1e3:.4f} mg/hr")
    print(f"      U_evap = (dm/dt)/(rho A) = {U_from_curve * DAY * 1e3:.3f} mm/day")
    print(f"      paper states 5.1 mm/day  -> ratio {U_from_curve * DAY * 1e3 / 5.1:.3f}")
    print("      (this uses rho = 1200 kg/m^3, so it is a second, independent")
    print("       consistency check on that back-calibrated value)")

    print("\n  (b) growth rates")
    s43 = float(np.sum(pub43[0] * pub43[1]) / np.sum(pub43[0] ** 2)) * 24.0
    s34 = float(np.sum(pub34[0] * pub34[1]) / np.sum(pub34[0] ** 2)) * 24.0
    print(f"      eps_s = 0.43:  published {s43:.4f}   ours {rate43 * DAY * 1e3:.4f} mm/day"
          f"   ratio {rate43 * DAY * 1e3 / s43:.4f}")
    print(f"      eps_s = 0.34:  published {s34:.4f}   ours {rate34 * DAY * 1e3:.4f} mm/day"
          f"   ratio {rate34 * DAY * 1e3 / s34:.4f}")
    print(f"      spread between the two porosities:"
          f"  published {s43 / s34:.4f}   ours {rate43 / rate34:.4f}")
    print("      -> our eps_s = 0.34 line matches well; the eps_s = 0.43 line is ~5% high,")
    print("         and the published pair are closer together than eq (32) predicts.")

    exp_slope = float(np.sum(t_ex * h_ex) / np.sum(t_ex**2)) * 24.0
    print(f"\n      experimental points give {exp_slope:.4f} mm/day through the origin")
    print(f"      (our eps_s = 0.34 prediction is "
          f"{abs(rate34 * DAY * 1e3 / exp_slope - 1) * 100:.1f}% from it)")


if __name__ == "__main__":
    main()
