"""Transferability of the Prat (2026) growth law across initial concentration.

Prat validates the law in Figure 8 against one experiment at one concentration:
the 6 Molal sand column of Norouzi Rad, Shokri & Sahimi, Phys. Rev. E 88,
032404 (2013), here ref 18. That paper also reports 4.0 M and 3.5 M columns,
which Prat does not use. This script evaluates the law on all three.

All model quantities come from prat2026/core.py and params.py unmodified. All
experimental quantities come from the digitised ref 18 panels held in
new_model/validation/data, which are read rather than copied so that a
re-digitisation propagates automatically.

Run:  ../../.venv/bin/python analyze.py
"""

import base64
import csv
import io
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import core                                                       # noqa: E402
from params import C_EQ, RHO, RHO_CR, NorouziRad                  # noqa: E402

DATA = HERE.parent.parent / "new_model" / "validation" / "data"

# --- experimental configuration, common to all three columns -----------------
EPS_S = NorouziRad.EPS_S        # 0.43, salt structure porosity
D_SUB = NorouziRad.D_SUB        # 0.58 mm, mean sand grain diameter
AREA = 95.033e-6                # m^2, from the 11 mm column diameter stated in
                                # ref 18. prat2026/fig08 recovered 94.6 mm^2 by
                                # back-calibration without access to ref 18.

# ref 18 Figure 3, 4.0 M, segmented in new_model/validation/fig03_coverage.py
COV_T = np.array([3.5, 7.5, 10.0, 15.0])
COV = np.array([0.051795, 0.303474, 0.450372, 0.628275])

SERIES = [(3.5, "3p5", "#d62728", "^"),
          (4.0, "4p0", "#1a7a3a", "D"),
          (6.0, "6p0", "#3333aa", "o")]

# window used for the 4.0 M thickness comparison: the 3.5 hr image has
# coverage 0.052 and a precipitated mass indistinguishable from the Fig 4a
# intercept, so a thickness cannot be resolved there.
WINDOW = (7.5, 15.0)


def load(name):
    rows = [r for r in csv.reader(open(DATA / name))
            if r and not r[0].lstrip().startswith("#")]
    a = np.array([[float(v) for v in r] for r in rows[1:]])
    return a[:, 0], a[:, 1]


def fit(x, y, origin=False):
    """(slope, intercept, slope_stderr); origin=True forces intercept 0."""
    if origin:
        s = float(x @ y / (x @ x))
        r = y - s * x
        return s, 0.0, float(np.sqrt((r @ r) / (len(x) - 1) / (x @ x)))
    A = np.vstack([x, np.ones_like(x)]).T
    (s, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - (s * x + b)
    sxx = float(((x - x.mean()) ** 2).sum())
    return float(s), float(b), float(np.sqrt((r @ r) / (len(x) - 2) / sxx))


def flux_curvature(t, mw):
    """Departure of M_w(t) from a straight line, in standard errors.

    A constant evaporative flux makes M_w(t) linear, so the quadratic
    coefficient a2 vanishes. a2/se(a2) measures how far the record departs
    from that. Returns (a2, se, sigma, dMw_dt at the first and last point).
    """
    A = np.vstack([t ** 2, t, np.ones_like(t)]).T
    c, *_ = np.linalg.lstsq(A, mw, rcond=None)
    r = mw - A @ c
    cov = (r @ r) / (len(t) - 3) * np.linalg.inv(A.T @ A)
    se = float(np.sqrt(cov[0, 0]))
    ends = (float(2 * c[0] * t.min() + c[1]), float(2 * c[0] * t.max() + c[1]))
    return float(c[0]), se, float(c[0]) / se, ends, c


def water_weighted_coverage(times, cov, t_end, t_mw, mw):
    """Mean coverage weighted by evaporated water, not by time.

    eta = dM_s/dM_w is the global slope of M_s against M_w, so it equals
    int eta(t) dM_w / int dM_w. The two averages coincide only for a constant
    flux, which does not hold at 4.0 M.
    """
    t = np.concatenate(([0.0], times, [t_end]))
    c = np.concatenate(([0.0], cov, [cov[-1]]))
    g = np.linspace(0.0, t_end, 4001)
    cg = np.interp(g, t, c)
    quad = np.polyfit(t_mw, mw, 2)
    w = np.clip(np.polyval(np.polyder(quad), g), 0.0, None)
    return (float(np.trapezoid(cg * w, g) / np.trapezoid(w, g)),
            float(np.trapezoid(cg, g) / t_end))


def c0_mass_fraction(molal):
    m = molal * 58.44e-3
    return m / (1 + m)


def png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# --- the model, evaluated once: it carries no C_0 dependence -----------------
DA_STAR = float(core.Da_star(eps_s=EPS_S, d_sub=D_SUB))
F_KIN = float(core.F(DA_STAR))
CLOSURE = C_EQ * F_KIN / (1.0 - C_EQ)       # dMs/dMw at eta = 1
CEILING = C_EQ / (1.0 - C_EQ)               # the same at F = 1


def main():
    rec = {}
    for molal, tag, colour, marker in SERIES:
        mw, ms = load(f"nr2013_fig04a_{tag}.csv")
        t, mwt = load(f"nr2013_mass_loss_{tag}.csv")

        s_f, b_f, se_f = fit(mw, ms)
        s_o, _, se_o = fit(mw, ms, origin=True)
        dmw_dt, mw_0, _ = fit(t, mwt)                       # g/hr, run average
        a2, a2_se, a2_sig, ends, quad = flux_curvature(t, mwt)

        j = dmw_dt * 1e-3 / 3600.0 / AREA                   # kg/m^2/s of water
        dh_dt = float(core.growth_rate(j, eps_s=EPS_S, d_sub=D_SUB))
        c0 = c0_mass_fraction(molal)

        rec[molal] = dict(
            tag=tag, colour=colour, marker=marker, n=len(mw),
            C0=c0, C0_over_Ceq=c0 / C_EQ, budget=c0 / (1.0 - c0),
            slope_free=s_f, intercept_free=b_f, slope_stderr=se_f,
            slope_origin=s_o, slope_origin_stderr=se_o,
            eta_free=s_f / CLOSURE, eta_free_stderr=se_f / CLOSURE,
            eta_origin=s_o / CLOSURE,
            eta_budget_max=(c0 / (1.0 - c0)) / CLOSURE,
            dMw_dt=dmw_dt, Mw_intercept=mw_0, j=j, quad=quad,
            flux_a2=a2, flux_a2_sigma=a2_sig, flux_ends=ends,
            flux_ratio=ends[1] / ends[0],
            U_evap_mm_day=(j / RHO) * 86400e3,
            dh_dt_mm_day=dh_dt * 86400e3,
            t_span=[float(t.min()), float(t.max())],
            mw=mw, ms=ms, t=t, mwt=mwt,
        )

    r4, r6 = rec[4.0], rec[6.0]
    t_end = float(r4["t"].max())
    eta_water, eta_time = water_weighted_coverage(
        COV_T, COV, t_end, r4["t"], r4["mwt"])

    # --- eq (35) at 4.0 M, fed the measured coverage ------------------------
    def eq35(eta):
        return float(core.dMs_dt(r4["j"], AREA, eta,
                                 eps_s=EPS_S, d_sub=D_SUB)) * 3600e3   # g/hr
    dms_exp = r4["slope_free"] * r4["dMw_dt"]
    check = dict(eta_water=eta_water, eta_time=eta_time,
                 dMs_dt_exp_g_hr=dms_exp,
                 dMs_dt_model_water_g_hr=eq35(eta_water),
                 dMs_dt_model_time_g_hr=eq35(eta_time),
                 dev_water_pct=(eq35(eta_water) / dms_exp - 1) * 100,
                 dev_time_pct=(eq35(eta_time) / dms_exp - 1) * 100,
                 sigma_water=abs(eta_water - r4["eta_free"])
                 / r4["eta_free_stderr"],
                 sigma_time=abs(eta_time - r4["eta_free"])
                 / r4["eta_free_stderr"])

    # --- 4.0 M: window-matched model rate, and the thickening/spreading split
    # Ms = rho_cr (1-eps_s) A (eta h), so the data constrain the product eta*h.
    # d(eta h) = eta_bar dh + h_bar d_eta is exact for midpoint averages.
    g = np.linspace(*WINDOW, 400)
    j_window = float(np.polyval(np.polyder(r4["quad"]), g).mean()) \
        * 1e-3 / 3600.0 / AREA
    dh_dt_window = float(core.growth_rate(j_window, eps_s=EPS_S,
                                          d_sub=D_SUB)) * 86400e3

    def reconstruct(slope, intercept):
        mw_t = np.polyval(r4["quad"], COV_T)                    # g
        ms_t = slope * mw_t + intercept                         # g
        eta_h = ms_t * 1e-3 / (RHO_CR * (1 - EPS_S) * AREA)     # m
        return mw_t, ms_t, eta_h * 1e3, (eta_h / COV) * 1e3     # mm, mm

    recon = {}
    for lab, (s, b) in (("free", (r4["slope_free"], r4["intercept_free"])),
                        ("origin", (r4["slope_origin"], 0.0))):
        mw_t, ms_t, eta_h, h = reconstruct(s, b)
        k = COV_T >= WINDOW[0]
        slope_h = float(np.polyfit(COV_T[k], h[k], 1)[0] * 24.0)   # mm/day
        seg = []
        for i in range(1, len(COV_T)):
            d_etah = eta_h[i] - eta_h[i - 1]
            eb, hb = 0.5 * (COV[i] + COV[i - 1]), 0.5 * (h[i] + h[i - 1])
            thick, spread = eb * (h[i] - h[i - 1]), hb * (COV[i] - COV[i - 1])
            seg.append(dict(t0=float(COV_T[i - 1]), t1=float(COV_T[i]),
                            d_etah=float(d_etah), thicken=float(thick),
                            spread=float(spread),
                            spread_share=float(spread / d_etah)))
        tot_t = sum(s_["thicken"] for s_ in seg[1:])
        tot_s = sum(s_["spread"] for s_ in seg[1:])
        recon[lab] = dict(Mw=mw_t, Ms=ms_t, eta_h=eta_h, h=h,
                          slope_h_mm_day=slope_h,
                          model_over_exp=dh_dt_window / slope_h,
                          segments=seg,
                          spread_share=float(tot_s / (tot_s + tot_t)))

    figs = make_figures(rec, recon, dh_dt_window, eta_water, eta_time)
    metrics = dict(
        model=dict(Da_star=DA_STAR, F=F_KIN, closure=CLOSURE, ceiling=CEILING,
                   eps_s=EPS_S, d_sub_m=D_SUB, area_m2=AREA,
                   note="Da*, F and the closure carry no C_0 dependence"),
        series={f"{m:.1f}": {k: v for k, v in rec[m].items()
                             if k not in ("mw", "ms", "t", "mwt", "quad")}
                for m, *_ in SERIES},
        coverage_4M=dict(times_hr=COV_T.tolist(), coverage=COV.tolist(),
                         eta_water_weighted=eta_water,
                         eta_time_weighted=eta_time),
        eq35_check_4M=check,
        thickness_4M=dict(window_hr=list(WINDOW), j_window=j_window,
                          dh_dt_model_window_mm_day=dh_dt_window,
                          **{k: {kk: vv for kk, vv in v.items()
                                 if kk not in ("Mw", "Ms", "eta_h", "h")}
                             for k, v in recon.items()}),
    )
    (HERE / "metrics.json").write_text(json.dumps(metrics, indent=2))
    (HERE / "report.html").write_text(
        render(rec, recon, check, figs, dh_dt_window, eta_water, eta_time),
        encoding="utf-8")

    print(f"model:  Da* = {DA_STAR:.1f}  F = {F_KIN:.4f}  "
          f"closure = {CLOSURE:.4f}  (no C_0 anywhere)")
    print(f"{'molal':>6}{'C0/Ceq':>8}{'j run-avg':>11}{'a2 sigma':>10}"
          f"{'flux end/start':>15}{'dh/dt model':>12}{'eta req':>9}"
          f"{'eta budget':>11}")
    for m, *_ in SERIES:
        r = rec[m]
        print(f"{m:>6.1f}{r['C0_over_Ceq']:>8.3f}{r['j']:>11.3e}"
              f"{r['flux_a2_sigma']:>+10.2f}{r['flux_ratio']:>15.2f}"
              f"{r['dh_dt_mm_day']:>12.3f}{r['eta_free']:>9.3f}"
              f"{r['eta_budget_max']:>11.3f}")
    print(f"\n4.0 M imaged coverage, water-weighted {eta_water:.3f} "
          f"({check['sigma_water']:.2f} sigma from the Fig 4a slope); "
          f"time-weighted {eta_time:.3f} ({check['sigma_time']:.2f} sigma)")
    print(f"eq (35) fed the water-weighted coverage: "
          f"{check['dMs_dt_model_water_g_hr']:.5f} vs experiment "
          f"{check['dMs_dt_exp_g_hr']:.5f} g/hr "
          f"({check['dev_water_pct']:+.1f}%)")
    print(f"\n4.0 M thickness over {WINDOW[0]}-{WINDOW[1]} hr "
          f"(model, window-matched flux: {dh_dt_window:.2f} mm/day)")
    for lab, v in recon.items():
        print(f"  {lab:>6} fit: {v['slope_h_mm_day']:.2f} mm/day, "
              f"model/exp {v['model_over_exp']:.2f}x, spreading share "
              f"{v['spread_share'] * 100:.0f}%")
    print(f"\nreport:  {HERE / 'report.html'}")


# =============================================================================
# Figures. Labels are in English: no CJK font is available to matplotlib here.
# =============================================================================
def make_figures(rec, recon, dh_window, eta_water, eta_time):
    figs = {}

    # --- F1: the flux, and the salt yield -----------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.7))
    for m, tag, c, mk in SERIES:
        r = rec[m]
        a1.plot(r["t"], r["mwt"], mk, color=c, ms=3.4, ls="none", alpha=0.7)
        tt = np.linspace(max(0, r["t"].min()), r["t"].max(), 100)
        a1.plot(tt, np.polyval(r["quad"], tt), "-", color=c, lw=1.9,
                label=f"{m:.1f} M  quadratic fit,  $a_2$ = "
                      f"{r['flux_a2_sigma']:+.1f}$\\sigma$,  "
                      f"flux {r['flux_ratio']:.2f}$\\times$ over the run")
        a1.plot(tt, r["dMw_dt"] * tt + r["Mw_intercept"], ":", color=c, lw=1.1)
    a1.set_xlabel("time  [hr]")
    a1.set_ylabel("cumulative evaporated water  $M_w$  [g]")
    a1.set_title("(a) ref 18 Fig 2: only the 6.0 M flux is constant\n"
                 "dotted = linear fit, whose slope is the run average",
                 fontsize=10)
    a1.legend(fontsize=7.6, loc="upper left", frameon=False)
    a1.grid(alpha=0.3)

    xx = np.linspace(0, 0.85, 50)
    a2.plot(xx, CLOSURE * xx, "-", color="0.15", lw=2.6, zorder=3,
            label=f"model, $\\eta$ = 1 (slope {CLOSURE:.4f}); independent of $C_0$")
    a2.plot(xx, CEILING * xx, ":", color="0.55", lw=1.4,
            label=f"saturated-brine ceiling  {CEILING:.4f}")
    for m, tag, c, mk in SERIES:
        r = rec[m]
        a2.plot(r["mw"], r["ms"], mk, color=c, ms=5.5, ls="none", zorder=5,
                label=f"{m:.1f} M   slope {r['slope_free']:.4f}   "
                      f"$\\eta$ = {r['eta_free']:.3f}")
        a2.plot(xx, r["budget"] * xx, "--", color=c, lw=1.2, alpha=0.85,
                zorder=2)
        a2.annotate(f"$C_0/(1-C_0)$ = {r['budget']:.3f}",
                    xy=(0.845, r["budget"] * 0.845), fontsize=7.4, color=c,
                    ha="right", va="bottom")
    a2.set_xlim(0, 0.87)
    a2.set_ylim(0, 0.32)
    a2.set_xlabel("cumulative evaporated water  $M_w$  [g]")
    a2.set_ylabel("cumulative precipitated salt  $M_s$  [g]")
    a2.set_title("(b) ref 18 Fig 4a: one model line against three columns\n"
                 "dashed = that column's own salt inventory", fontsize=10)
    a2.legend(fontsize=7.4, loc="upper left", frameon=False)
    a2.grid(alpha=0.3)
    fig.tight_layout()
    figs["f1"] = png(fig)

    # --- F2: eta against C_0, with the conservation bound -------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.5))
    ratios = [rec[m]["C0_over_Ceq"] for m, *_ in SERIES]
    free = [rec[m]["eta_free"] for m, *_ in SERIES]
    ferr = [rec[m]["eta_free_stderr"] for m, *_ in SERIES]
    orig = [rec[m]["eta_origin"] for m, *_ in SERIES]

    xs_c = np.linspace(0.60, 1.02, 100)
    cx = xs_c * C_EQ
    a1.fill_between(xs_c, np.minimum((cx / (1 - cx)) / CLOSURE, 1.35), 1.35,
                    color="#ffcdd2", alpha=0.55, zorder=0,
                    label="excluded: more salt than the column contains")
    a1.plot(xs_c, (cx / (1 - cx)) / CLOSURE, "-", color="#c62828", lw=1.8,
            zorder=2, label="conservation bound on $\\eta$")
    a1.axhline(1.0, color="0.15", lw=2.4, zorder=3,
               label="model, $\\eta$ = 1 (structural, not fitted)")
    a1.errorbar(ratios, free, yerr=ferr, fmt="o-", color="#1f77b4", lw=2,
                ms=9, capsize=4, zorder=5, label="required by Fig 4a, free slope")
    a1.plot(ratios, orig, "s--", color="#9467bd", lw=1.5, ms=7, alpha=0.85,
            zorder=5, label="required by Fig 4a, through origin")
    for (m, tag, c, mk), rr, ff in zip(SERIES, ratios, free):
        a1.annotate(f"{m:.1f} M", xy=(rr, ff), xytext=(0, 13),
                    textcoords="offset points", ha="center", fontsize=9,
                    color=c, weight="bold")
    a1.set_xlim(0.60, 1.02)
    a1.set_ylim(0, 1.35)
    a1.set_xlabel("initial saturation  $C_0/C_{eq}$")
    a1.set_ylabel("$\\eta$")
    a1.set_title("(a) the coverage required by the data varies with $C_0$",
                 fontsize=10)
    a1.legend(fontsize=7.2, loc="lower right", frameon=False)
    a1.grid(alpha=0.3)

    xs = np.arange(3)
    lo = [CLOSURE / rec[m]["slope_free"] for m, *_ in SERIES]
    hi = [CLOSURE / rec[m]["slope_origin"] for m, *_ in SERIES]
    pr = [0.76 * CLOSURE / rec[m]["slope_free"] for m, *_ in SERIES]
    a2.bar(xs - 0.26, lo, 0.25, color="#1f77b4", label="$\\eta$=1, free slope")
    a2.bar(xs, hi, 0.25, color="#9467bd", label="$\\eta$=1, through origin")
    a2.bar(xs + 0.26, pr, 0.25, color="#e07b39",
           label="$\\eta$=0.76, the value fitted at 6 M")
    for i, (p, q, s) in enumerate(zip(lo, hi, pr)):
        for dx, v in ((-0.26, p), (0.0, q), (0.26, s)):
            a2.text(i + dx, v + 0.06, f"{v:.2f}", ha="center", fontsize=7.6)
    a2.axhline(1.0, color="0.15", lw=2, zorder=3)
    a2.axhspan(0.9, 1.1, color="0.75", alpha=0.4, zorder=0, label="$\\pm$10%")
    a2.set_xticks(xs)
    a2.set_xticklabels([f"{m:.1f} M\n$C_0/C_{{eq}}$={rec[m]['C0_over_Ceq']:.2f}"
                        for m, *_ in SERIES], fontsize=9)
    a2.set_ylabel("model / experiment,  salt per evaporated water")
    a2.set_ylim(0, 4.4)
    a2.set_title("(b) ratio of predicted to observed salt yield", fontsize=10)
    a2.legend(fontsize=7.4, loc="upper right", frameon=False)
    a2.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    figs["f2"] = png(fig)

    # --- F3: 4.0 M coverage, thickness, and the split -----------------------
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15.6, 4.5))
    a1.plot(COV_T, COV, "o-", color="#1a7a3a", ms=8, lw=2,
            label="imaged coverage $\\eta(t)$, ref 18 Fig 3")
    a1.axhline(eta_water, color="#1a7a3a", ls="--", lw=1.8,
               label=f"run average, water-weighted {eta_water:.3f}")
    a1.axhline(eta_time, color="#1a7a3a", ls=":", lw=1.3, alpha=0.8,
               label=f"run average, time-weighted {eta_time:.3f}\n"
                     f"(equal only for a constant flux)")
    a1.axhspan(0.336 - 0.017, 0.336 + 0.017, color="#1f77b4", alpha=0.22,
               label="$\\eta$ required by Fig 4a, 0.336 $\\pm$ 0.017")
    a1.set_xlabel("time  [hr]")
    a1.set_ylabel("areal coverage  $\\eta$")
    a1.set_ylim(0, 0.72)
    a1.set_title("(a) 4.0 M coverage rises throughout the run", fontsize=10)
    a1.legend(fontsize=7.2, loc="upper left", frameon=False)
    a1.grid(alpha=0.3)

    h0 = recon["free"]["h"][1]
    tt = np.linspace(WINDOW[0], 16, 10)
    a2.plot(tt, dh_window * (tt - WINDOW[0]) / 24 + h0, "-", color="0.15",
            lw=2.4, label=f"model, eq (32), window-matched flux "
                          f"{dh_window:.2f} mm/day")
    for lab, st, c in (("free", "o-", "#1a7a3a"), ("origin", "s--", "#7fb069")):
        v = recon[lab]
        k = COV_T >= WINDOW[0]
        a2.plot(COV_T[k], v["h"][k], st, color=c, ms=8, lw=1.8,
                label=f"reconstructed, {lab} fit  "
                      f"{v['slope_h_mm_day']:.2f} mm/day")
    t6, h6 = load("fig08b_exp.csv")
    a2.plot(t6, h6, "o", color="#3333aa", ms=4.5, alpha=0.6, ls="none",
            label="6.0 M $h_{equ}$, measured (ref 18 Fig 4b)")
    a2.set_xlim(4, 22)
    a2.set_ylim(0, 2.2)
    a2.set_xlabel("time  [hr]")
    a2.set_ylabel("thickness  $h=M_s/[\\rho_{cr}(1-\\epsilon_s)\\eta A]$  [mm]")
    a2.set_title("(b) reconstructed 4.0 M thickness, against 6.0 M and the model",
                 fontsize=10)
    a2.legend(fontsize=7.2, loc="upper left", frameon=False)
    a2.grid(alpha=0.3)

    labels, thick, spread = [], [], []
    for s_ in recon["free"]["segments"][1:]:
        labels.append(f"{s_['t0']:.1f}–{s_['t1']:.0f} hr")
        thick.append(s_["thicken"])
        spread.append(s_["spread"])
    xs = np.arange(len(labels))
    a3.bar(xs, thick, 0.55, color="#1f77b4",
           label="thickening  $\\bar\\eta\\,\\Delta h$")
    a3.bar(xs, spread, 0.55, bottom=thick, color="#e07b39",
           label="lateral spreading  $\\bar h\\,\\Delta\\eta$")
    for i, (p, q) in enumerate(zip(thick, spread)):
        a3.text(i, p + q + 0.004, f"spreading {q / (p + q) * 100:.0f}%",
                ha="center", fontsize=8.5)
    a3.set_xticks(xs)
    a3.set_xticklabels(labels, fontsize=9)
    a3.set_ylabel("$\\Delta(\\eta h)$  [mm],  proportional to $\\Delta M_s$")
    a3.set_title("(c) 4.0 M: partition of the newly precipitated salt\n"
                 "$\\Delta(\\eta h)=\\bar\\eta\\Delta h+\\bar h\\Delta\\eta$, "
                 "an identity", fontsize=10)
    a3.legend(fontsize=8.4, loc="upper left", frameon=False)
    a3.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    figs["f3"] = png(fig)
    return figs


# =============================================================================
# Report
# =============================================================================
def render(rec, recon, check, figs, dh_window, eta_water, eta_time):
    def row(m):
        r = rec[m]
        return (f"<tr><td>{m:.1f} M</td><td>{r['C0']:.4f}</td>"
                f"<td>{r['C0_over_Ceq']:.3f}</td><td>{r['dMw_dt']:.5f}</td>"
                f"<td>{r['flux_a2_sigma']:+.2f}</td>"
                f"<td>{r['flux_ends'][0]:.4f} → {r['flux_ends'][1]:.4f}</td>"
                f"<td>{r['j']:.3e}</td><td>{r['U_evap_mm_day']:.2f}</td>"
                f"<td>{r['dh_dt_mm_day']:.3f}</td></tr>")

    def row2(m):
        r = rec[m]
        cls = "over" if r["eta_budget_max"] < 1 else ""
        return (f"<tr><td>{m:.1f} M</td><td>{r['n']}</td>"
                f"<td>{r['slope_free']:.4f}</td>"
                f"<td>{r['eta_free']:.3f} ± {r['eta_free_stderr']:.3f}</td>"
                f"<td>{r['eta_origin']:.3f}</td><td>{r['budget']:.4f}</td>"
                f"<td class='{cls}'>{r['eta_budget_max']:.3f}</td>"
                f"<td>{CLOSURE / r['slope_free']:.2f}</td></tr>")

    f, o = recon["free"], recon["origin"]
    r4, r6 = rec[4.0], rec[6.0]
    t6, h6 = load("fig08b_exp.csv")
    fit6 = np.polyfit(t6, h6, 1)
    sigma64 = ((r6["eta_free"] - r4["eta_free"])
               / float(np.hypot(r6["eta_free_stderr"], r4["eta_free_stderr"])))

    return TEMPLATE.format(
        f1=figs["f1"], f2=figs["f2"], f3=figs["f3"],
        rows="\n".join(row(m) for m, *_ in SERIES),
        rows2="\n".join(row2(m) for m, *_ in SERIES),
        da=DA_STAR, F=F_KIN, closure=CLOSURE, ceiling=CEILING,
        area=AREA * 1e6, eps=EPS_S, dsub=D_SUB * 1e3,
        eta6=r6["eta_free"], eta4=r4["eta_free"], eta35=rec[3.5]["eta_free"],
        ratio64=r6["eta_free"] / r4["eta_free"], sigma64=sigma64,
        dh35=rec[3.5]["dh_dt_mm_day"], dh4=r4["dh_dt_mm_day"],
        dh6=r6["dh_dt_mm_day"],
        jspread=(max(rec[m]["j"] for m, *_ in SERIES)
                 / min(rec[m]["j"] for m, *_ in SERIES) - 1) * 100,
        sig35=rec[3.5]["flux_a2_sigma"], sig4=r4["flux_a2_sigma"],
        sig6=r6["flux_a2_sigma"], fr35=rec[3.5]["flux_ratio"],
        fr4=r4["flux_ratio"], fr6=r6["flux_ratio"],
        u6=r6["U_evap_mm_day"],
        etaw=eta_water, etat=eta_time,
        sig_w=check["sigma_water"], sig_t=check["sigma_time"],
        dms_e=check["dMs_dt_exp_g_hr"], dms_w=check["dMs_dt_model_water_g_hr"],
        dev_w=check["dev_water_pct"], dev_t=check["dev_time_pct"],
        w0=WINDOW[0], w1=WINDOW[1], dhw=dh_window,
        hf=f["slope_h_mm_day"], ho=o["slope_h_mm_day"],
        rf=f["model_over_exp"], ro=o["model_over_exp"],
        h4_15=f["h"][-1], h4o_15=o["h"][-1],
        h6_15=float(np.polyval(fit6, 15.0)), dh6_exp=float(fit6[0] * 24.0),
        r6_ratio=r6["dh_dt_mm_day"] / float(fit6[0] * 24.0),
        sp_f=f["spread_share"] * 100, sp_o=o["spread_share"] * 100,
        budget4=r4["eta_budget_max"], budget35=rec[3.5]["eta_budget_max"],
        budget6=r6["eta_budget_max"], budget6raw=r6["budget"],
        degen=(1 - CLOSURE / r6["budget"]) * 100,
        sep4=(CLOSURE / r4["budget"] - 1) * 100,
        cross=CLOSURE / (1 + CLOSURE) / C_EQ,
        crossm=(CLOSURE / (1 + CLOSURE)) / (1 - CLOSURE / (1 + CLOSURE))
        / 58.44e-3,
        eta6o=r6["eta_origin"], eta6flat=r6["eta_origin"] / r6["eta_free"],
        covfin=COV[-1], cov_ratio=COV[-1] / eta_water,
    )


TEMPLATE = """<meta charset="utf-8">
<title>Prat (2026) 生长律的浓度可迁移性</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ max-width: 1040px; margin: 0 auto; padding: 2.4rem 1.6rem 6rem;
  font: 16px/1.8 -apple-system, "Segoe UI", "Noto Sans CJK SC", sans-serif; }}
h1 {{ font-size: 1.7rem; line-height: 1.35; margin-bottom: .3rem; }}
h2 {{ margin-top: 2.8rem; padding-bottom: .35rem; border-bottom: 1px solid #8884;
  font-size: 1.22rem; }}
h3 {{ margin-top: 1.9rem; font-size: 1.02rem; font-weight: 600; }}
.sub {{ color: #7a7a7a; margin-top: 0; font-size: .94rem; }}
table {{ border-collapse: collapse; width: 100%; margin: 1.1rem 0;
  font-size: .86rem; }}
th, td {{ border: 1px solid #8883; padding: .38rem .55rem; text-align: right; }}
th {{ background: #8881; font-weight: 600; text-align: center; }}
td:first-child, th:first-child {{ text-align: left; }}
.over {{ color: #c62828; }}
figure {{ margin: 1.6rem 0; }}
figure img {{ width: 100%; border: 1px solid #8883; border-radius: 3px; }}
figcaption {{ font-size: .84rem; color: #7a7a7a; margin-top: .5rem;
  line-height: 1.65; }}
.note {{ border-left: 3px solid #8886; padding: .1rem 0 .1rem 1rem;
  margin: 1.2rem 0; color: #666; font-size: .93rem; }}
code {{ background: #8882; padding: .1rem .3rem; border-radius: 3px;
  font-size: .87em; }}
.eq {{ text-align: center; margin: 1.1rem 0; }}
ol li, ul li {{ margin: .35rem 0; }}
.tc {{ text-align: center; }}
</style>

<h1>Prat (2026) 生长律的浓度可迁移性</h1>
<p class="sub">对 Norouzi Rad, Shokri &amp; Sahimi (Phys. Rev. E <b>88</b>,
032404, 2013) 三个砂柱（3.5 / 4.0 / 6.0 Molal NaCl）的重新评估</p>

<h2>摘要</h2>
<p>Prat (2026) 的蒸发驱动生长律在其 Figure 8 中只对一个浓度（6.0 M）作过检验。
本文把同一条律应用于该实验的全部三个浓度。生长律的表达式中只含溶解度
C<sub>eq</sub>，不含初始浓度 C₀，因此它预言单位蒸发水的产盐量与 C₀ 无关。
实测值与此不符：三个柱要求的表面覆盖率 η 分别为 {eta35:.3f}、{eta4:.3f} 和
{eta6:.3f}，6.0 M 与 4.0 M 相差 {sigma64:.1f} 个标准误。柱内盐存量给出一条与模型
无关的上界，该上界排除了 C₀/C<sub>eq</sub> &lt; {cross:.3f} 时的 η = 1。在
6.0 M 上，模型闭合系数与该守恒上界相差 {degen:.0f}%，两者在本数据上不可分辨。
4.0 M 的成像覆盖率在全程持续上升，新增盐中约有一半进入横向铺展项，而生长律只描述
增厚项。结论是：局部生长律本身不含浓度依赖，可以跨浓度使用；Figure 8 采用的整柱积分
形式则需要逐浓度提供 η，模型对此没有预言能力。</p>

<h2>1. 背景与问题</h2>
<p>Prat (2026) 将盐结构界面的法向生长速度写为局部蒸发通量的线性函数：</p>
<p class="eq">dh/dt = C<sub>eq</sub> j / [ρ<sub>cr</sub>(1−ε<sub>s</sub>)(1−C<sub>eq</sub>)] · F(Da*)</p>
<p>其 Figure 8 以 ref 18 的 6.0 Molal 砂柱作为检验对象，报告的一致性为 0.08%
（相对于论文自身给出的模型线）。ref 18 同时报告了 4.0 M 和 3.5 M 两个柱，Prat
未使用。本文检验该律在这两个浓度上是否同样成立。</p>

<h2>2. 数据与方法</h2>
<h3>2.1 实验数据</h3>
<p>全部实验数据由脚本从 600 dpi 位图数字化，不是原始数据。本文读取
<code>new_model/validation/data/</code> 下的产物。各面板的浓度覆盖并不对称，
这一点决定了后文可以检验什么、不能检验什么：</p>
<table>
<tr><th>面板</th><th class="tc">内容</th><th class="tc">原文覆盖的浓度</th><th class="tc">已数字化</th><th class="tc">本文用途</th></tr>
<tr><td>Fig 2</td><td class="tc">M<sub>w</sub>(t)，累计蒸发水</td><td class="tc">3.5 / 4.0 / 6.0 M</td><td class="tc">全部</td><td class="tc">确定蒸发通量 j</td></tr>
<tr><td>Fig 4(a)</td><td class="tc">M<sub>s</sub> 对 M<sub>w</sub></td><td class="tc">3.5 / 4.0 / 6.0 M</td><td class="tc">全部</td><td class="tc">确定实测产盐率，反推 η</td></tr>
<tr><td>Fig 4(b)</td><td class="tc">M<sub>s</sub>(t)，对数纵轴</td><td class="tc">6.0 / 4.0 M，3.5 M 在 inset</td><td class="tc">仅 6.0 M</td><td class="tc">6.0 M 的 η 时域交叉校验</td></tr>
<tr><td>Fig 3</td><td class="tc">盐层厚度图，4 个时刻</td><td class="tc">仅 4.0 M</td><td class="tc">是</td><td class="tc">直接成像的覆盖率 η(t)</td></tr>
</table>
<p class="note">Prat 的 Figure 8(b) 画的 h<sub>equ</sub>(t) 是由 M<sub>s</sub>
换算而来的导出量，只对 6.0 M 给出。ref 18 的 Fig 4(b) 本身则含 4.0 M 与 3.5 M，
本仓库尚未提取，见 5.6 节。</p>

<h3>2.2 浓度口径</h3>
<p>ref 18 的浓度以质量摩尔浓度给出（"3.5 M (moles of NaCl/kg of water), 4 M,
and 6 M"），对应质量分数 C₀ = 0.1698 / 0.1895 / 0.2596。原文采用的溶解度为
25 °C 下的 6.12 molal，折合质量分数 0.2634，与 Prat 取的
C<sub>eq</sub> = 0.264 相差 0.2%，两者可直接互用。本文的
C₀/C<sub>eq</sub> 均按质量分数计；若按质量摩尔浓度计则为
0.572 / 0.654 / 0.980。</p>

<h3>2.3 数字化的自校验</h3>
<ul>
<li>Fig 8a 与 8b 是同一实验的两种单位。两次独立提取的 13 个时间点最大偏差
0.049 hr，由二者反推的 ηA 在 13 个点上相对标准差为 0.29%。</li>
<li>6.0 M 的 η 可由两条路径求得：Fig 4(a) 的全局斜率给 {eta6:.3f}；分别拟合
M<sub>w</sub>(t) 与 M<sub>s</sub>(t) 再相除的时域路径给 0.759。两者相差 0.9%，
说明全局斜率法可以外推到没有 Fig 4(b) 面板的 3.5 M 和 4.0 M。</li>
</ul>

<h3>2.4 柱截面积</h3>
<p>式(35) 需要柱截面积 A，Prat 未给出。prat2026 在尚未取得 ref 18 的条件下，由
Fig 8a 与 8b 的比值反推得 94.6 mm²。ref 18 原文的实验条件为
"cylindrical plastic columns, 35 mm in height and 11 mm in diameter"，
对应 A = {area:.3f} mm²，与反推值相差 0.4%。本文采用原文值；该反推现已由
一手来源确认，不再是推断。</p>

<h3>2.5 模型侧</h3>
<p>全部模型量由 <code>prat2026/core.py</code> 原样计算，参数取自
<code>prat2026/params.py</code>。三个柱使用同一砂（d<sub>sub</sub> = {dsub:.2f} mm）、
同一结构孔隙率（ε<sub>s</sub> = {eps}）和同一温度，因此
Da* = {da:.1f} 和 F = {F:.4f} 三档相同。由此得到闭合系数</p>
<p class="eq">dM<sub>s</sub>/dM<sub>w</sub> = η · C<sub>eq</sub>F(Da*)/(1−C<sub>eq</sub>) = η × {closure:.4f}</p>
<p>该式中不出现 C₀。生长律的物理设定是盐结构由基质经毛细供给<b>饱和</b>盐水，
表面过饱和度为 O(10⁻⁵)，本体溶液成分不进入界面物质平衡。因此浓度依赖只可能由
η 承载，问题化为：η 是否为材料常数。</p>

<h2>3. 结果</h2>

<h3>3.1 蒸发通量：只有 6.0 M 是恒定的</h3>
<p>对 M<sub>w</sub>(t) 作二次拟合。蒸发通量恒定时曲线应为直线，二次系数 a₂ 为零；
a₂ 与其标准误之比给出偏离直线的显著性。</p>
<table>
<tr><th>浓度</th><th>C₀</th><th>C₀/C<sub>eq</sub></th><th>dM<sub>w</sub>/dt 跑均 [g/hr]</th><th>a₂ 显著性</th><th>dM<sub>w</sub>/dt 起 → 止 [g/hr]</th><th>j 跑均 [kg/m²/s]</th><th>U<sub>evap</sub> [mm/day]</th><th>模型 dh/dt [mm/day]</th></tr>
{rows}
</table>
<p>三个柱的<b>跑均</b>通量相差 {jspread:.0f}%，但只有 6.0 M 的通量在全程恒定
（a₂ 为 {sig6:+.2f}σ，与零一致）。3.5 M 在实验过程中减速 {fr35:.2f} 倍，
4.0 M 加速 {fr4:.2f} 倍。因此表中的模型 dh/dt 是跑均量，只有 6.0 M 一行可以
直接与全程的线性生长比较。</p>
<p>ref 18 原文称"The nearly constant slope of the curves indicates that in all
cases the process was in the stage-1 evaporation"，但同时说明各柱之间蒸发量的差别
"might be due to the lower atmospheric demand&hellip; the relative humidity and
temperature inside the x-ray chamber could not be controlled precisely"。
柱内的曲率很可能同源，即环境漂移而非浓度效应。这不影响 3.4 节的加权，
因为那里需要的是实测 dM<sub>w</sub>，与曲率成因无关；但它意味着不应把上表中
a₂ 的符号差异读作浓度的物理效应。</p>
<p>6.0 M 的 U<sub>evap</sub> 为 {u6:.2f} mm/day，与 Prat 正文引用的 7.7 mm/day
一致，说明数字化到 j 的整条计算链自洽。</p>

<figure>
<img src="data:image/png;base64,{f1}">
<figcaption><b>图 1.</b>（a）三个柱的蒸发曲线与二次拟合；虚线为线性拟合，其斜率即
跑均通量。（b）同一条模型线对三组实验数据；各柱的虚线为该柱自身盐存量给出的守恒界
C₀/(1−C₀)。</figcaption>
</figure>

<h3>3.2 要求的覆盖率随浓度变化</h3>
<table>
<tr><th>浓度</th><th>点数</th><th>实测 dM<sub>s</sub>/dM<sub>w</sub></th><th>要求的 η（自由拟合）</th><th>η（过原点）</th><th>C₀/(1−C₀)</th><th>守恒允许的 η<sub>max</sub></th><th>模型(η=1)/实验</th></tr>
{rows2}
</table>
<p>6.0 M 与 4.0 M 之间的差为 {sigma64:.1f}σ。若沿用 Prat 在 6.0 M 上拟合的
η = 0.76，模型对 4.0 M 的产盐率高估约 {ratio64:.1f} 倍。3.5 M 与 4.0 M 之间的差为
−0.043 ± 0.032，仅 1.32σ 且符号与"η 随 C₀ 单调上升"相反，两者在统计上打平。</p>

<h3>3.3 一条与模型无关的守恒界</h3>
<p>全程干燥后柱内的盐全部析出、水全部蒸发，因此跑均产盐率满足</p>
<p class="eq">dM<sub>s</sub>/dM<sub>w</sub> ≤ C₀/(1−C₀)</p>
<p>模型的闭合系数 {closure:.4f} 与 C₀ 无关。两者相交给出：η = 1 在
C₀/C<sub>eq</sub> &lt; {cross:.3f}（molality &lt; {crossm:.2f} m）时超出柱内盐总量。
对应的上界为 4.0 M 的 η ≤ {budget4:.3f}、3.5 M 的 η ≤ {budget35:.3f}；6.0 M 的
上界 {budget6:.3f} 大于 1，不构成约束。</p>
<p>在 6.0 M 上，模型闭合系数 {closure:.4f} 与该柱的守恒界 {budget6raw:.4f} 相差
{degen:.0f}%。两条预言在该浓度的数据精度下不可分辨。在 4.0 M 上两者相差
{sep4:.0f}%，具备分辨能力。</p>

<figure>
<img src="data:image/png;base64,{f2}">
<figcaption><b>图 2.</b>（a）红色区域为守恒排除区；模型的水平线 η = 1 在
C₀/C<sub>eq</sub> &lt; {cross:.3f} 处进入该区域。（b）预测与实测产盐率之比。</figcaption>
</figure>

<h3>3.4 4.0 M：η 与成像覆盖率的比较</h3>
<p>η = dM<sub>s</sub>/dM<sub>w</sub> 是全局斜率，等于按<b>蒸发水</b>加权的平均覆盖率，
而非按时间加权。两者只在通量恒定时相同，而 4.0 M 的通量在全程升高
{fr4:.2f} 倍（3.1 节）。</p>
<table>
<tr><th>加权方式</th><th>成像覆盖率跑均</th><th>与 Fig 4a 的 η = {eta4:.3f} ± 0.017 的偏离</th></tr>
<tr><td>按蒸发水（正确口径）</td><td>{etaw:.3f}</td><td>{sig_w:.2f}σ</td></tr>
<tr><td>按时间（此前采用）</td><td>{etat:.3f}</td><td>{sig_t:.2f}σ</td></tr>
</table>
<p>因此 Prat 将 η 认同为几何覆盖率这一定义，在 ref 18 自身数据上的支持程度为
约 {sig_w:.0f}σ 的张力，而非此前报告的 0.1σ。该位移完全由 a₂ 驱动：改用线性
M<sub>w</sub>(t) 会精确返回时间加权值。</p>
<p>把成像覆盖率代入式(35)：水加权值给 {dms_w:.5f} g/hr，实测为 {dms_e:.5f} g/hr，
偏离 {dev_w:+.1f}%。总量层面二者仍大致相符。</p>

<h3>3.5 4.0 M：新增盐在增厚与铺展之间的分配</h3>
<p>上节的 η 是面积份额，与厚度无关。厚度可由两条彼此独立的路径得到。</p>
<p><b>路径一，直接测量。</b>Fig 3 的色标以 mm 标定，因此该面板本身就是一个厚度场。
把原生栅格（自 PDF 提取，1020×934）的色标反演即得厚度。色标的口径由与 Fig 4(a) 的
质量闭合判定：按多孔壳包络高度读（ε<sub>s</sub> = 0.43）只有 Fig 4(a) 质量的 0.53，
按等效实心盐厚度读则为 0.93，故色标记录的是等效实心厚度，Prat 口径的结构高度为
h = T/(1−ε<sub>s</sub>)。详见 <code>new_model/validation/fig03_thickness.py</code>。</p>
<p><b>路径二，三面板还原。</b>由 Fig 4(a) 的 M<sub>s</sub>、Fig 2 的
M<sub>w</sub>(t) 和 Fig 3 的 η(t)，按 Prat 自身的定义
h = M<sub>s</sub>/[ρ<sub>cr</sub>(1−ε<sub>s</sub>)ηA] 反算。</p>
<p>比较窗口取 {w0:.1f}–{w1:.0f} hr：3.5 hr 一帧的覆盖率仅 0.043，其
M<sub>s</sub> 落在 Fig 4(a) 的截距噪声内，厚度不可分辨。模型侧的通量按同一窗口
重取（{dhw:.2f} mm/day），以避免用跑均通量对比窗口内速率。</p>
<table>
<tr><th>量</th><th>4.0 M 直接测量</th><th>4.0 M 还原（自由）</th><th>4.0 M 还原（过原点）</th><th>6.0 M 直接实测</th></tr>
<tr><td>15 hr 时的厚度 [mm]</td><td>0.491</td><td>{h4_15:.3f}</td><td>{h4o_15:.3f}</td><td>{h6_15:.3f}</td></tr>
<tr><td>增厚速率 [mm/day]</td><td><b>0.75</b></td><td>{hf:.2f}</td><td>{ho:.2f}</td><td>{dh6_exp:.2f}</td></tr>
<tr><td>模型 / 实验</td><td class="over"><b>3.31</b></td><td class="over">{rf:.2f}</td><td class="over">{ro:.2f}</td><td>{r6_ratio:.2f}</td></tr>
</table>
<p>直接测量落在还原路径给出的 {ho:.2f}–{hf:.2f} mm/day 区间内，两条路径互证。
因此 4.0 M 的增厚速率被模型高估 <b>3.3 倍</b>，这是实测结论而非推导结论。</p>
<p>式(35) 把 η 处理为常数。由 M<sub>s</sub> = ρ<sub>cr</sub>(1−ε<sub>s</sub>)A·ηh
求导可得一个恒等式：</p>
<p class="eq">dM<sub>s</sub>/dt = ρ<sub>cr</sub>(1−ε<sub>s</sub>)A [ η·dh/dt + h·dη/dt ]</p>
<p>生长律只对应第一项。在 4.0 M 的 {w0:.1f}–{w1:.0f} hr 窗口内，新增盐的
{sp_f:.0f}%（过原点拟合为 {sp_o:.0f}%）进入第二项。</p>
<p>6.0 M 的 η(t) 无法直接测得（见 5.3 节），但 ref 18 原文对这两类行为给出了
明确区分。对 3.5 M 与 4 M：<i>"the salt precipitation rate (the slope of the
line) initially changes with time, followed by a nearly constant value&hellip;
The transition point corresponds to the time at which all evaporation spots
serve as the precipitation spots"</i>；对 6 M：<i>"the salt deposition rate
remains nearly constant during the entire course of the experiment because the
initial salt concentration is very close to the solubility limit of NaCl, and
all the evaporation spots are quickly (within a few minutes) transformed into
precipitation spots after the onset of the experiment"</i>。这与上式的分解一致：
6.0 M 的 dη/dt 在实验开始数分钟内即归零，低浓度两档则存在一个贯穿实验的转变期。
Fig 4(a) 的自由与过原点拟合之比 {eta6o:.3f}/{eta6:.3f} = {eta6flat:.2f}
与该陈述方向一致。</p>

<figure>
<img src="data:image/png;base64,{f3}">
<figcaption><b>图 3.</b>（a）4.0 M 的成像覆盖率及其两种加权跑均。
（b）还原的 4.0 M 厚度、直接实测的 6.0 M 厚度与窗口匹配的模型线。
（c）新增盐在增厚项与铺展项之间的分配；分解式为恒等式而非拟合。</figcaption>
</figure>

<h2>4. 讨论</h2>
<p>三个浓度的差异不能归因于驱动通量：三者的跑均 j 相差 {jspread:.0f}%，模型据此
预言的增厚速率为 {dh35:.2f}、{dh4:.2f} 和 {dh6:.2f} mm/day。差异也不能归因于
Da* 或 F，因为三者共用同一砂、同一孔隙率和同一温度。在模型的现有表述中，
只剩 η 可以承载浓度依赖，而 η 是外部输入而非模型输出。</p>
<p>6.0 M 之所以能通过检验，可以由两条彼此独立的观察解释。其一，该浓度的
C₀/C<sub>eq</sub> = 0.983，模型闭合系数与柱内盐守恒界只差 {degen:.0f}%，两条预言
在该数据上简并，Figure 8 的一致性不能区分二者。其二，按 ref 18 原文的描述，
6 M 柱的全部蒸发斑在实验开始数分钟内即转为析出斑，dη/dt 随即归零，
"定面积增厚"因而是一个合理近似；3.5 M 与 4 M 则存在一个转变期，该近似不成立。
两条解释都指向同一结论：Figure 8 的检验条件比生长律的一般适用条件更严，
其结果不能外推。</p>
<p>就孔网模型耦合而言，这一分析支持把该律作为局部界面闭合使用，但要求由孔网模型
自身判定结壳底部是否确实处于 C<sub>eq</sub>，而不能预设。基质经达西流送达的是本体
浓度的溶液，3.3 节的守恒界表明低浓度下该预设难以长期维持。浓度依赖的物理落点是
成核与铺展的时空统计，即 dη/dt，而非生长律的系数；ref 18 原文也将 η 的浓度依赖
归因于成核转变。</p>

<h2>5. 局限与适用范围</h2>
<h3>5.1 误差数字必须连同浓度报告</h3>
<p>0.08% 与 −5.17% 均为 6.0 M 专属。同一闭合在 4.0 M 上的偏差为
{sep4:.0f}% 量级或更大，取决于所用的 η。</p>
<h3>5.2 4.0 M 厚度的两条路径与残余不确定度</h3>
<p>还原路径对 Fig 4(a) 的截距敏感：自由与过原点拟合给出的增厚速率相差 2 倍
（{hf:.2f} 与 {ho:.2f} mm/day）。直接测量路径不受此影响，给出 0.75 mm/day，
落在该区间内。因此"模型高估约 3.3 倍"可以引用；此前版本给出的
{ro:.1f}–{rf:.1f} 区间已由直接测量收紧。</p>
<p>直接测量的残余不确定度有三项：色标最大值由印刷刻度读出，是唯一不能自动化的输入；
薄盐颜色接近白色的空气相，彩度阈 25→80 时覆盖率从 0.046 变到 0.039；
9 mm 方形视野在几何上放不进 11 mm 圆柱（最大内接方 7.8 mm），而 ref 18 的 Fig 1
又把截面标作 10 mm，正文说 11 mm。质量按"视野代表全柱"缩放，改用 81 mm²
会让全部质量变 17%，不影响色标口径的判别。</p>
<h3>5.3 6.0 M 的 η(t) 无从直接测得</h3>
<p>成像面板 Fig 3 只覆盖 4.0 M。Prat 的 Figure 8(b) 是 Fig 8(a) 的单位换算，
h<sub>equ</sub> 的定义中已写入 η = 0.76，不构成独立观测。因此"6.0 M 覆盖率早期
饱和"这一判断的依据是 ref 18 正文的定性陈述（斑点在数分钟内全部转化）与
Fig 4(a) 截距的间接支持，而非覆盖率的直接测量。</p>
<h3>5.4 三点单调性不成立</h3>
<p>3.5 M 与 4.0 M 之间为 1.32σ 且符号相反，具备统计效力的只有 6.0 M 对低浓度两档的
{sigma64:.1f}σ 一步。这是两点关系，不应表述为三点单调趋势。</p>
<h3>5.5 水加权的敏感性</h3>
<p>3.4 节中 0.05 的 η 位移完全由 M<sub>w</sub>(t) 的二次项驱动，而该项本身是对
数字化面板的拟合（{sig4:+.2f}σ）。覆盖率一侧只有 4 个时刻，中间按线性插值、
15 hr 之后按常数外推。因此"约 {sig_w:.0f}σ 张力"应理解为方向明确（成像覆盖率高于
斜率反推的 η）、幅度待定。</p>
<h3>5.6 ref 18 Fig 4(b) 的低浓度面板尚未提取</h3>
<p>ref 18 的 Fig 4(b) 含 6.0 M 与 4.0 M 两条 M<sub>s</sub>(t) 曲线，3.5 M 在
inset 中，本仓库只数字化了 6.0 M。提取 4.0 M 一条可以把 3.5 节的厚度还原从
"Fig 4(a) 斜率 × Fig 2 的 M<sub>w</sub>(t)"改为直接读 M<sub>s</sub>(t)，
去掉一个环节，并直接检验原文所述的转变时刻。该面板为对数纵轴，提取难度高于
其余面板。这是目前收紧 3.5 节结论最直接的途径。</p>

<h2>6. 结论</h2>
<ol>
<li>生长律的局部形式（给定 j、单位已成结壳面积的 dh/dt）不含 C₀，可以跨浓度使用。
把实测覆盖率代入后，4.0 M 的产盐总量与预测相符到 {dev_w:+.1f}%。</li>
<li>Figure 8 采用的整柱积分形式不可跨浓度迁移，因为它需要一个逐浓度确定的 η，
而模型对 η 没有预言能力。沿用 6.0 M 拟合值会使 4.0 M 高估约 {ratio64:.1f} 倍。</li>
<li>柱内盐守恒给出一条与模型无关的上界，排除了 C₀/C<sub>eq</sub> &lt; {cross:.3f}
时的 η = 1；同一条界在 6.0 M 上与模型闭合系数简并到 {degen:.0f}%。</li>
<li>低浓度下新增盐的约一半进入横向铺展项，该项不在生长律的表述范围内。浓度依赖的
物理落点是 dη/dt。</li>
</ol>

<h2>附录：复现</h2>
<p><code>cd prat2026/multi_c &amp;&amp; ../../.venv/bin/python analyze.py</code>
重算全部数字、重绘三张图、重写本文件并输出 <code>metrics.json</code>。</p>
<p>模型侧取自 <code>prat2026/core.py</code> 与 <code>params.py</code>，未作修改；
实验侧取自 <code>new_model/validation/data/</code> 的数字化产物。3.4 节的水加权
覆盖率与 <code>new_model/validation/fig03_coverage.py</code> 使用同一算法。</p>
<p class="sub">M. Prat, <i>Langmuir</i> 2026, DOI 10.1021/acs.langmuir.6c02211 ·
Norouzi Rad, Shokri &amp; Sahimi, <i>Phys. Rev. E</i> <b>88</b>, 032404 (2013)</p>
"""


if __name__ == "__main__":
    main()
