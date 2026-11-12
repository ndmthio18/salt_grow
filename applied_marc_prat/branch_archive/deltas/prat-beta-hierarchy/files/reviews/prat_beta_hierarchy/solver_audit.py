"""Task 2b: audit of the moving-frame finite-volume solver at discontinuous layer properties, and
qualification of the eps_final infeasibility claim.

A. layer 2 -> 3 transition of the 27.5 um run under grid, time-step, stretching and face-averaging changes
B. width of the recorded eps_final(z) jump across a layer boundary (numerical diffusion of the eps advection)
C. dense birth-porosity scan for the eps_final interpretation on fine crystals, two grids, both variants

    $PY solver_audit.py        # ~4 min
"""
from __future__ import annotations

import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (CONSISTENT, FIG, HYBRID, P0, U_EVAP, dump, homogeneous_steady, layer_plateaus, layered_run,
                    prat_analytic, transition_duration)

LINES = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); LINES.append(s)


def section_a():
    say("\n## A. layer 2 -> 3 transition and closure under discretisation changes (27.5 um, eps_new 0.43, consistent)")
    t0 = time.time(); out = {}; runs = {}
    settings = [("n_band=20", dict(n_band=20)), ("n_band=40 (reference)", dict(n_band=40)), ("n_band=80", dict(n_band=80)),
                ("n_band=160", dict(n_band=160)), ("dt x0.2", dict(dt_over_band=0.1)), ("dt x4", dict(dt_over_band=2.0)),
                ("stretch 1.015", dict(stretch=1.015)), ("L = 120 delta_C", dict(L_over=120.0)),
                ("harmonic face average", dict(face_average="harmonic")), ("hybrid variant", dict(variant=HYBRID)),
                ("dt x4, converged check", dict(dt_over_band=2.0, closure_check="converged")),
                ("dt x8, converged check", dict(dt_over_band=4.0, closure_check="converged")),
                ("dt x8, iterate check", dict(dt_over_band=4.0))]
    say("  setting                 N    F2        F3        tau_23 [dC/v]  tau_34   closure layer  h [um]    t [h]    max Picard  cum res salt/solid")
    for lab, kw in settings:
        run = layered_run(27.5e-6, 0.43, n_layers=5, record_every=1, **kw)
        pl = layer_plateaus(run); tr3 = transition_duration(run, 3); tr4 = transition_duration(run, 4); cl = run["closure"]
        out[lab] = dict(N=run["N"], plateaus=pl, tr3=tr3, tr4=tr4, closure=cl, picard_max=int(run["picard"].max()), cum=run["cumulative"])
        runs[lab] = run
        say(f"  {lab:22s} {run['N']:4d}  {pl[2]['F']:.5f}  {pl[3]['F']:.5f}  {tr3['tau_over_band'] if tr3 else np.nan:7.2f}        "
            f"{tr4['tau_over_band'] if tr4 else np.nan:6.2f}   {cl['layer'] if cl else '-':^13}  {cl['h']*1e6 if cl else np.nan:8.2f}  {cl['t']/3600 if cl else np.nan:6.3f}   "
            f"{int(run['picard'].max()):3d}        {run['cumulative']['res_salt']:.1e} / {run['cumulative']['res_solid']:.1e}")
    say("  Reading: the plateau F of a layer is grid-, dt- and averaging-independent to the digits shown and the transient after")
    say("  a boundary lasts a few delta_C/v (seconds to minutes) against layer transit times of hours.  DEFECT of the original")
    say("  solver: closure is tested on every Picard iterate, and the first iterate uses the lagged sink, so a step of")
    say("  >= 2 delta_C/v overshoots and closes one layer EARLY (layer 4 instead of 5).  Testing the converged step (extended copy,")
    say("  closure_check='converged') makes the closure layer dt-independent up to 4 delta_C/v.  All results of the previous")
    say("  review used dt = 0.5 delta_C/v, where both checks agree; with the converged check the closure height moves by")
    say("  < one band thickness per factor 4 in dt.")
    dump("audit_transition", out)
    # zoom figure
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for lab, c in (("n_band=20", "0.6"), ("n_band=40 (reference)", "r"), ("n_band=160", "b"), ("harmonic face average", "g"), ("dt x4", "m")):
        run = runs[lab]; t_b = run["t"][np.argmax(run["layer_top"] == 3)]
        sel = (run["t"] > t_b - 400) & (run["t"] < t_b + 1500)
        ax[0].plot(run["t"][sel] - t_b, run["v"][sel] / U_EVAP, color=c, label=lab, lw=1)
    ax[0].set_xlabel("t - t(layer 3 reaches the interface) [s]"); ax[0].set_ylabel("v / U_evap"); ax[0].legend(fontsize=7); ax[0].set_title("layer 2 -> 3 transient")
    for lab, c in (("n_band=20", "0.6"), ("n_band=40 (reference)", "r"), ("stretch 1.015", "b"), ("n_band=160", "k")):
        run = runs[lab]; ok = run["z_m_final"] > 0
        ax[1].plot(run["z_m_final"][ok] * 1e6, run["eps_final"][ok], color=c, label=lab, lw=1)
    for zb in runs["n_band=40 (reference)"]["layer_bottoms"][1:4]:
        ax[1].axvline(zb * 1e6, color="0.8", lw=0.6)
    ax[1].set_xlabel("z [um] (material)"); ax[1].set_ylabel("recorded eps_final"); ax[1].legend(fontsize=7); ax[1].set_title("eps_final(z): jump smearing")
    fig.tight_layout(); fig.savefig(FIG / "audit_transition.png", dpi=120); plt.close(fig)
    say(f"  [A {time.time()-t0:.0f} s]")
    return runs


def section_b(runs):
    say("\n## B. width of the recorded eps_final(z) jump at the layer 1 -> 2 boundary (10 %-90 % of the jump)")
    out = {}
    for lab in ("n_band=20", "n_band=40 (reference)", "n_band=80", "n_band=160", "stretch 1.015", "L = 120 delta_C"):
        run = runs[lab]; z = run["z_m_final"]; e = run["eps_final"]; zb = run["layer_bottoms"][1]
        sel = (z > zb - 15e-6) & (z < zb + 15e-6)
        if sel.sum() < 5:
            continue
        zz, ee = z[sel], e[sel]
        e_hi, e_lo = np.median(ee[zz < zb - 8e-6]), np.median(ee[zz > zb + 8e-6])
        lvl = lambda f: e_hi - f * (e_hi - e_lo)
        z10 = zz[np.argmax(ee <= lvl(0.1))]; z90 = zz[np.argmax(ee <= lvl(0.9))]
        # cell size at the window bottom (where eps_final is read) and physical band thickness
        out[lab] = dict(width_um=(z90 - z10) * 1e6, e_hi=e_hi, e_lo=e_lo, L_um=run["L"] * 1e6, delta_C1_um=run["delta_C"][0] * 1e6)
        say(f"  {lab:22s} jump {e_hi:.4f} -> {e_lo:.4f}, 10-90 % width {(z90 - z10)*1e6:6.2f} um  (delta_C,1 = {run['delta_C'][0]*1e6:.3f} um, window L = {run['L']*1e6:.1f} um)")
    say("  Reading: the physical transition is ~delta_C (the band); the recorded width is set by the coarse cells near the")
    say("  window bottom (implicit upwind advection of eps over L): it shrinks with the stretching, not with n_band, and it")
    say("  does not affect per-layer medians or closure.  Figures s2_* of the previous review show this smeared width.")
    dump("audit_eps_jump_width", out)


def section_c():
    say("\n## C. dense birth-porosity scan for the eps_final interpretation (target 0.43), both variants")
    t0 = time.time(); out = {}
    one_minus = np.logspace(np.log10(0.57), np.log10(1e-3), 25)           # 1 - eps_new from 0.57 down to 0.001
    eps_news = 1 - one_minus
    for d_c in (0.481e-6, 0.241e-6, 0.120e-6):
        d_sub = d_c / P0.dc_over_dsub
        for lab, var in (("hybrid", HYBRID), ("consistent n_band 40", CONSISTENT), ("consistent n_band 80", dict(CONSISTENT, n_band=80))):
            kw = dict(var); nb = kw.pop("n_band", 40)
            ef = []
            for e in eps_news:
                try:
                    r = homogeneous_steady(d_sub, float(e), U_EVAP, P0, n_band=nb, **kw); ef.append(r["eps_final"])
                except RuntimeError:
                    ef.append(0.0)
            ef = np.array(ef); feasible = ef >= 0.43
            out[f"{d_c*1e6:.3f}_{lab}"] = dict(eps_new=eps_news, eps_final=ef, max_eps_final=float(ef.max()), argmax_eps_new=float(eps_news[ef.argmax()]),
                                               feasible_range=[float(eps_news[feasible].min()), float(eps_news[feasible].max())] if feasible.any() else None)
            say(f"  d_c = {d_c*1e6:.3f} um {lab:22s}: max eps_final {ef.max():.4f} at eps_new {eps_news[ef.argmax()]:.4f}; "
                + (f"eps_final >= 0.43 for eps_new in [{eps_news[feasible].min():.3f}, {eps_news[feasible].max():.3f}]" if feasible.any() else "NO eps_new in [0.43, 0.999] reaches 0.43")
                + f"; closed (eps -> 0) at {int((ef == 0).sum())} of 25 points")
    say("  Qualification: the scan covers 25 points in 1 - eps_new from 0.57 to 1e-3 on two grids.  Up to eps_new = 0.99 the")
    say("  consistent model reaches eps_final = 0.43 on neither grid for d_c <= 0.24 um, while the linearised model does")
    say("  (eps_new ~ 0.93).  Beyond 0.99 a grid-DEPENDENT sliver appears (eps_new >= 0.994-0.998) where a_v = 5(1-eps)/d_c -> 0")
    say("  makes the band inert and eps_final ~ eps_new: Prat's packing closure has no meaning at 99 % void, and the numerics")
    say("  are not converged there.  The previous review's 'infeasible' claim is therefore restated as: no birth porosity")
    say("  <= 0.99 reaches the target; it rests on the local-eps closure and is a finite-scan statement, not a theorem.")
    dump("audit_infeasibility_scan", out)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for k, d_c in enumerate((0.481e-6, 0.241e-6, 0.120e-6)):
        for lab, c in (("hybrid", "k"), ("consistent n_band 40", "r"), ("consistent n_band 80", "b")):
            q = out[f"{d_c*1e6:.3f}_{lab}"]; ax[k].plot(q["eps_new"], q["eps_final"], c + ".-", label=lab, lw=1)
        ax[k].axhline(0.43, color="0.6", ls=":"); ax[k].set_xlabel("eps_new"); ax[k].set_ylabel("eps_final"); ax[k].set_title(f"d_c = {d_c*1e6:.2f} um"); ax[k].legend(fontsize=7)
    fig.suptitle("eps_final(eps_new): where the target 0.43 is reachable (0 = band closed)"); fig.tight_layout()
    fig.savefig(FIG / "audit_infeasibility.png", dpi=120); plt.close(fig)
    say(f"  [C {time.time()-t0:.0f} s]")


if __name__ == "__main__":
    runs = section_a(); section_b(runs); section_c()
    (FIG.parent / "results" / "solver_audit.md").write_text("# solver_audit.py output\n" + "\n".join(LINES) + "\n")
