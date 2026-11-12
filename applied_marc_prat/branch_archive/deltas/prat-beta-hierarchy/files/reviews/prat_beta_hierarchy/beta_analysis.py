"""Task 1: the geometric coefficient beta of Prat (2026) eq 20.

A. beta from the Appendix growth cell (eqs A1-A5) and from alternative cells / habits.
B. Analytic sensitivity of F, C*_top, dh/dt, internal share and the crust limit to beta.
C. Numerical sweep of beta on the reference solver (homogeneous and layered, both variants).

    $PY beta_analysis.py        # ~3 min
"""
from __future__ import annotations

import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (CONSISTENT, FIG, HYBRID, P0, U_EVAP, closure_layer_analytic, d_c_crit, dump, homogeneous_steady,
                    layer_plateaus, layered_run, prat_analytic, replace)

LINES = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); LINES.append(s)


# ============================================================================ A. derivation
def beta_of_cell(area_fn, volume_fn, d=1.0, dd=1e-6):
    """beta = A(d) / (dV/dd), Prat's A3-A5: the layer thickness advances at the rate the crystal
    edge grows, dh/dt = d(d_t)/dt = A k_r (rho/rho_cr) dC / (dV/dd)."""
    dV = (volume_fn(d + dd) - volume_fn(d - dd)) / (2 * dd)
    return area_fn(d) / dV


def section_a():
    say("\n## A. beta from the growth cell")
    cells = {
        "Prat Appendix cell (A1, A2): A = 24 d^2, V = 4 d^3": (lambda d: 24 * d ** 2, lambda d: 4 * d ** 3),
        "one free cube, six faces growing": (lambda d: 6 * d ** 2, lambda d: d ** 3),
        "cube attached by its base (five faces grow)": (lambda d: 5 * d ** 2, lambda d: d ** 3),
        "cube with base + two lateral contacts (three faces grow)": (lambda d: 3 * d ** 2, lambda d: d ** 3),
        "cube with only its top face growing (layer-by-layer overgrowth)": (lambda d: 1 * d ** 2, lambda d: d ** 3),
        "hollow (hopper) cube, wall 0.2 d, inner faces wetted": (lambda d: 6 * d ** 2 + 6 * (0.6 * d) ** 2, lambda d: d ** 3 - (0.6 * d) ** 3),
    }
    rows = []
    for name, (A, V) in cells.items():
        b1 = beta_of_cell(A, V, 1.0); b2 = beta_of_cell(A, V, 7.3)         # dilation invariance
        rows.append((name, b1, b2))
        say(f"  {name:70s} beta = {b1:.4f} (at d = 7.3: {b2:.4f})")
    # a needle d x d x L whose LENGTH sets the layer thickness: only the length grows, A = 4 d L + 2 d^2, dV/dL = d^2
    for L_over_d in (3.0, 10.0):
        say(f"  {'needle d x d x %g d, h along the needle, only L grows' % L_over_d:70s} beta = {4 * L_over_d + 2:.4f} (= 4 L/d + 2)")
        rows.append((f"needle L = {L_over_d} d", 4 * L_over_d + 2, 4 * L_over_d + 2))
    say("  Dilation d -> lambda d leaves beta unchanged for every cell (A ~ d^2, V ~ d^3 => A/(dV/dd) ~ d^0).")
    say("  The gap e (porosity) does not enter A1-A4 at all: beta = 2 is independent of eps_s by construction.")
    say("  Replicating the cell laterally (N cells) multiplies A and dV/dd by N: beta unchanged.")
    say("  What changes beta: faces that cannot grow (contact, partial wetting) -> 2 - n_blocked/3; a habit whose")
    say("  growing area per volume differs from a cube (hollow/hopper cubes of Wijnhorst 2025 micro-CT, needles);")
    say("  face-dependent kinetics (effectively a different k_r per face); a layer-thickness convention other than")
    say("  h = d_t + e (A5).  None of these is a similarity transformation of the cubic cell.")
    # mass-balance view: growing area per unit macroscopic area of a complete crystal layer
    for eps in (0.43, 0.20):
        a_v_dc = 5 * (1 - eps)                       # eq 16: a_v d_c = growing area of one crystal layer per macro area
        say(f"  mass-balance view, eps = {eps}: a_v d_c = {a_v_dc:.2f} growing-area per macro area in one crystal layer,"
            f" i.e. beta_mb = a_v d_c/(1-eps) = 5.00 if the whole top crystal layer grew on all its wetted faces;"
            f" Prat's beta = 2 instead counts cubes growing from small nuclei (left sketch of Figure 3).")
    return rows


# ============================================================================ B. analytic sensitivity
def section_b():
    say("\n## B. analytic sensitivity to beta (eqs 27, 29, 31)")
    betas = np.array([0.5, 1.0, 2 ** 0.5, 5 / 3, 2.0, 2 ** 1.5, 3.0, 4.0, 5.0, 6.0, 8.0])
    out = {}
    for d_sub in (27.5e-6, 112.5e-6, 1080e-6):
        rows = []
        for b in betas:
            p = replace(P0, beta=float(b))
            a = prat_analytic(p.dc_over_dsub * d_sub, 0.43, U_EVAP, p)
            n_star = closure_layer_analytic(d_sub, 0.43, 0.5, p)
            rows.append(dict(beta=float(b), Da_star=a["Da_star"], F=a["F"], v_over_U=a["v"] / U_EVAP, C_top_star=a["C_top_star"],
                             internal_share=1 - a["F"], eps_final_closure=1 - (1 - 0.43) / a["F"], d_c_crit_um=d_c_crit(0.43, p) * 1e6,
                             closure_layer_r05=n_star, elasticity_v=1 - a["F"]))
        out[f"{d_sub*1e6:.1f}"] = rows
        say(f"\n  d_sub = {d_sub*1e6:.1f} um, eps_s = 0.43, r = 0.5")
        say("  beta   Da*      F      v/U     C*_top    1-F    eps_final  d_c,crit[um]  n*  dlnv/dlnbeta")
        for r in rows:
            say(f"  {r['beta']:4.2f}  {r['Da_star']:7.3f}  {r['F']:.4f}  {r['v_over_U']:.4f}  {r['C_top_star']:.3e}  {r['internal_share']:.3f}  "
                f"{r['eps_final_closure']:+.3f}    {r['d_c_crit_um']:8.3f}     {r['closure_layer_r05']}   {r['elasticity_v']:.3f}")
    say("  Scalings: Da* ~ beta^2;  F = sqrt(Da*)/(1+sqrt(Da*));  dlnF/dlnbeta = 1 - F;  C*_top ~ 1/[beta(1-eps)(1-C_eq)k_r + const];")
    say("  d_c,crit ~ beta^-2, hence with r = 0.5 the closure layer shifts by ONE for every factor sqrt(2) in beta:")
    say("  n* = 2 + floor(ln(0.14 d_sub / d_c,crit) / ln 2) = 2 + floor(ln(0.14 d_sub beta^2 / K) / ln 2).")
    dump("beta_analytic", out)
    return out


# ============================================================================ C. numerical sweep
def section_c():
    say("\n## C. numerical sweep (reference solver)")
    t0 = time.time()
    betas = [1.0, 2 ** 0.5, 2.0, 2 ** 1.5, 4.0, 6.0]
    hom = {}
    for d_sub in (27.5e-6, 112.5e-6):
        say(f"\n  homogeneous steady state, d_sub = {d_sub*1e6:.1f} um, eps_new = 0.43")
        say("  beta  | eq29 v/U | hybrid v/U (rel)   C*_top   1-F_share  eps_final | consistent v/U   C*_top   1-F_share  eps_final")
        for b in betas:
            p = replace(P0, beta=float(b)); a = prat_analytic(p.dc_over_dsub * d_sub, 0.43, U_EVAP, p)
            rh = homogeneous_steady(d_sub, 0.43, U_EVAP, p, n_band=40, **HYBRID)
            try:
                rc = homogeneous_steady(d_sub, 0.43, U_EVAP, p, n_band=40, **CONSISTENT); cc = dict(v=rc["v"], C=rc["C_top_star"], s=1 - rc["F_num"], e=rc["eps_final"])
            except RuntimeError as ex:
                cc = None
            hom[f"{d_sub*1e6:.1f}_{b:.3f}"] = dict(beta=b, analytic=a, hybrid={k: rh[k] for k in ("v", "C_top_star", "F_num", "F_prat", "eps_final")},
                                                   consistent=None if cc is None else cc)
            say(f"  {b:4.2f}  | {a['v']/U_EVAP:.4f}   | {rh['v']/U_EVAP:.4f} ({rh['v']/a['v']-1:+.1e})  {rh['C_top_star']:.3e}  {1-rh['F_num']:.3f}  {rh['eps_final']:+.4f}   | "
                + ("closed (band porosity -> 0)" if cc is None else f"{cc['v']/U_EVAP:.4f}   {cc['C']:.3e}  {cc['s']:.3f}  {cc['e']:+.4f}"))
    dump("beta_numeric_homogeneous", hom)
    lay = {}
    say("\n  layered 27.5 um (r 0.5, alpha 2, c_h 1, eps_new 0.43, 5 layers): closure layer / h / t  vs analytic n*")
    say("  beta  | analytic n* | hybrid                     | consistent")
    for b in betas:
        p = replace(P0, beta=float(b))
        row = dict(beta=b, n_star=closure_layer_analytic(27.5e-6, 0.43, 0.5, p))
        for lab, var in (("hybrid", HYBRID), ("consistent", CONSISTENT)):
            run = layered_run(27.5e-6, 0.43, p=p, variant=var, n_layers=5, record_every=5)
            cl = run["closure"]
            row[lab] = dict(closure=cl, plateaus=layer_plateaus(run), t_end=float(run["t"][-1]), h_end=float(run["h"][-1]))
        lay[f"{b:.3f}"] = row
        f = lambda c: "none in 5 layers" if c is None else f"layer {c['layer']}, h {c['h']*1e6:6.1f} um, t {c['t']/3600:5.2f} h"
        say(f"  {b:4.2f}  | {row['n_star']:^11} | {f(row['hybrid']['closure']):28s} | {f(row['consistent']['closure'])}")
    dump("beta_numeric_layered", lay)
    # figure
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    bb = np.linspace(0.5, 8, 200)
    for d_sub, c in ((27.5e-6, "r"), (112.5e-6, "b"), (1080e-6, "k")):
        A = [prat_analytic(P0.dc_over_dsub * d_sub, 0.43, U_EVAP, replace(P0, beta=float(b))) for b in bb]
        ax[0].plot(bb, [a["v"] / U_EVAP for a in A], c, label=f"{d_sub*1e6:.0f} um"); ax[1].plot(bb, [a["C_top_star"] for a in A], c)
        ax[2].plot(bb, [1 - a["F"] for a in A], c)
    for d_sub, c in ((27.5e-6, "r"), (112.5e-6, "b")):
        for lab, m in (("hybrid", "o"), ("consistent", "s")):
            pts = [(hom[f"{d_sub*1e6:.1f}_{b:.3f}"]["beta"], hom[f"{d_sub*1e6:.1f}_{b:.3f}"][lab]) for b in betas]
            pts = [(b, q) for b, q in pts if q is not None]
            ax[0].plot([b for b, _ in pts], [q["v"] / U_EVAP for _, q in pts], c + m, mfc="none" if lab == "hybrid" else c, ms=5)
            ax[1].plot([b for b, _ in pts], [q["C_top_star" if lab == "hybrid" else "C"] for _, q in pts], c + m, mfc="none" if lab == "hybrid" else c, ms=5)
            ax[2].plot([b for b, _ in pts], [(1 - q["F_num"]) if lab == "hybrid" else q["s"] for _, q in pts], c + m, mfc="none" if lab == "hybrid" else c, ms=5)
    for a_, yl in zip(ax, ("v / U_evap (eq 29 lines; o hybrid, s consistent)", "C*_top", "internal precipitation share 1 - F")):
        a_.set_xlabel("beta"); a_.set_ylabel(yl); a_.axvline(2.0, color="0.7", lw=0.8)
    ax[0].legend(fontsize=7); ax[1].set_yscale("log")
    fig.suptitle("beta sensitivity, eps_s = 0.43, U_evap = 1 cm/day"); fig.tight_layout(); fig.savefig(FIG / "beta_sweep.png", dpi=120); plt.close(fig)
    say(f"  [C {time.time()-t0:.0f} s]")


if __name__ == "__main__":
    section_a(); section_b(); section_c()
    (FIG.parent / "results" / "beta_analysis.md").write_text("# beta_analysis.py output\n" + "\n".join(LINES) + "\n")
