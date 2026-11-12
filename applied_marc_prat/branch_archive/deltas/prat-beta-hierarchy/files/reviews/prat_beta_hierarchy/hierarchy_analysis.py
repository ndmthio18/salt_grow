"""Task 2a: the imposed hierarchy (r, alpha, c_h): what the reaction-transport problem decides and what
the prescribed layer geometry decides.

A. analytic maps: closure layer n*(r, d_sub), trigger height h*(alpha, c_h), the exponent of h* ~ d_sub^m.
B. numerical sweeps on the reference solver (27.5 um, eps_new 0.43): r, alpha, c_h varied one at a time;
   closure layer / height / time, per-layer plateaus, transition transients.

    $PY hierarchy_analysis.py       # ~4 min
"""
from __future__ import annotations

import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (CONSISTENT, FIG, HYBRID, P0, U_EVAP, closure_layer_analytic, d_c_crit, dump, layer_bottom,
                    layer_plateaus, layered_run, prat_analytic, transition_duration)

LINES = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s); LINES.append(s)


def section_a():
    say("\n## A. analytic maps")
    say(f"  d_c,crit (beta 2, Prat params): eps 0.43 -> {d_c_crit(0.43)*1e6:.3f} um; eps 0.20 -> {d_c_crit(0.20)*1e6:.3f} um")
    say("  Reaction-transport decides n* through d_c(n) < d_c,crit only; alpha and c_h never enter.")
    rs = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]); ds = np.array([5, 10, 27.5, 56, 112.5, 250, 500, 1080]) * 1e-6
    maps = {}
    for eps in (0.43, 0.20):
        M = np.array([[closure_layer_analytic(d, eps, r) for d in ds] for r in rs], float)
        maps[f"n_star_eps{eps}"] = M
        say(f"\n  closure layer n*, eps_s = {eps} (rows r = {rs.tolist()}, cols d_sub[um] = {(ds*1e6).tolist()})")
        for r, row in zip(rs, M):
            say(f"  r={r:.1f}: " + "  ".join(f"{int(x):2d}" if np.isfinite(x) else " -" for x in row))
    # trigger height and its exponent
    alphas = np.array([1.25, 1.5, 2.0, 3.0, 4.0]); chs = np.array([0.3, 0.5, 1.0, 1.4, 2.0, 3.0])
    n_star = closure_layer_analytic(27.5e-6, 0.43, 0.5)
    H = np.array([[layer_bottom(n_star, 27.5e-6, a, c) for c in chs] for a in alphas])
    maps["h_star_27um"] = H
    say(f"\n  trigger height h* = layer bottom of n* = {n_star} for 27.5 um, eps 0.43, r 0.5 [um] (rows alpha, cols c_h):")
    for a, row in zip(alphas, H):
        say(f"  alpha={a:.2f}: " + "  ".join(f"{x*1e6:7.1f}" for x in row))
    say("  h* = c_h d_sub (alpha^(n*-1) - 1)/(alpha - 1): linear in c_h, geometric in alpha; the layer index n* is untouched.")
    say("  Exponent of h* ~ d_sub^m with n*(d_sub) from reaction-transport: m = 1 + ln(alpha)/ln(1/r) (plan/13), e.g.")
    for r in (0.4, 0.5, 0.6, 0.7):
        say("   r = %.1f: " % r + ", ".join(f"alpha {a:.2f}: m = {1 + np.log(a)/np.log(1/r):.2f}" for a in alphas))
    say("  m is pure geometry (alpha, r); only the prefactor (through d_c,crit ~ D_eff(1-eps)/(beta^2 eps (1-C_eq) k_r)) is physics.")
    dump("hierarchy_analytic_maps", dict(rs=rs, ds=ds, alphas=alphas, chs=chs, **maps))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    im = ax[0].imshow(maps["n_star_eps0.43"], origin="lower", aspect="auto", cmap="viridis")
    ax[0].set_xticks(range(len(ds))); ax[0].set_xticklabels([f"{d*1e6:g}" for d in ds], fontsize=7); ax[0].set_yticks(range(len(rs))); ax[0].set_yticklabels(rs)
    ax[0].set_xlabel("d_sub [um]"); ax[0].set_ylabel("r"); ax[0].set_title("closure layer n* (eps_s 0.43, beta 2)"); plt.colorbar(im, ax=ax[0])
    for (i, j), val in np.ndenumerate(maps["n_star_eps0.43"]):
        ax[0].text(j, i, f"{int(val)}" if np.isfinite(val) else "-", ha="center", va="center", color="w", fontsize=7)
    im2 = ax[1].imshow(np.log10(H * 1e6), origin="lower", aspect="auto", cmap="magma")
    ax[1].set_xticks(range(len(chs))); ax[1].set_xticklabels(chs); ax[1].set_yticks(range(len(alphas))); ax[1].set_yticklabels(alphas)
    ax[1].set_xlabel("c_h"); ax[1].set_ylabel("alpha"); ax[1].set_title("log10 h* [um] at n* = %d (27.5 um)" % n_star); plt.colorbar(im2, ax=ax[1])
    fig.tight_layout(); fig.savefig(FIG / "hierarchy_maps.png", dpi=120); plt.close(fig)


def section_b():
    say("\n## B. numerical sweeps (27.5 um, eps_new 0.43, consistent unless stated)")
    t0 = time.time(); out = {}
    def report(tag, run, n_star):
        cl = run["closure"]; pl = layer_plateaus(run)
        tr = {n: transition_duration(run, n) for n in range(2, run["n_layers"] + 1)}
        trs = ", ".join(f"{n}: {t['tau_over_band']:.1f}" for n, t in tr.items() if t)
        say(f"  {tag:34s} closure {'none' if cl is None else f'layer {cl['layer']} h {cl['h']*1e6:7.1f} um t {cl['t']/3600:5.2f} h'} | analytic n* {n_star} | "
            f"F plateaus " + " ".join(f"{pl[n]['F']:.3f}" for n in sorted(pl)) + f" | transient [delta_C/v]: {trs}")
        out[tag] = dict(closure=cl, plateaus=pl, transients=tr, n_star=n_star, N=run["N"], cum=run["cumulative"])
    say("\n  r sweep (alpha 2, c_h 1, up to 7 layers or 1.5 mm)")
    for r in (0.4, 0.5, 0.6, 0.7, 0.8):
        for lab, var in (("hybrid", HYBRID), ("consistent", CONSISTENT)):
            run = layered_run(27.5e-6, 0.43, r=r, n_layers=7, variant=var, record_every=5, h_max=1.5e-3)
            report(f"r={r} {lab}", run, closure_layer_analytic(27.5e-6, 0.43, r))
    say("\n  alpha sweep (r 0.5, c_h 1, 5 layers)")
    for a in (1.25, 1.5, 2.0, 3.0, 4.0):
        run = layered_run(27.5e-6, 0.43, alpha=a, n_layers=5, record_every=5)
        report(f"alpha={a} consistent", run, closure_layer_analytic(27.5e-6, 0.43, 0.5))
    say("\n  c_h sweep (r 0.5, alpha 2, 5 layers)")
    for c in (0.3, 0.5, 1.0, 1.4, 2.0, 3.0):
        run = layered_run(27.5e-6, 0.43, c_h=c, n_layers=5, record_every=5)
        report(f"c_h={c} consistent", run, closure_layer_analytic(27.5e-6, 0.43, 0.5))
    dump("hierarchy_numeric_sweeps", out)
    # figure: closure h, t vs alpha and c_h; closure layer vs r
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    rr = [0.4, 0.5, 0.6, 0.7, 0.8]
    for lab, m in (("hybrid", "o"), ("consistent", "s")):
        ax[0].plot(rr, [out[f"r={r} {lab}"]["closure"]["layer"] if out[f"r={r} {lab}"]["closure"] else np.nan for r in rr], m + "-", label=lab)
    ax[0].plot(rr, [out[f"r={r} hybrid"]["n_star"] for r in rr], "k:", label="analytic n*")
    ax[0].set_xlabel("r"); ax[0].set_ylabel("closure layer"); ax[0].legend(fontsize=7)
    aa = [1.25, 1.5, 2.0, 3.0, 4.0]
    ax[1].plot(aa, [out[f"alpha={a} consistent"]["closure"]["h"] * 1e6 for a in aa], "s-", label="numerical")
    ax[1].plot(aa, [layer_bottom(5, 27.5e-6, a, 1.0) * 1e6 for a in aa], "k:", label="layer-5 bottom")
    ax[1].set_xlabel("alpha"); ax[1].set_ylabel("closure height [um]"); ax[1].set_yscale("log"); ax[1].legend(fontsize=7)
    cc = [0.3, 0.5, 1.0, 1.4, 2.0, 3.0]
    ax[2].plot(cc, [out[f"c_h={c} consistent"]["closure"]["h"] * 1e6 for c in cc], "s-", label="numerical")
    ax[2].plot(cc, [layer_bottom(5, 27.5e-6, 2.0, c) * 1e6 for c in cc], "k:", label="layer-5 bottom")
    ax[2].set_xlabel("c_h"); ax[2].set_ylabel("closure height [um]"); ax[2].legend(fontsize=7)
    fig.suptitle("27.5 um, eps_new 0.43: the closure LAYER follows r; alpha and c_h only place it"); fig.tight_layout()
    fig.savefig(FIG / "hierarchy_sweeps.png", dpi=120); plt.close(fig)
    say(f"  [B {time.time()-t0:.0f} s]")


if __name__ == "__main__":
    section_a(); section_b()
    (FIG.parent / "results" / "hierarchy_analysis.md").write_text("# hierarchy_analysis.py output\n" + "\n".join(LINES) + "\n")
