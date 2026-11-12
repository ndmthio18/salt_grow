"""Runs every comparison of the validation and writes results/ and figures/.

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
    cd reviews/prat_layered_validation && $PY run_validation.py            # ~3-6 min

Sections (each writes its own json/csv and prints a markdown table into results/summary.md):
  S1  homogeneous verification: numerical reference vs Prat eqs 19, 29, 31 and vs the project's
      birth-porosity closure, for the four model variants (V0 hybrid, V2 constant eps, V3 no liquid
      consumption, P consistent) and spatial / temporal convergence
  S2  layered structures: (A) homogeneous Prat with layer-1 properties, (B) layer-wise Prat, (C) the
      numerical reference; eps_new and eps_final interpretations; eps_s = 0.43 and 0.20
  S3  layer-thickness sensitivity (c_h, alpha) and window / grid convergence of the layered run
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from layered_reference import (DAY, Layering, PratParams, ReferenceSolver, eps_new_for_final,
                               homogeneous_steady, prat_analytic, prat_profile)

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; FIG = HERE / "figures"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
P = PratParams()
U_EVAP = 1e-2 / DAY                      # 1 cm/day, the paper's example (Figure 6 text)
SUMMARY = []


def md(title, header, rows):
    s = [f"\n### {title}\n", "| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        s.append("| " + " | ".join(str(x) for x in r) + " |")
    SUMMARY.extend(s); print("\n".join(s))


def dump(name, obj):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        return o
    (RES / f"{name}.json").write_text(json.dumps(conv(obj), indent=1))


# ============================================================================ S1 homogeneous
def s1_homogeneous():
    t0 = time.time()
    variants = {
        "V0 hybrid: eps evolves (2)-(3), transport at birth eps (Prat's linearisation)": dict(evolve_eps=True, local_props=False, transport_eps_fixed=True),
        "V2 constant eps everywhere (self-replication as literally stated)": dict(evolve_eps=False, local_props=False),
        "V3 as V0 but the liquid consumption of (3)-(4) dropped": dict(evolve_eps=True, local_props=False, transport_eps_fixed=True, liquid_consumption=False),
        "P  consistent: eps evolves in every term, a_v and D_eff local": dict(evolve_eps=True, local_props=True),
    }
    rows = []; out = {}
    for d in (10e-6, 25e-6, 56e-6, 112.5e-6, 1080e-6):
        a = prat_analytic(P.dc_over_dsub * d, 0.43, U_EVAP, P)
        eps_closure = 1 - (1 - 0.43) / a["F"]
        rows.append([f"{d*1e6:.1f}", "Prat eq 29/31 (A)", f"{a['v']/U_EVAP:.5f}", "—", f"{a['F']:.4f}", f"{a['C_top_star']:.4e}", f"{eps_closure:.4f}", "—", "—"])
        out[f"{d*1e6:.1f}"] = dict(analytic=a, eps_closure=eps_closure, variants={})
        for lab, kw in variants.items():
            r = homogeneous_steady(d, 0.43, U_EVAP, P, n_band=40, **kw)
            out[f"{d*1e6:.1f}"]["variants"][lab] = {k: r[k] for k in ("v", "F_prat", "F_num", "C_top_star", "eps_final", "res_salt", "res_solid", "N", "W_i", "W_L")}
            out[f"{d*1e6:.1f}"]["variants"][lab].update(xi=r["xi"], C=r["C"], eps=r["eps"])
            rows.append([f"{d*1e6:.1f}", lab.split(":")[0], f"{r['v']/U_EVAP:.5f}", f"{r['v']/a['v']-1:+.2e}", f"{r['F_prat']:.4f}",
                         f"{r['C_top_star']:.4e}", f"{r['eps_final']:.4f}", f"{r['res_salt']:.1e}", f"{r['res_solid']:.1e}"])
    md("S1a Homogeneous steady state, eps_new = 0.43, U_evap = 1 cm/day, Prat parameters (n_band = 40)",
       ["d_sub [um]", "model", "v/U_evap", "rel. to eq 29", "F (Prat def.)", "C*_top", "eps_final", "res salt", "res solid"], rows)
    dump("s1a_homogeneous", out)

    # profiles vs eq 19 (10 um): V0 and P
    d = 10e-6; a = prat_analytic(P.dc_over_dsub * d, 0.43, U_EVAP, P)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for lab, kw in (("V0 hybrid", dict(evolve_eps=True, local_props=False, transport_eps_fixed=True)),
                    ("P consistent", dict(evolve_eps=True, local_props=True))):
        r = homogeneous_steady(d, 0.43, U_EVAP, P, n_band=80, **kw)
        ax[0].plot(r["xi"] * 1e6, (r["C"] - P.C_eq) / a["C_top_star"], label=f"{lab} (numerical)")
        ax[1].plot(r["xi"] * 1e6, r["eps"], label=lab)
        err = np.max(np.abs((r["C"] - P.C_eq) - prat_profile(r["xi"], a["delta_C"], a["C_top_star"]))) / a["C_top_star"]
        SUMMARY.append(f"\nprofile {lab}: max |C*_num - eq19| / C*_top(eq31) = {err:.3e} over the window; eps_final = {r['eps_final']:.4f}")
    xx = np.linspace(0, 6 * a["delta_C"], 300)
    ax[0].plot(xx * 1e6, prat_profile(xx, a["delta_C"], a["C_top_star"]) / a["C_top_star"], "k--", label="eq 19 (exp, Da>>1)")
    ax[0].set_xlim(0, 6 * a["delta_C"] * 1e6); ax[0].set_xlabel("xi = h - z [um]"); ax[0].set_ylabel("C* / C*_top(eq 31)"); ax[0].legend(fontsize=7)
    ax[1].set_xlim(0, 6 * a["delta_C"] * 1e6); ax[1].set_xlabel("xi [um]"); ax[1].set_ylabel("porosity"); ax[1].legend(fontsize=7)
    ax[1].axhline(1 - 0.57 / a["F"], color="k", ls=":", label="closure 1-(1-eps)/F")
    fig.suptitle("d_sub = 10 um, eps_new = 0.43: band profiles"); fig.tight_layout(); fig.savefig(FIG / "s1_profiles_10um.png", dpi=130); plt.close(fig)

    # convergence: grid, window, time step (10 um, V0 and P)
    rows = []; conv = {}
    for lab, kw in (("V0", dict(evolve_eps=True, local_props=False, transport_eps_fixed=True)), ("P", dict(evolve_eps=True, local_props=True))):
        for nb in (10, 20, 40, 80, 160):
            r = homogeneous_steady(d, 0.43, U_EVAP, P, n_band=nb, **kw)
            conv[f"{lab}_nband{nb}"] = dict(v=r["v"], eps_final=r["eps_final"], N=r["N"], C_top_star=r["C_top_star"])
            rows.append([lab, f"n_band={nb}", r["N"], f"{r['v']/U_EVAP:.6f}", f"{r['v']/a['v']-1:+.2e}", f"{r['eps_final']:.6f}"])
        for Lf in (30.0, 60.0, 120.0):
            r = homogeneous_steady(d, 0.43, U_EVAP, P, n_band=40, L_over_deltaC=Lf, **kw)
            conv[f"{lab}_L{Lf:.0f}"] = dict(v=r["v"], eps_final=r["eps_final"], N=r["N"])
            rows.append([lab, f"L={Lf:.0f} delta_C", r["N"], f"{r['v']/U_EVAP:.6f}", f"{r['v']/a['v']-1:+.2e}", f"{r['eps_final']:.6f}"])
        for dtf in (0.1, 0.5, 2.0):
            r = homogeneous_steady(d, 0.43, U_EVAP, P, n_band=40, dt_over_transit=dtf, **kw)
            conv[f"{lab}_dt{dtf}"] = dict(v=r["v"], eps_final=r["eps_final"])
            rows.append([lab, f"dt={dtf} delta_C/v", r["N"], f"{r['v']/U_EVAP:.6f}", f"{r['v']/a['v']-1:+.2e}", f"{r['eps_final']:.6f}"])
    md("S1b Convergence of the homogeneous steady state (d_sub = 10 um)", ["model", "setting", "N cells", "v/U_evap", "rel. to eq 29", "eps_final"], rows)
    dump("s1b_convergence", conv)
    print(f"[S1 {time.time()-t0:.0f} s]")


# ============================================================================ S2 layered
def layered_case(d_sub, eps_s, interpretation, r=0.5, alpha=2.0, c_h=1.0, n_layers=5, n_band=40, L_over=60.0,
                 stretch=1.03, record_every=5, dt_over_band=0.5, variant=None):
    """One layered run.  interpretation: 'eps_new' (every layer born at eps_s) or 'eps_final' (birth
    porosity of each layer chosen so that its steady final porosity is eps_s).  Returns (C) with the
    analytic (A) and (B) velocities evaluated along its own h(t)."""
    variant = dict(evolve_eps=True, local_props=True) if variant is None else variant
    lay = Layering(d_sub=d_sub, r=r, alpha=alpha, c_h=c_h)
    # birth porosity table
    eps_new = {}; shooting = {}; infeasible_layer = None
    for n in range(1, n_layers + 2):
        if interpretation == "eps_new":
            eps_new[n] = eps_s
        else:
            e, hist, rr = eps_new_for_final(d_sub * r ** (n - 1), eps_s, U_EVAP, P, n_band=n_band, **variant)
            shooting[n] = hist
            if e is None:                                   # no birth porosity reaches eps_final = eps_s on this layer
                infeasible_layer = n
                n_layers = n - 1
                break
            eps_new[n] = e
    if n_layers == 0:
        return dict(d_sub=d_sub, eps_s=eps_s, interpretation=interpretation, r=r, alpha=alpha, c_h=c_h,
                    n_layers=0, infeasible_layer=infeasible_layer, eps_new=eps_new, shooting=shooting, closure=None)
    # window and first cell from the finest layer that will be reached
    a1 = prat_analytic(lay.d_c(1, P), eps_new[1], U_EVAP, P)
    an = prat_analytic(lay.d_c(n_layers, P), eps_new[n_layers], U_EVAP, P)
    s = ReferenceSolver(lay, U_EVAP, eps_new, P, L=L_over * a1["delta_C"], dxi0=an["delta_C"] / n_band, stretch=stretch, **variant)
    h_end = lay.layer_bottom(n_layers + 1)
    dt = dt_over_band * an["delta_C"] / max(an["v"], 1e-30)
    t_end = 5 * h_end / max(an["v"] * 0.3, 1e-30)
    s.run(t_end, dt, h_end=h_end, record_every=record_every)
    H = s.history
    t = np.array([q["t"] for q in H]); h = np.array([q["h"] for q in H]); v = np.array([q["v"] for q in H])
    n_top = np.array([q["layer_top"] for q in H])
    # (A): layer-1 Prat everywhere; (B): Prat of the layer at the interface
    vA = np.full_like(h, a1["v"]); vB = np.empty_like(h); FB = np.empty_like(h); eps_fin_B = np.empty_like(h)
    for k, n in enumerate(n_top):
        if interpretation == "eps_new":
            b = prat_analytic(lay.d_c(n, P), eps_s, U_EVAP, P)
            vB[k] = b["v"]; FB[k] = b["F"]; eps_fin_B[k] = 1 - (1 - eps_s) / b["F"]
        else:
            b = prat_analytic(lay.d_c(n, P), eps_s, U_EVAP, P)          # Da* at the measured (final) porosity
            en = 1 - b["F"] * (1 - eps_s)                                # birth porosity implied by F
            vB[k] = P.rho * P.C_eq * U_EVAP / (P.rho_cr * (1 - en) * (1 - P.C_eq)) * b["F"]   # = v_max(eps_s): F cancels
            FB[k] = b["F"]; eps_fin_B[k] = eps_s
    hA = np.cumsum(np.concatenate([[0], np.diff(t)]) * vA); hB = np.cumsum(np.concatenate([[0], np.diff(t)]) * vB)
    # material that left the window carries its final porosity: z_m = h - L at that time
    out = dict(d_sub=d_sub, eps_s=eps_s, interpretation=interpretation, r=r, alpha=alpha, c_h=c_h, n_layers=n_layers,
               n_band=n_band, L_over=L_over, N=int(s.xi.size), dt=dt, eps_new=eps_new, shooting=shooting,
               t=t, h=h, v=v, vA=vA, vB=vB, hA=hA, hB=hB, FB=FB, eps_fin_B=eps_fin_B,
               F_prat=np.array([q["F_prat"] for q in H]), F_share=np.array([q["F_num"] for q in H]),
               C_top_star=np.array([q["C_top_star"] for q in H]), S_int=np.array([q["S_int"] for q in H]),
               P_i=np.array([q["P_i"] for q in H]), eps_final=np.array([q["eps_final"] for q in H]),
               z_m_final=h - s.L, layer_top=n_top, res_salt=np.array([q["res_salt"] for q in H]),
               res_solid=np.array([q["res_solid"] for q in H]), picard=np.array([q["picard"] for q in H]),
               closure=s.closure_event, infeasible_layer=infeasible_layer, cumulative=s.cumulative_balance(),
               layer_bottoms=[lay.layer_bottom(n) for n in range(1, n_layers + 2)],
               final_profile=dict(xi=s.xi, C=s.C, eps=s.eps))
    return out


def s2_layered():
    t0 = time.time()
    cases = []
    for d in (27.5e-6, 56e-6, 112.5e-6, 1080e-6):
        for eps_s in (0.43, 0.20):
            for interp in ("eps_new", "eps_final"):
                n_layers = 5 if d < 500e-6 else 3
                tc = time.time()
                c = layered_case(d, eps_s, interp, n_layers=n_layers)
                c["wall_s"] = time.time() - tc
                cases.append(c)
                cl = c["closure"]
                if c["n_layers"] == 0:
                    print(f"  d={d*1e6:.1f} eps_s={eps_s} {interp}: eps_final interpretation infeasible already on layer 1"); continue
                print(f"  d={d*1e6:.1f} eps_s={eps_s} {interp}: N={c['N']} steps~{len(c['t'])*5} h_end={c['h'][-1]*1e6:.0f} um infeasible_layer={c['infeasible_layer']} "
                      f"closure={None if cl is None else (cl['layer'], f'{cl['h']*1e6:.1f} um', f'{cl['t']/3600:.2f} h')} "
                      f"cum res salt {c['cumulative']['res_salt']:.1e} solid {c['cumulative']['res_solid']:.1e} ({c['wall_s']:.0f} s)")
    def thin(c, n_keep=400):
        """time series thinned to <= n_keep rows for the json (the png figures use the full series)"""
        out = {}
        for k, v in c.items():
            if k == "final_profile":
                continue
            if isinstance(v, np.ndarray) and v.ndim == 1 and v.size > n_keep:
                v = v[:: int(np.ceil(v.size / n_keep))]
            out[k] = v
        return out
    dump("s2_layered_cases", [thin(c) for c in cases])
    # table: per layer, F from (B) and from (C), eps_final, closure
    rows = []
    for c in cases:
        if c["n_layers"] == 0:
            rows.append([f"{c['d_sub']*1e6:.1f}", c["eps_s"], c["interpretation"], 1, "none", "", "", "", "", "", "eps_final = eps_s unreachable for any birth porosity"])
            continue
        for n in range(1, c["n_layers"] + 1):
            sel = c["layer_top"] == n
            if not sel.any():
                rows.append([f"{c['d_sub']*1e6:.1f}", c["eps_s"], c["interpretation"], n, f"{c['eps_new'][n]:.3f}", "not reached", "", "", "", "", ""])
                continue
            # steady part of the layer: last 20 % of its records
            idx = np.where(sel)[0]; tail = idx[int(0.8 * len(idx)):]
            # final porosity of material born in this layer = eps_final read when z_m_final lies in the layer
            zb = c["layer_bottoms"]; mat = (c["z_m_final"] >= zb[n - 1]) & (c["z_m_final"] < zb[n])
            ef = f"{np.median(c['eps_final'][mat]):.4f}" if mat.any() else "n/a"
            vC = np.mean(c["v"][tail]); vB = np.mean(c["vB"][tail]); vA = c["vA"][0]
            rows.append([f"{c['d_sub']*1e6:.1f}", c["eps_s"], c["interpretation"], n, f"{c['eps_new'][n]:.3f}",
                         f"{c['FB'][tail][0]:.3f}", f"{np.mean(c['F_prat'][tail]):.3f}", f"{vA/U_EVAP:.4f}", f"{vB/U_EVAP:.4f}", f"{vC/U_EVAP:.4f}",
                         f"{c['eps_fin_B'][tail][0]:.3f} / {ef}"])
        cl = c["closure"]
        note = "none within run" if cl is None else f"layer {cl['layer']} at h = {cl['h']*1e6:.1f} um, t = {cl['t']/3600:.2f} h"
        if c.get("infeasible_layer"):
            note += f"; eps_final = eps_s unreachable from layer {c['infeasible_layer']} on (run ends at its bottom)"
        rows.append([f"{c['d_sub']*1e6:.1f}", c["eps_s"], c["interpretation"], "closure", "", "", "", "", "", "", note])
    md("S2 Layered structures (r = 0.5, alpha = 2, c_h = 1), U_evap = 1 cm/day: per-layer quasi-steady values",
       ["d_sub [um]", "eps_s", "interp.", "layer", "eps_new", "F (B)", "F (C, Prat def.)", "v_A/U", "v_B/U", "v_C/U", "eps_final B / C"], rows)
    # figures: h(t) and v(t) for 27.5 um eps_new 0.43 and 0.20; eps_final(z_m)
    for c in cases:
        if c["n_layers"] == 0:
            continue
        tag = f"{c['d_sub']*1e6:.0f}um_eps{c['eps_s']}_{c['interpretation']}"
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        th = c["t"] / 3600
        ax[0].plot(th, c["hA"] * 1e6, "k:", label="A layer-1 Prat"); ax[0].plot(th, c["hB"] * 1e6, "b--", label="B layer-wise Prat"); ax[0].plot(th, c["h"] * 1e6, "r-", label="C numerical")
        for zb in c["layer_bottoms"][1:]:
            ax[0].axhline(zb * 1e6, color="0.8", lw=0.6)
        ax[0].set_xlabel("t [h]"); ax[0].set_ylabel("h [um]"); ax[0].legend(fontsize=7)
        ax[1].plot(th, c["vA"] / U_EVAP, "k:"); ax[1].plot(th, c["vB"] / U_EVAP, "b--"); ax[1].plot(th, c["v"] / U_EVAP, "r-")
        ax[1].set_xlabel("t [h]"); ax[1].set_ylabel("v / U_evap")
        ok = c["z_m_final"] > 0
        ax[2].plot(c["z_m_final"][ok] * 1e6, c["eps_final"][ok], "r-", label="C eps_final(z)")
        ax[2].plot(c["h"] * 1e6, c["eps_fin_B"], "b--", label="B 1-(1-eps_new)/F_n")
        ax[2].set_xlabel("z [um] (material)"); ax[2].set_ylabel("final porosity"); ax[2].legend(fontsize=7)
        if c["closure"] is not None:
            ax[2].axvline(c["closure"]["h"] * 1e6, color="r", ls=":", lw=0.8)
        fig.suptitle(f"d_sub = {c['d_sub']*1e6:.1f} um, eps_s = {c['eps_s']}, {c['interpretation']}"); fig.tight_layout()
        fig.savefig(FIG / f"s2_{tag}.png", dpi=120); plt.close(fig)
    print(f"[S2 {time.time()-t0:.0f} s]")
    return cases


# ============================================================================ S3 sensitivity / convergence
def s3_sensitivity():
    t0 = time.time()
    rows = []; out = {}
    d = 27.5e-6
    base = dict(d_sub=d, eps_s=0.43, interpretation="eps_new", n_layers=5)      # layer 5 is where (C) closes
    for c_h in (0.5, 1.0, 1.4, 2.0):
        c = layered_case(**base, c_h=c_h)
        cl = c["closure"]; out[f"c_h={c_h}"] = dict(h=c["h"][::20], t=c["t"][::20], closure=cl, layer_bottoms=c["layer_bottoms"])
        rows.append([f"c_h={c_h}", f"{c['h'][-1]*1e6:.1f}", f"{c['t'][-1]/3600:.2f}", "none" if cl is None else f"layer {cl['layer']}, h={cl['h']*1e6:.1f} um (layer bottom {c['layer_bottoms'][cl['layer']-1]*1e6:.1f}), t={cl['t']/3600:.2f} h", f"{c['cumulative']['res_salt']:.1e}"])
    for alpha in (1.5, 3.0):
        c = layered_case(**base, alpha=alpha)
        cl = c["closure"]; out[f"alpha={alpha}"] = dict(h=c["h"][::20], t=c["t"][::20], closure=cl, layer_bottoms=c["layer_bottoms"])
        rows.append([f"alpha={alpha}", f"{c['h'][-1]*1e6:.1f}", f"{c['t'][-1]/3600:.2f}", "none" if cl is None else f"layer {cl['layer']}, h={cl['h']*1e6:.1f} um (layer bottom {c['layer_bottoms'][cl['layer']-1]*1e6:.1f}), t={cl['t']/3600:.2f} h", f"{c['cumulative']['res_salt']:.1e}"])
    md("S3a Layer-thickness sensitivity (27.5 um, eps_new = 0.43, up to 5 layers)", ["setting", "h_end [um]", "t_end [h]", "closure", "cum. salt residual"], rows)
    rows = []
    ref = None
    for nb, Lo, st in ((20, 60.0, 1.03), (40, 60.0, 1.03), (80, 60.0, 1.03), (40, 120.0, 1.03), (40, 60.0, 1.015)):
        c = layered_case(**base, n_band=nb, L_over=Lo, stretch=st)
        # F of layer 3 (quasi-steady tail) and closure
        sel = np.where(c["layer_top"] == 3)[0]; tail = sel[int(0.8 * len(sel)):] if len(sel) else sel
        F3 = np.mean(c["F_prat"][tail]) if len(tail) else np.nan
        zb = c["layer_bottoms"]; mat = (c["z_m_final"] >= zb[1]) & (c["z_m_final"] < zb[2])
        ef2 = np.median(c["eps_final"][mat]) if mat.any() else np.nan
        cl = c["closure"]
        rows.append([f"n_band={nb}, L={Lo:.0f} dC, stretch={st}", c["N"], f"{F3:.5f}", f"{ef2:.5f}", f"{c['t'][-1]/3600:.3f}",
                     "none" if cl is None else f"layer {cl['layer']}, h={cl['h']*1e6:.2f} um, t={cl['t']/3600:.3f} h"])
        out[f"grid_nb{nb}_L{Lo:.0f}_s{st}"] = dict(N=c["N"], F3=F3, eps_final_layer2=ef2, t_end=c["t"][-1], closure=cl)
    md("S3b Grid / window convergence of the layered run (27.5 um, eps_new = 0.43)", ["grid", "N", "F layer 3", "eps_final layer 2", "t_end [h]", "closure"], rows)
    dump("s3_sensitivity", out)
    print(f"[S3 {time.time()-t0:.0f} s]")


if __name__ == "__main__":
    T0 = time.time()
    s1_homogeneous()
    s2_layered()
    s3_sensitivity()
    (RES / "summary.md").write_text("# Generated tables (run_validation.py)\n" + "\n".join(SUMMARY) + "\n")
    print(f"total {time.time()-T0:.0f} s; tables in results/summary.md")
