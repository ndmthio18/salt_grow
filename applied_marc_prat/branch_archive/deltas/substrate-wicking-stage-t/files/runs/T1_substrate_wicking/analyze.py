"""Post-processing of Stage-T batches: paired B2 / C differences against the time-step error.

    python runs/T1_substrate_wicking/analyze.py OUT_DIR BATCH_DIR [BATCH_DIR ...]

Reads the records of every BATCH_DIR (a case name must occur in one batch only) and writes
OUT_DIR/pairs.md and OUT_DIR/fig_*.png.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

C_SAT = 0.264
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def load(dirs):
    S, A = {}, {}
    for d in map(Path, dirs):
        for p in sorted(d.glob("*.json")):
            if p.stem in ("manifest", "PROVENANCE"):
                continue
            if p.stem in S:
                raise SystemExit(f"{p.stem} occurs in two batches")
            S[p.stem] = json.loads(p.read_text()); S[p.stem]["batch"] = d.name
            if (d / f"{p.stem}.npz").exists():
                A[p.stem] = np.load(d / f"{p.stem}.npz")
    return S, A


def tsat(s):
    return np.nan if s is None or s["censored"] else s["t_sat_first_h"]


def tcross_stats(s):
    t = np.array([np.nan if v is None else v for v in s["t_cross_h"]], float)
    return np.nanmin(t), np.nanmedian(t), np.nanmax(t), int(np.sum(np.isfinite(t))), t.size


def main(out, *dirs):
    d = Path(out); d.mkdir(parents=True, exist_ok=True)
    S, A = load(dirs)
    L = []
    P = L.append
    P("# Stage-T analysis of batches " + ", ".join(Path(x).name for x in dirs) + "\n")
    # ---- time-step convergence of the key observables (hotspot, both variants)
    P("## Time-step refinement (hotspot, 5x5x50, conservative)\n")
    P("| variant | dt [s] | t_sat first [h] | var w_surf @1h | max w_surf @1h | hot salt lateral @1h [kg/s] |")
    P("|---|---|---|---|---|---|")
    conv = {}
    for v in ("B2", "C"):
        for dt in (60, 30, 15, 7.5):
            name = f"{v}_hotspot" if dt == 30 else f"{v}_hotspot_dt{dt:g}"
            s = S.get(name)
            if s is None:
                continue
            s1 = s["snapshots"]["t1h"]
            conv.setdefault(v, []).append((dt, tsat(s), s1["w_surf_var"], s1["w_surf_max"]))
            P(f"| {v} | {dt:g} | {tsat(s):.5f} | {s1['w_surf_var']:.6e} | {s1['w_surf_max']:.7f} | "
              f"{s1['hot_salt_lateral']:+.4e} |")
    P("")
    err = {}
    for v, rows in conv.items():
        rows = sorted(rows)                       # dt ascending
        if len(rows) >= 3:
            (d1, t1, v1, m1), (d2, t2, v2, m2), (d3, t3, v3, m3) = rows[:3]
            p_ord = np.log(abs((t3 - t2) / (t2 - t1))) / np.log(2.0) if t2 != t1 else np.nan
            rich = t1 + (t1 - t2) / (2 ** p_ord - 1) if np.isfinite(p_ord) and p_ord > 0 else np.nan
            e30 = abs(dict((r[0], r[1]) for r in rows)[30] - rich)
            err[v] = dict(order=p_ord, t_rich=rich, err_dt30=e30,
                          var_err30=abs(dict((r[0], r[2]) for r in rows)[30] - v1))
            P(f"- {v}: observed order of t_sat in dt = {p_ord:.2f}; Richardson t_sat = {rich:.5f} h; "
              f"|t_sat(dt=30) - Richardson| = {e30*3600:.2f} s; |var(dt30) - var(dt7.5)| @1h = "
              f"{err[v]['var_err30']:.3e}")
    P("")
    # ---- paired differences
    P("## B2 vs C pairs (dt = 30 s, conservative)\n")
    P("| pair | t_sat B2 [h] | t_sat C [h] | C - B2 [s] | (C-B2)/B2 | dt-error bound [s] | "
      "var@1h B2 | var@1h C | var ratio C/B2 | sub-lat liquid/E C @1h | crossings B2 min/med/max [h] | C min/med/max [h] |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|")
    ebound = max((e["err_dt30"] for e in err.values()), default=np.nan)
    pairs = [("uniform", "B2_uniform", "C_uniform"), ("stripe", "B2_stripe", "C_stripe"),
             ("hotspot", "B2_hotspot", "C_hotspot"), ("cv25", "B2_cv25", "C_cv25"),
             ("hotspot L10 (7 mm)", "B2_hotspot_L10", "C_hotspot_L10"),
             ("hotspot L25 (17.5 mm)", "B2_hotspot_L25", "C_hotspot_L25"),
             ("hotspot n10", "B2_hotspot_n10", "C_hotspot_n10"),
             ("hotspot n20", "B2_hotspot_n20", "C_hotspot_n20"),
             ("hotspot production", "B2_hotspot_production", "C_hotspot_production")]
    for lab, b, c in pairs:
        if b not in S or c not in S:
            P(f"| {lab} | missing | | | | | | | | | | |")
            continue
        sb, sc = S[b], S[c]
        tb, tc = tsat(sb), tsat(sc)
        vb, vc = sb["snapshots"]["t1h"]["w_surf_var"], sc["snapshots"]["t1h"]["w_surf_var"]
        cb, cc = tcross_stats(sb), tcross_stats(sc)
        P(f"| {lab} | {tb:.5f} | {tc:.5f} | {(tc-tb)*3600:+.1f} | {(tc-tb)/tb:+.3%} | {ebound*3600:.2f} | "
          f"{vb:.3e} | {vc:.3e} | {vc/vb if vb > 0 else float('nan'):.4f} | "
          f"{sc['snapshots']['t1h']['subsurface_lateral_liquid_over_E']:.3e} | "
          f"{cb[0]:.3f}/{cb[1]:.3f}/{cb[2]:.3f} ({cb[3]}/{cb[4]}) | {cc[0]:.3f}/{cc[1]:.3f}/{cc[2]:.3f} ({cc[3]}/{cc[4]}) |")
    P("")
    # ---- controls and conservation
    P("## Controls, conservation, supply (dt = 30 s)\n")
    P("| case | mode | t_sat [h] | E_cum to end [kg] | liquid in [kg] | salt in [kg] | res salt [kg] | "
      "rel res salt (of salt gained) | rel res water (of evaporated) | max step res salt [kg] | "
      "p_visc max @t_sat [Pa] | p_req G1 / G2 max [Pa] | P_open / P_throat,min [Pa] | margin G1 / G2 |")
    P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name in sorted(S):
        s = S[name]
        sf = s["snapshots"].get("first_saturation", s["snapshots"]["end"])
        cu = s["cumulative"]
        P(f"| {name} | {s['mode']} | {tsat(s):.4f} | {cu['evaporated']:.4e} | {cu['liquid_in']:.4e} | "
          f"{cu['salt_in']:.4e} | {cu['res_salt']:+.2e} | {s['rel_res_salt']:+.2e} | {s['rel_res_water']:+.2e} | "
          f"{s['max_step_residual']['salt']:.1e} | {sf['p_visc_max']:.2e} | {sf['p_req_G1_max']:.1f} / "
          f"{sf['p_req_G2_max']:.1f} | {s['capillary']['P_open_min']:.0f} / {s['capillary']['P_throat_min']:.0f} | "
          f"{sf['margin_open_G1']:.3f} / {sf['margin_open_G2']:.3f} |")
    P("")
    # ---- supply distribution at 1 h
    P("## Per-pore supply at t = 1 h (hotspot pore and spread over the surface pores)\n")
    P("| case | hot demand [kg/s] | hot liquid from below | hot liquid lateral | hot salt from below | "
      "hot salt lateral | CV liquid from below | CV salt from below | CV reservoir liquid | CV reservoir salt |")
    P("|---|---|---|---|---|---|---|---|---|---|")
    for name in ("A_hotspot", "B1_hotspot", "B2_hotspot", "C_hotspot", "B2_stripe", "C_stripe",
                 "B2_cv25", "C_cv25", "A_uniform", "B1_uniform", "B2_uniform", "C_uniform"):
        if name not in S:
            continue
        snaps = S[name]["snapshots"]
        lab = "t1h" if "t1h" in snaps else "t0.25h"     # A saturates everywhere before 1 h
        s1 = snaps[lab]
        P(f"| {name}{'' if lab == 't1h' else ' (@0.25 h)'} | {s1['hot_demand']:.3e} | {s1['hot_liquid_from_below']:.3e} | "
          f"{s1['hot_liquid_lateral']:+.3e} | {s1['hot_salt_from_below']:.3e} | {s1['hot_salt_lateral']:+.3e} | "
          f"{s1['liquid_from_below_cv']:.4f} | {s1['salt_from_below_cv']:.4f} | "
          f"{s1['reservoir_liquid_cv']:.3e} | {s1['reservoir_salt_cv']:.3e} |")
    # ---- physical depth (E) and lateral domain (D)
    P("\n## E: physical depth at fixed pore scale and lateral area (hotspot, 5x5)\n")
    P("| layers | depth [mm] | t_sat B2 [h] | t_sat C [h] | C - B2 [s] | var@1h B2 | var@1h C | p_req G1 max [Pa] |")
    P("|---|---|---|---|---|---|---|---|")
    for nl in (5, 10, 25, 50, 100):
        b, c = (f"{v}_hotspot" + ("" if nl == 50 else f"_L{nl}") for v in ("B2", "C"))
        if b in S and c in S:
            sb, sc = S[b], S[c]
            P(f"| {nl} | {nl*0.7:.1f} | {tsat(sb):.5f} | {tsat(sc):.5f} | {(tsat(sc)-tsat(sb))*3600:+.1f} | "
              f"{sb['snapshots']['t1h']['w_surf_var']:.3e} | {sc['snapshots']['t1h']['w_surf_var']:.3e} | "
              f"{sc['snapshots']['first_saturation']['p_req_G1_max']:.1f} |")
    P("\n## D: lateral domain at fixed pore scale (50 layers)\n")
    P("| forcing | n | t_sat B2 [h] | t_sat C [h] | C - B2 [s] | (C-B2)/B2 | var@1h B2 | var@1h C | wall C [s] |")
    P("|---|---|---|---|---|---|---|---|---|")
    for pat in ("hotspot", "hotspotbg"):
        for n in (5, 10, 20):
            if n == 5:
                b, c = "B2_hotspot", "C_hotspot"
            else:
                b, c = f"B2_{pat}_n{n}", f"C_{pat}_n{n}"
            if b in S and c in S:
                sb, sc = S[b], S[c]
                P(f"| {pat} | {n} | {tsat(sb):.5f} | {tsat(sc):.5f} | {(tsat(sc)-tsat(sb))*3600:+.1f} | "
                  f"{(tsat(sc)-tsat(sb))/tsat(sb):+.3%} | {sb['snapshots']['t1h']['w_surf_var']:.3e} | "
                  f"{sc['snapshots']['t1h']['w_surf_var']:.3e} | {sc['wall_s']:.0f} |")
    P("\n## Pore-size heterogeneity cv = 0.25, uniform demand: realisations and dt check\n")
    P("| case | t_sat B2 [h] | t_sat C [h] | C - B2 [s] | var@1h B2 | var@1h C | C/B2 | first pore B2 / C | r_first B2 / C [um] |")
    P("|---|---|---|---|---|---|---|---|---|")
    signs = []
    for suf in ("", "_dt15", "_s2", "_s3", "_s4", "_s5"):
        b, c = f"B2_cv25{suf}", f"C_cv25{suf}"
        if b in S and c in S:
            sb, sc = S[b], S[c]
            vb, vc = sb["snapshots"]["t1h"]["w_surf_var"], sc["snapshots"]["t1h"]["w_surf_var"]
            if suf != "_dt15":
                signs.append((tsat(sc) - tsat(sb), vc / vb))
            P(f"| cv25{suf or ' (seed 1)'} | {tsat(sb):.5f} | {tsat(sc):.5f} | {(tsat(sc)-tsat(sb))*3600:+.1f} | "
              f"{vb:.3e} | {vc:.3e} | {vc/vb:.3f} | {sb['first_pore']} / {sc['first_pore']} | "
              f"{_r_first(A, b, sb):.1f} / {_r_first(A, c, sc):.1f} |")
    if signs:
        dts = np.array([x[0] for x in signs]) * 3600; vr = np.array([x[1] for x in signs])
        P(f"\n- {len(signs)} realisations: C - B2 in t_sat = {dts.mean():+.1f} s mean, range "
          f"[{dts.min():+.1f}, {dts.max():+.1f}] s, C earlier in {int(np.sum(dts < 0))}/{len(dts)}; "
          f"var ratio C/B2 at 1 h mean {vr.mean():.3f}, range [{vr.min():.3f}, {vr.max():.3f}]")
    P("\n## Extra dt check (15 s vs 30 s) of the stripe and cv25 pairs\n")
    for pat in ("stripe", "cv25"):
        r = []
        for v in ("B2", "C"):
            if f"{v}_{pat}" in S and f"{v}_{pat}_dt15" in S:
                r.append((v, tsat(S[f"{v}_{pat}"]), tsat(S[f"{v}_{pat}_dt15"])))
        if len(r) == 2:
            d30 = (r[1][1] - r[0][1]) * 3600; d15 = (r[1][2] - r[0][2]) * 3600
            P(f"- {pat}: C - B2 = {d30:+.1f} s at dt 30, {d15:+.1f} s at dt 15; per-variant shift "
              f"{(r[0][2]-r[0][1])*3600:+.2f} / {(r[1][2]-r[1][1])*3600:+.2f} s")
    P("\n## Subsurface lateral liquid flux per subsurface layer at 1 h (sum |flux| / (layers - 1) / E)\n")
    for name in ("C_stripe", "C_hotspot", "C_cv25", "C_hotspot_L10", "C_hotspot_L100"):
        if name in S:
            s1 = S[name]["snapshots"]["t1h"]; nl = S[name]["n_layers"]
            P(f"- {name}: total {s1['subsurface_lateral_liquid_over_E']:.3e}, per layer "
              f"{s1['subsurface_lateral_liquid_over_E'] / (nl - 1):.3e}")
    (d / "pairs.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))

    # ---- figures
    want = [("B2_stripe", "C_stripe"), ("B2_hotspot", "C_hotspot"), ("B2_cv25", "C_cv25")]
    want = [w for w in want if all(k in A for k in w)]
    if want:
        fig, ax = plt.subplots(len(want), 3, figsize=(10, 3.1 * len(want)), squeeze=False)
        for i, (b, c) in enumerate(want):
            n = S[b]["n_side"]
            wb = A[b]["snap_t1h_w"][:n * n].reshape(n, n); wc = A[c]["snap_t1h_w"][:n * n].reshape(n, n)
            lo, hi = min(wb.min(), wc.min()), max(wb.max(), wc.max())
            for j, (fld, ttl) in enumerate(((wb, b), (wc, c))):
                im = ax[i, j].imshow(fld.T, origin="lower", vmin=lo, vmax=hi, cmap="viridis")
                ax[i, j].set_title(f"{ttl}: w_surf @1 h"); fig.colorbar(im, ax=ax[i, j], fraction=0.046)
            im = ax[i, 2].imshow((wc - wb).T, origin="lower", cmap="RdBu_r",
                                 vmin=-abs(wc - wb).max(), vmax=abs(wc - wb).max())
            ax[i, 2].set_title("C - B2"); fig.colorbar(im, ax=ax[i, 2], fraction=0.046)
        fig.tight_layout(); fig.savefig(d / "fig_surface_maps_1h.png", dpi=110); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for name, ls in (("B2_hotspot", "--"), ("C_hotspot", "-"), ("B2_stripe", "--"), ("C_stripe", "-"),
                     ("A_hotspot", ":"), ("B1_hotspot", "-.")):
        if name not in A:
            continue
        t = A[name]["rec_t"] / 3600.0
        ax[0].plot(t, A[name]["rec_w_surf"].max(axis=1), ls, label=name)
        ax[1].plot(t, A[name]["rec_var_surf"], ls, label=name)
        E = np.gradient(A[name]["rec_evaporated"], A[name]["rec_t"])
        ax[2].plot(t, A[name]["rec_liq_sub_abs"] / E, ls, label=name)
    ax[0].axhline(C_SAT, color="k", lw=0.8); ax[0].set_ylabel("max surface w")
    ax[1].set_ylabel("var surface w"); ax[1].set_yscale("log")
    ax[2].set_ylabel("sum |subsurface lateral liquid flux| / E")
    for a in ax:
        a.set_xlabel("t [h]"); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(d / "fig_timeseries.png", dpi=110); plt.close(fig)


def _r_first(A, name, s):
    """Inscribed radius [um] of the first surface pore to saturate.  The record holds no radii;
    the network is deterministic in (variant, n_side, n_layers, cv, seed), so it is rebuilt."""
    if s["first_pore"] is None:
        return float("nan")
    from saltv2.cases import stage_t_wicking as st
    layer = st.build_network(s["variant"], s["n_side"], s["n_layers"], s["cv"], s["seed"])
    r = layer.r[layer.surface]
    return float(r[s["first_pore"]] * 1e6)


if __name__ == "__main__":
    main(sys.argv[1], *sys.argv[2:])
