"""Step 1 launcher: evaporation from a PRESCRIBED wet salt pillar (saltv2.cases.step1_frozen_salt).

    OMP_NUM_THREADS=1 $PY runs/S_surface/step1_frozen_salt/run.py [--solver direct|cg]

Primary matrix (20 um voxels, 2 mm box): A0; A1 and A2 at h = 0.1, 0.3, 0.5 mm (seven runs), plus
A1/A2 at h = 0 (continuity, same faces as A0 by construction).  Checks: the face-selective assembly
against solve_vapor (50 um and 20 um), the 50 um grid at A0 / A1 / A2 h = 0.3 mm, and A2 h = 0.5 mm
in a 4 mm box.  Writes results/step1_results.json, step1_table.csv, step1_fields.npz, step1_summary.md.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from saltv2.cases import step0_gas as s0                                       # noqa: E402
from saltv2.cases import step1_frozen_salt as s1                               # noqa: E402

RES = Path(__file__).resolve().parent / "results"
LINES = ["# run.py output (Step 1, prescribed wet salt pillar, gas only)", ""]
COLS = ("run", "cond", "h_m", "a_m", "h_box_m", "E_total", "E_top", "E_side", "E_farfield", "A_projected",
        "A_open_phys", "A_open_raster", "A_top_wet", "A_side_wet_phys", "A_side_wet", "A_total_wet",
        "j_projected", "j_top_mean", "j_side_mean", "j_wet_mean", "imbalance", "residual", "n_unknown",
        "top_max_over_mean", "side_max_over_mean", "wall_s")


def say(s=""):
    print(s, flush=True); LINES.append(s)


def one(name, cond, h, a, h_box, solver, fields=None, keep=()):
    t0 = time.time()
    out, X, j_face, masks = s1.run_case(cond, h, a=a, h_box=h_box, solver=solver)
    out["wall_s"] = time.time() - t0; out["run"] = name
    if fields is not None and name in keep:
        jc = X.shape[1] // 2
        fields[f"X_{name}"] = X[:, jc, :]                                  # section through the pore centre
        fields[f"pillar_{name}"] = masks["pillar"][:, jc, :]
        fields[f"jtop_{name}"] = j_face[(2, -1)][..., int(round(h / a))]     # lateral map of the top-face flux
        for key, m in masks["faces"].items():
            if key[0] != 2:
                fields[f"jside_{name}_{key[0]}_{key[1]:+d}"] = np.where(m, j_face[key], np.nan)[:, jc, :]
        fields[f"wet_{name}"] = np.stack([masks["faces"].get(k, np.zeros(X.shape, bool))[:, jc, :]
                                          for k in ((2, -1), (0, +1), (0, -1))])
    say(f"  {name:14s} {cond}  h {h*1e3:4.2f} mm  a {a*1e6:2.0f} um  box {h_box*1e3:3.1f} mm  "
        f"E_tot {out['E_total']:.10e}  top {out['E_top']:.6e}  side {out['E_side']:.6e}  far {out['E_farfield']:.10e}  "
        f"imb {out['imbalance']:+.1e}  res {out['residual']:.1e}  n {out['n_unknown']:6d}  {out['wall_s']:6.1f} s")
    return out


def main(solver):
    t_all = time.time()
    results, fields, rows = dict(solver=solver, settings=dict(T=295.15, RH=0.07, C=0.264, pore_side=s1.PORE_SIDE,
                                 L_box=s1.L_BOX, delta_eff=s1.DELTA_EFF, a=s1.A_VOXEL, a_coarse=s1.A_VOXEL_COARSE,
                                 h_box=s1.H_BOX, heights=list(s1.HEIGHTS))), {}, []
    keep = ("A1_h0.3", "A2_h0.3")
    say("## 0. face-selective assembly against the production kernel (A2 = SALT_WET pillar; h = 0 = bottom opening)")
    for a, sol in ((s1.A_VOXEL_COARSE, "direct"), (s1.A_VOXEL, solver)):
        t0 = time.time(); cc = s1.cross_check(h=0.3e-3, a=a, solver=sol)
        results[f"cross_check_{int(a*1e6)}um"] = cc
        say(f"  a = {a*1e6:2.0f} um ({sol}): A2 E rel {cc['A2_E_rel']:+.1e}, X max rel {cc['A2_X_rel_max']:.1e}, "
            f"E_top rel {cc['A2_E_top_rel']:+.1e}, kernel base flux {cc['A2_E_base_kernel']:.1e}; "
            f"h=0 E rel {cc['A0_E_rel']:+.1e}, X max rel {cc['A0_X_rel_max']:.1e}  [{time.time()-t0:.1f} s]")
    say("\n## 1. primary matrix: 20 um voxels, 2 mm box (far field 20.6 mm)")
    prim = {}
    prim["A0"] = one("A0", "A0", 0.0, s1.A_VOXEL, s1.H_BOX, solver)
    for cond in ("A1", "A2"):
        for h in s1.HEIGHTS:
            prim[f"{cond}_h{h*1e3:.1f}"] = one(f"{cond}_h{h*1e3:.1f}", cond, h, s1.A_VOXEL, s1.H_BOX, solver, fields, keep)
    results["primary"] = prim
    E0 = prim["A0"]["E_total"]
    say("\n  ratios (E0 = A0 bare opening)")
    say("  h [mm]   R_height = E_A1/E0   R_sides = E_A2/E_A1   R_total = E_A2/E0   1-D control d/(d-h)   A2 side share   j_top_mean A1 / A2   j_side_mean A2")
    ratios = {}
    for h in s1.HEIGHTS:
        k = f"h{h*1e3:.1f}"; a1, a2 = prim[f"A1_{k}"], prim[f"A2_{k}"]
        r = dict(h=h, R_height=a1["E_total"] / E0, R_sides=a2["E_total"] / a1["E_total"], R_total=a2["E_total"] / E0,
                 control_1d=a1["ratio_1d_control"], side_share=a2["E_side"] / a2["E_total"],
                 A_total_wet_over_A0=a2["A_total_wet"] / a2["A_open_raster"])
        ratios[k] = r
        say(f"  {h*1e3:4.2f}     {r['R_height']:.6f}             {r['R_sides']:.6f}             {r['R_total']:.6f}         "
            f"{r['control_1d']:.6f}            {r['side_share']:.4f}        {a1['j_top_mean']:.4e} / {a2['j_top_mean']:.4e}   "
            f"{a2['j_side_mean'] if a2['j_side_mean'] else 0:.4e}")
    results["ratios"] = ratios
    say("\n## 2. grid comparison: 50 um (the Step 0 grid, 10 x 10 cells = 500 um) at A0, A1 and A2 h = 0.3 mm")
    coarse = {}
    for name, cond, h in (("A0", "A0", 0.0), ("A1_h0.3", "A1", 0.3e-3), ("A2_h0.3", "A2", 0.3e-3)):
        coarse[name] = one(name + "_50um", cond, h, s1.A_VOXEL_COARSE, s1.H_BOX, "direct")
        f, c = prim[name], coarse[name]
        say(f"      vs 20 um: E_total {(c['E_total'] - f['E_total']) / f['E_total']:+.3e}, j_top_mean {(c['j_top_mean'] - f['j_top_mean']) / f['j_top_mean']:+.3e}"
            + (f", j_side_mean {(c['j_side_mean'] - f['j_side_mean']) / f['j_side_mean']:+.3e}" if c["j_side_mean"] else "")
            + f"; A_open_raster {c['A_open_raster']:.3e} vs {f['A_open_raster']:.3e} (phys {f['A_open_phys']:.3e})")
    results["coarse"] = coarse
    E_step0 = s0.solve_pore(s1.H_BOX)[0]["E"]                  # the Step 0 calculation itself (50 um, direct)
    results["step0_A0_50um_reference"] = E_step0
    say(f"  Step 0 regression: A0 50 um E = {coarse['A0']['E_total']:.10e} vs Step 0 {E_step0:.10e} "
        f"({(coarse['A0']['E_total'] - E_step0) / E_step0:+.1e})")
    say("\n## 3. domain truncation: A2 h = 0.5 mm in a 4 mm box (L_ext 16.6 mm) against the 2 mm box")
    big = one("A2_h0.5_box4", "A2", 0.5e-3, s1.A_VOXEL, s1.H_BOX_CHECK, solver, fields, ("A2_h0.5_box4",))
    ref = prim["A2_h0.5"]
    results["box_check"] = dict(box4=big, rel=dict((k, (big[k] - ref[k]) / ref[k]) for k in ("E_total", "E_top", "E_side", "j_top_mean", "j_side_mean")))
    say("  (4 mm - 2 mm) / 2 mm: " + ", ".join(f"{k} {v:+.2e}" for k, v in results["box_check"]["rel"].items()))
    # local flux distributions: the top-face map and the side faces of the two boxes
    o5, X5, j5, m5 = s1.run_case("A2", 0.5e-3, a=s1.A_VOXEL, h_box=s1.H_BOX, solver=solver)
    o5b, X5b, j5b, m5b = s1.run_case("A2", 0.5e-3, a=s1.A_VOXEL, h_box=s1.H_BOX_CHECK, solver=solver)
    k5 = int(round(0.5e-3 / s1.A_VOXEL))
    d_top = np.abs(j5[(2, -1)][..., k5] - j5b[(2, -1)][..., k5]).max() / j5[(2, -1)][..., k5].max()
    d_side = max(np.abs(j5[k][m5["faces"][k]] - j5b[k][m5b["faces"][k]]).max() / j5[k][m5["faces"][k]].max()
                 for k in m5["faces"] if k[0] != 2)
    results["box_check"]["local_flux_max_rel_diff"] = dict(top=float(d_top), side=float(d_side))
    say(f"  local flux, max |4 mm - 2 mm| / max: top face {d_top:.2e}, side faces {d_side:.2e}")
    for d in (prim, coarse, {"A2_h0.5_box4": big}):
        for o in d.values():
            rows.append({c: o.get(c) for c in COLS})
    say(f"\n[{time.time() - t_all:.0f} s total]")
    RES.mkdir(exist_ok=True)
    (RES / "step1_results.json").write_text(json.dumps(results, indent=1, default=float))
    np.savez_compressed(RES / "step1_fields.npz", **fields)
    with (RES / "step1_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader(); w.writerows(rows)
    (RES / "step1_summary.md").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main("cg" if "--solver" in sys.argv and sys.argv[sys.argv.index("--solver") + 1] == "cg" else "direct")
