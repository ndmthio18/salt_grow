"""Step 1B: the 20 um -> 10 um refinement check of A1 and A2 at h = 0.3 mm (gas only).

    OMP_NUM_THREADS=1 $PY runs/S_surface/step1_grid_refinement/run.py        # ~20 s, ~0.6 GB

The two 10 um calculations use saltv2.cases.step1_frozen_salt.run_case unchanged (CG, rtol 1e-12).
The 20 um references are READ from runs/S_surface/step1_frozen_salt/results/step1_results.json
(the Step 1 record, commit 57903db); nothing at 20 um is recomputed.  Writes
results/step1_refinement.json and results/step1_refinement.csv.
"""
from __future__ import annotations

import csv
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from saltv2.cases import step1_frozen_salt as s1                               # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
REF_FILE = HERE.parent / "step1_frozen_salt" / "results" / "step1_results.json"
A_FINE = 10e-6
H = 0.3e-3
KEYS = ("E_total", "E_top", "E_side", "E_farfield", "A_top_wet", "A_side_wet", "A_total_wet",
        "j_top_mean", "j_side_mean", "j_wet_mean", "imbalance", "residual", "n_unknown")


def geometry_asserts(a):
    """The 10 um grid must represent the same physical opening, height and faces as 20 um."""
    grid, foot, L_ext = s1.build(s1.H_BOX, a)
    n_side = int(round(s1.PORE_SIDE / a)); n_h = int(round(H / a))
    assert grid.shape == (int(round(s1.L_BOX / a)),) * 2 + (int(round(s1.H_BOX / a)),), grid.shape
    assert int(foot.sum()) == n_side ** 2, foot.sum()
    rows = np.where(foot.any(axis=1))[0]; cols = np.where(foot.any(axis=0))[0]
    assert rows.size == n_side and cols.size == n_side and foot[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1].all()
    pillar, faces = s1.wet_faces(grid, foot, H, sides=True)
    assert int(pillar.sum()) == n_side ** 2 * n_h
    assert int(faces[(2, -1)].sum()) == n_side ** 2
    n_side_faces = sum(int(m.sum()) for k, m in faces.items() if k[0] != 2)
    assert n_side_faces == 4 * n_side * n_h
    A_top = n_side ** 2 * a * a; A_side = n_side_faces * a * a; V = n_side ** 2 * n_h * a ** 3
    assert abs(A_top - s1.PORE_SIDE ** 2) <= 1e-12 * s1.PORE_SIDE ** 2
    assert abs(A_side - 4 * s1.PORE_SIDE * H) <= 1e-12 * 4 * s1.PORE_SIDE * H
    assert abs(V - s1.PORE_SIDE ** 2 * H) <= 1e-12 * s1.PORE_SIDE ** 2 * H
    assert abs(L_ext - 18.6e-3) < 1e-15
    return dict(shape=list(grid.shape), n_cells=int(np.prod(grid.shape)), n_side=n_side, n_h=n_h,
                n_pillar_voxels=int(pillar.sum()), n_top_faces=n_side ** 2, n_side_faces=n_side_faces,
                A_top=A_top, A_side=A_side, V_pillar=V, L_ext=L_ext)


def main():
    ref_all = json.loads(REF_FILE.read_text())
    ref = {c: ref_all["primary"][f"{c}_h0.3"] for c in ("A1", "A2")}
    out = dict(reference=dict(file=str(REF_FILE.relative_to(HERE.parents[2])), sha256=hashlib.sha256(REF_FILE.read_bytes()).hexdigest(),
                              commit="57903db886478dd5de4085acdf500a9ae820cc8e", solver=ref_all["solver"], a=ref_all["settings"]["a"]),
               settings=dict(T=295.15, RH=0.07, C=0.264, pore_side=s1.PORE_SIDE, L_box=s1.L_BOX, delta_eff=s1.DELTA_EFF,
                             h_box=s1.H_BOX, h=H, a_fine=A_FINE, cg_rtol=s1.CG_RTOL),
               geometry={}, fine={}, coarse={k: {q: ref[k][q] for q in KEYS + ("top_max_over_mean", "side_max_over_mean")} for k in ref})
    for a in (ref_all["settings"]["a"], A_FINE):
        out["geometry"][f"{a*1e6:.0f}um"] = geometry_asserts(a)
    print("geometry 10 um:", out["geometry"]["10um"])
    for cond in ("A1", "A2"):
        t0 = time.time(); r0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        o, X, j, m = s1.run_case(cond, H, a=A_FINE, h_box=s1.H_BOX, solver="cg")
        o["wall_s"] = time.time() - t0; o["peak_rss_MB"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e3
        assert o["E_total"] == o["E_top"] + o["E_side"]
        assert np.all(np.isnan(X[m["foot"], 0]))                           # the covered opening is not a gas cell
        if cond == "A1":
            assert o["E_side"] == 0.0 and o["n_side_faces"] == 0
        out["fine"][cond] = {k: v for k, v in o.items() if k not in ("shape",)}
        print(f"{cond} 10 um: E_total {o['E_total']:.12e} top {o['E_top']:.12e} side {o['E_side']:.12e} far {o['E_farfield']:.12e} "
              f"imb {o['imbalance']:+.2e} res {o['residual']:.2e} CG it {o['cg_iterations']} n {o['n_unknown']} {o['wall_s']:.1f} s {o['peak_rss_MB']:.0f} MB")
    # paired metrics at both resolutions and the change 20 -> 10
    def paired(d):
        a1, a2 = d["A1"], d["A2"]
        return dict(dE=a2["E_total"] - a1["E_total"], R_sides=a2["E_total"] / a1["E_total"],
                    f_side=a2["E_side"] / a2["E_total"], R_sides_minus_1_pct=100 * (a2["E_total"] / a1["E_total"] - 1))
    out["paired"] = dict(coarse=paired(out["coarse"]), fine=paired(out["fine"]))
    comp = {}
    for cond in ("A1", "A2"):
        for k in ("E_total", "E_top", "E_side", "j_top_mean", "j_side_mean"):
            c, f = out["coarse"][cond][k], out["fine"][cond][k]
            if c is None or f is None or c == 0.0:
                continue
            comp[f"{cond}.{k}"] = dict(coarse=c, fine=f, rel_change=(f - c) / c)
    for k in ("dE", "R_sides", "f_side"):
        c, f = out["paired"]["coarse"][k], out["paired"]["fine"][k]
        comp[f"paired.{k}"] = dict(coarse=c, fine=f, abs_change=f - c, rel_change=(f - c) / c)
    comp["paired.R_sides_minus_1"] = dict(coarse=out["paired"]["coarse"]["R_sides"] - 1, fine=out["paired"]["fine"]["R_sides"] - 1,
                                          rel_change=(out["paired"]["fine"]["R_sides"] - out["paired"]["coarse"]["R_sides"]) / (out["paired"]["coarse"]["R_sides"] - 1))
    out["comparison"] = comp
    print("\n20 um -> 10 um")
    for k, v in comp.items():
        print(f"  {k:22s} 20um {v['coarse']:.12e}  10um {v['fine']:.12e}  rel change {v['rel_change']:+.3e}"
              + (f"  abs change {v['abs_change']:+.3e}" if "abs_change" in v else ""))
    RES.mkdir(exist_ok=True)
    (RES / "step1_refinement.json").write_text(json.dumps(out, indent=1, default=float))
    with (RES / "step1_refinement.csv").open("w", newline="") as fcsv:
        w = csv.writer(fcsv); w.writerow(["case", "a_um"] + list(KEYS) + ["cg_iterations", "wall_s", "peak_rss_MB"])
        for cond in ("A1", "A2"):
            w.writerow([cond, 20] + [out["coarse"][cond][k] for k in KEYS] + ["", "", ""])
            w.writerow([cond, 10] + [out["fine"][cond][k] for k in KEYS] + [out["fine"][cond]["cg_iterations"], out["fine"][cond]["wall_s"], out["fine"][cond]["peak_rss_MB"]])
        w.writerow([]); w.writerow(["paired", "a_um", "dE", "R_sides", "f_side"])
        for lab, a in (("coarse", 20), ("fine", 10)):
            w.writerow(["A2-A1", a, out["paired"][lab]["dE"], out["paired"][lab]["R_sides"], out["paired"][lab]["f_side"]])


if __name__ == "__main__":
    main()
