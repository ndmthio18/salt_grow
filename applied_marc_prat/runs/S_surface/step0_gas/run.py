"""Step 0 launcher: the four gas-only calculations of saltv2.cases.step0_gas plus the
supplementary rasterisation diagnostic.

    OMP_NUM_THREADS=1 $PY runs/S_surface/step0_gas/run.py        # seconds

Writes results/step0_results.json (every quantitative output), results/step0_table.csv (one row per
calculation), results/step0_fields.npz (the vapour fields for plot.py) and results/step0_summary.md.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from saltv2.cases import step0_gas as g                                        # noqa: E402

RES = Path(__file__).resolve().parent / "results"
LINES = ["# run.py output (Step 0, gas-only single-pore evaporation)", ""]
COLUMNS = ("calculation", "h_box_m", "L_ext_m", "delta_eff_m", "E_kg_s", "E_top_kg_s",
           "j_projected_kg_m2_s", "j_wet_kg_m2_s", "A_box_m2", "A_phys_m2", "A_raster_m2",
           "area_rel_err", "imbalance", "discrete_residual", "L_spread_m")


def say(s=""):
    print(s, flush=True); LINES.append(s)


def row(name, out):
    return dict(calculation=name, h_box_m=out["h_box"], L_ext_m=out["L_ext"], delta_eff_m=out["delta_eff"],
                E_kg_s=out["E"], E_top_kg_s=out["E_top"], j_projected_kg_m2_s=out["j_projected"],
                j_wet_kg_m2_s=out["j_wet"], A_box_m2=out["A_box"], A_phys_m2=out.get("A_phys", out["A_box"]),
                A_raster_m2=out["A_wet"], area_rel_err=out.get("area_rel_err", 0.0),
                imbalance=out["imbalance"], discrete_residual=out["discrete_residual"], L_spread_m=out["L_spread"])


def main():
    t0 = time.time()
    air, props = g.gas_properties()
    settings = dict(T=g.T, RH=g.RH, C_brine=g.C_BRINE, pore_side=g.PORE_SIDE, L_box=g.L_BOX,
                    delta_eff=g.DELTA_EFF, a_voxel=g.A_VOXEL, h_box=g.H_BOX,
                    D_v=air.D_v, rho_g=air.rho_g, W_sat=air.W_sat, rho_D=air.rho_D,
                    a_w=props.X_wet / air.W_sat, X_wet=props.X_wet, X_top=props.X_top)
    say("## settings")
    for k, v in settings.items():
        say(f"  {k:10s} {v}")
    results = dict(settings=settings, plane={}, pore={}, pore_shifted={})
    fields = {}
    rows = []
    say("\n## the four calculations (+ supplementary shifted footprint)")
    say("  case  wet       h_box  L_ext   E [kg/s]          j_proj [kg/m2/s]  j_wet [kg/m2/s]   A_wet [m2]  area err   imbalance  residual   L_spread [um]")
    for case, h in g.H_BOX.items():
        p, fp = g.solve_plane(h)
        q, fq = g.solve_pore(h)
        s, fs = g.solve_pore(h, shift=0.5 * g.A_VOXEL)
        results["plane"][case] = p; results["pore"][case] = q; results["pore_shifted"][case] = s
        fields[f"X_plane_{case}"] = fp.X; fields[f"X_pore_{case}"] = fq.X
        fields[f"base_flux_pore_{case}"] = fq.base_flux; fields[f"foot_{case}"] = fq.base_flux > 0.0
        for name, out in ((f"plane_{case}", p), (f"pore_{case}", q), (f"pore_shifted_{case}", s)):
            rows.append(row(name, out))
            say(f"  {case}     {name[:-2]:9s} {out['h_box']*1e3:4.1f}   {out['L_ext']*1e3:5.1f}   {out['E']:.10e}  "
                f"{out['j_projected']:.10e}  {out['j_wet']:.10e}  {out['A_wet']:.3e}  {out.get('area_rel_err', 0.0):+8.4f}  "
                f"{out['imbalance']:+.1e}  {out['discrete_residual']:.1e}  {out['L_spread']*1e6:10.4f}")
    say("\n## V0.1 planar analytical reference")
    for case in g.H_BOX:
        p = results["plane"][case]
        say(f"  {case}: j_1d {p['j_1d']:.10e}  j_num {p['j_projected']:.10e}  rel err {p['rel_err_1d']:+.2e}")
    say("\n## V0.2 vapour mass conservation (E - E_top) / E_top and discrete residual ||Ax-b||/||b||")
    for kind in ("plane", "pore"):
        for case in g.H_BOX:
            o = results[kind][case]
            say(f"  {kind:5s} {case}: imbalance {o['imbalance']:+.2e}  residual {o['discrete_residual']:.2e}")
    say("\n## V0.3 domain-height independence of the single-pore result (B = 4 mm as reference)")
    qa, qb = results["pore"]["A"], results["pore"]["B"]
    comp = {}
    for key in ("E", "j_projected", "j_wet", "L_spread"):
        comp[f"{key}_rel_diff_AB"] = (qa[key] - qb[key]) / qb[key]
        say(f"  {key:12s} A {qa[key]:.10e}  B {qb[key]:.10e}  (A - B)/B {comp[f'{key}_rel_diff_AB']:+.2e}")
    for case in g.H_BOX:
        bound = float(np.exp(-2.0 * np.pi * g.H_BOX[case] / g.L_BOX))
        comp[f"top_nonuniformity_{case}"] = results["pore"][case]["top_nonuniformity"]
        comp[f"first_mode_bound_{case}"] = bound
        say(f"  lateral non-uniformity of X at the box top, {case}: {results['pore'][case]['top_nonuniformity']:.2e}"
            f"  (first lateral mode exp(-2 pi h / L_box) = {bound:.2e})")
    pa = results["plane"]["A"]
    comp["E_pore_over_E_plane"] = qa["E"] / pa["E"]
    comp["E_pore_over_E_plane_area_scaled"] = qa["E"] / (pa["E"] * qa["A_wet"] / qa["A_box"])
    say(f"  E_pore / E_plane = {comp['E_pore_over_E_plane']:.6f}  (delta_eff / (delta_eff + L_spread) = "
        f"{g.DELTA_EFF / (g.DELTA_EFF + qa['L_spread']):.6f}); per wet area the opening evaporates "
        f"{comp['E_pore_over_E_plane_area_scaled']:.4f} x the plane flux")
    say("\n## V0.4 pore-opening rasterisation (production footprint rule of coupling.interface)")
    say(f"  physical side {g.PORE_SIDE*1e6:.0f} um, area {qa['A_phys']:.4e} m2; voxel {g.A_VOXEL*1e6:.0f} um -> "
        f"{int(round(np.sqrt(qa['n_cells'])))} x {int(round(np.sqrt(qa['n_cells'])))} cells = {qa['side_raster']*1e6:.0f} um, "
        f"area {qa['A_wet']:.4e} m2, relative error {qa['area_rel_err']:+.4f}")
    sa = results["pore_shifted"]["A"]
    comp["E_rel_diff_shifted_footprint"] = (sa["E"] - qa["E"]) / qa["E"]
    say(f"  supplementary: pore centre shifted by half a voxel -> {int(round(np.sqrt(sa['n_cells'])))} x "
        f"{int(round(np.sqrt(sa['n_cells'])))} cells, area error {sa['area_rel_err']:+.4f}, E changes by "
        f"{comp['E_rel_diff_shifted_footprint']:+.2e} (the far-field resistance delta_eff dominates the "
        f"spreading resistance L_spread, so E is insensitive to the footprint area)")
    say(f"  edge enhancement of the base flux on the wet cells (max / mean): {qa['base_flux_max_over_mean']:.3f}")
    results["comparison"] = comp
    say(f"\n[{time.time() - t0:.1f} s]")
    RES.mkdir(exist_ok=True)
    (RES / "step0_results.json").write_text(json.dumps(results, indent=1, default=float))
    np.savez_compressed(RES / "step0_fields.npz", **fields)
    with (RES / "step0_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS); w.writeheader(); w.writerows(rows)
    (RES / "step0_summary.md").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
