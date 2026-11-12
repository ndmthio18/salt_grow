"""Step 2A: independent verification of the Prat (2026) single-scale growth law, one case.

    OMP_NUM_THREADS=1 $PY runs/S_surface/step2a_prat/run.py        # ~1 min

Task A  analytic law from prat2026/core.py (via saltv2.growth.prat), with Prat's rounded a_v
        coefficient 36 and the unrounded 5 / 0.14 both reported.
Task B  independent 1-D reaction-transport reference: ``layered_reference.py`` in this directory,
        a verbatim copy (sha256 5f001da3...) of reviews/prat_layered_validation/layered_reference.py at
        commit d80e791 of branch review/prat-layered-validation; variant V0 (Prat's linearisation:
        porosity evolves by eq 2 and enters the liquid balance eq 3, the transport operator keeps the
        birth porosity, a_v and D_eff at the birth values).  Eq 29 is never used to drive it.
Task C  salt, solid and liquid-volume balances of the same solution.
Writes results/step2a_results.json and results/step2a_table.csv.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from saltv2.growth.prat import prat                                            # noqa: E402  prat2026/core.py
import params as prat_params                                                   # noqa: E402  prat2026/params.py

_spec = importlib.util.spec_from_file_location("layered_reference", HERE / "layered_reference.py")
lr = importlib.util.module_from_spec(_spec); sys.modules["layered_reference"] = lr; _spec.loader.exec_module(lr)

RES = HERE / "results"
LINES = ["# run.py output (Step 2A, Prat single-scale growth law)", ""]

# ---- case values (the Prat 2026 reproduction convention, ~21 C) --------------------------------
C_EQ = 0.264          # [-, L] paper Materials and Methods
RHO_CR = 2170.0       # [kg/m^3, L] below eq 2
RHO = 1200.0          # [kg/m^3, B] back-calibrated (prat2026/params.py)
EPS_NEW = 0.43        # birth porosity [-, A] (paper quotes Dai's measured efflorescence porosity)
K_R = 2.3e-3          # [m/s, L] Naillon
D_S = 1.3e-9          # [m^2/s, L] paper text (NOT the PNM value 1.54e-9)
BETA = 2.0            # [-, L]
J_W = 1.6e-5          # prescribed water evaporation flux [kg m^-2 s^-1, A]
D_SUB = 112.5e-6      # nominal substrate grain [m, A]
DC_OVER_DSUB = 0.14   # [-, L] Wijnhorst; an empirical ratio, not a validated local mapping
D_C = DC_OVER_DSUB * D_SUB
H0 = 0.3e-3           # initial structure height [m, A]
T_END = 24 * 3600.0   # transient run length [s, N]
U_EVAP = J_W / RHO    # evaporation velocity [m/s]

for mine, theirs, name in ((C_EQ, prat_params.C_EQ, "C_EQ"), (RHO_CR, prat_params.RHO_CR, "RHO_CR"), (RHO, prat_params.RHO, "RHO"),
                           (EPS_NEW, prat_params.EPS_S, "EPS_S"), (K_R, prat_params.KR_T21, "KR_T21"), (D_S, prat_params.DS_T21, "DS_T21"),
                           (BETA, prat_params.BETA, "BETA"), (DC_OVER_DSUB, prat_params.DC_OVER_DSUB, "DC_OVER_DSUB")):
    assert mine == theirs, name
P = lr.PratParams()
assert (P.C_eq, P.rho, P.rho_cr, P.beta, P.k_r, P.D_s, P.dc_over_dsub) == (C_EQ, RHO, RHO_CR, BETA, K_R, D_S, DC_OVER_DSUB)
V0 = dict(evolve_eps=True, liquid_consumption=True, local_props=False, transport_eps_fixed=True)


def say(s=""):
    print(s, flush=True); LINES.append(s)


def analytic(av):
    """Eq 17 (a_v given), Ds_eff, delta_C, Da* (27), F (29), C*_top (31), v (29/32), partition."""
    D_eff = float(prat.Ds_eff(EPS_NEW, D_S))
    delta_C = float(prat.boundary_layer_thickness(EPS_NEW, kr=K_R, Ds=D_S, av=av))
    Da_s = float(prat.Da_star(EPS_NEW, kr=K_R, Ds=D_S, beta=BETA, av=av))
    F = float(prat.F(Da_s))
    S_tot = C_EQ * J_W / (1.0 - C_EQ)                       # salt delivered per unit area and time, saturated supply
    v = S_tot / (RHO_CR * (1.0 - EPS_NEW)) * F              # eq 29
    C_top_star = U_EVAP * C_EQ / (BETA * (1 - EPS_NEW) * (1 - C_EQ) * K_R + np.sqrt((1 - C_EQ) * av * K_R * EPS_NEW * D_eff))  # eq 31
    return dict(a_v=av, D_eff=D_eff, delta_C=delta_C, Da_star=Da_s, F=F, v=v, v_over_U=v / U_EVAP,
                C_top_star=float(C_top_star), S_tot=S_tot, M_front=F * S_tot, M_int=(1 - F) * S_tot,
                eps_final_closure=1.0 - (1.0 - EPS_NEW) / F, U_over_Uevap=float(prat.U_over_Uevap(EPS_NEW, D_SUB, K_R, D_S, BETA)))


def main():
    t_all = time.time()
    out = dict(case=dict(C_eq=C_EQ, rho_cr=RHO_CR, rho=RHO, eps_new=EPS_NEW, k_r=K_R, D_s=D_S, beta=BETA, j_w=J_W, U_evap=U_EVAP,
                         d_sub=D_SUB, d_c=D_C, h0=H0, t_end=T_END),
               reference_solver=dict(file="layered_reference.py", sha256=hashlib.sha256((HERE / "layered_reference.py").read_bytes()).hexdigest(),
                                     source="reviews/prat_layered_validation/layered_reference.py @ d80e791c1af966df63d88655e2cc6f252a92aba8 (review/prat-layered-validation)",
                                     variant="V0", switches=V0))
    # ---- Task A
    say("## A. analytic Prat law (prat2026/core.py)")
    av36 = float(prat.specific_surface(EPS_NEW, D_SUB)); av_exact = 5.0 * (1 - EPS_NEW) / D_C
    A36, Aex = analytic(av36), analytic(av_exact)
    v_core = float(prat.growth_rate(J_W, EPS_NEW, D_SUB, K_R, D_S, BETA))
    assert abs(v_core - A36["v"]) <= 1e-12 * v_core                    # core's own path reproduces eq 29 with a_v = 36(1-eps)/d_sub
    Da_h0 = float(prat.Da(H0, EPS_NEW, D_SUB, K_R, D_S)); Pe_h0 = float(prat.Pe(U_EVAP, H0, EPS_NEW, D_S))
    v_eq28 = float(prat.growth_rate(J_W, EPS_NEW, D_SUB, K_R, D_S, BETA, h=H0))
    out["A"] = dict(a_v_36=A36, a_v_exact=Aex, v_core_growth_rate=v_core, Da_h0=Da_h0, Pe_h0=Pe_h0, tanh_sqrt_Da_h0=float(np.tanh(np.sqrt(Da_h0))),
                    v_eq28_h0=v_eq28, eq28_vs_eq29=(v_eq28 - A36["v"]) / A36["v"],
                    coefficient_36_vs_exact=dict(a_v=(av36 - av_exact) / av_exact, Da_star=(A36["Da_star"] - Aex["Da_star"]) / Aex["Da_star"],
                                                 F=(A36["F"] - Aex["F"]) / Aex["F"], v=(A36["v"] - Aex["v"]) / Aex["v"]))
    for lab, A in (("a_v = 36 (1-eps)/d_sub (core, eq 17 rounded)", A36), ("a_v = 5 (1-eps)/d_c, d_c = 0.14 d_sub (unrounded)", Aex)):
        say(f"  {lab}")
        say(f"    a_v {A['a_v']:.6e} 1/m  D_eff {A['D_eff']:.6e} m2/s  delta_C {A['delta_C']*1e6:.5f} um  Da* {A['Da_star']:.6f}  F {A['F']:.7f}")
        say(f"    C*_top {A['C_top_star']:.6e}  v {A['v']:.8e} m/s = {A['v']*1e6*3600:.5f} um/h  v/U_evap {A['v_over_U']:.7f}  U/U_evap {A['U_over_Uevap']:.6f}")
        say(f"    salt delivered {A['S_tot']:.6e} kg/m2/s -> front {A['M_front']:.6e} (F), internal {A['M_int']:.6e} (1-F); eps_final closure {A['eps_final_closure']:.6f}")
    say(f"  36 vs 5/0.14: a_v {out['A']['coefficient_36_vs_exact']['a_v']:+.3e}, Da* {out['A']['coefficient_36_vs_exact']['Da_star']:+.3e}, "
        f"F {out['A']['coefficient_36_vs_exact']['F']:+.3e}, v {out['A']['coefficient_36_vs_exact']['v']:+.3e}")
    say(f"  large-Da check at h0 = {H0*1e3:.1f} mm: Da = {Da_h0:.3e}, Pe = {Pe_h0:.3e}, tanh(sqrt Da) = {out['A']['tanh_sqrt_Da_h0']:.16f}, eq 28 vs eq 29: {out['A']['eq28_vs_eq29']:+.1e}")
    say(f"  reference trajectory: h(t) = {H0*1e3:.1f} mm + {Aex['v']*1e6*3600:.5f} um/h x t (unrounded a_v; the reference solver uses this a_v)")
    # ---- Task B: steady V0 with refinement, independent of eq 29
    say("\n## B. independent reference (V0), steady travelling wave, reaction-zone refinement")
    say("  n_band  N    v/U_evap     vs eq 29     F_prat (v/v_max)  F_num (P_i/(P_i+S))  C*_top       eps_final   closure 1-(1-eps)/F_num  res salt   res solid  [s]")
    B = {}
    for nb in (20, 40, 80):
        t0 = time.time()
        r = lr.homogeneous_steady(D_SUB, EPS_NEW, U_EVAP, P, n_band=nb, **V0)
        a = r["analytic"]
        rec = dict(n_band=nb, N=r["N"], v=r["v"], v_over_U=r["v"] / U_EVAP, rel_v=(r["v"] - a["v"]) / a["v"], F_prat=r["F_prat"], F_num=r["F_num"],
                   rel_F_num=(r["F_num"] - a["F"]) / a["F"], C_top_star=r["C_top_star"], rel_C_top=(r["C_top_star"] - a["C_top_star"]) / a["C_top_star"],
                   eps_final=r["eps_final"], eps_closure_num=1 - (1 - EPS_NEW) / r["F_num"], eps_closure_analytic=1 - (1 - EPS_NEW) / a["F"],
                   P_i=r["P_i"], S_int=r["S_int"], W_i=r["W_i"], W_L=r["W_L"], res_salt=r["res_salt"], res_solid=r["res_solid"],
                   dt=r["dt"], t=r["t"], h=r["h"], wall_s=time.time() - t0, analytic=a,
                   profile=dict(xi=r["xi"].tolist(), C_star=(r["C"] - C_EQ).tolist(), eps=r["eps"].tolist()))
        prof_eq19 = lr.prat_profile(r["xi"], a["delta_C"], a["C_top_star"])
        rec["profile_max_dev_over_Ctop"] = float(np.max(np.abs((r["C"] - C_EQ) - prof_eq19)) / a["C_top_star"])
        rec["eps_min"] = float(r["eps"].min()); rec["eps_max"] = float(r["eps"].max())
        B[nb] = rec
        say(f"  {nb:4d}   {r['N']:3d}  {rec['v_over_U']:.7f}  {rec['rel_v']:+.3e}   {r['F_prat']:.7f}         {r['F_num']:.7f}            {r['C_top_star']:.5e}  "
            f"{r['eps_final']:.6f}    {rec['eps_closure_num']:.6f}               {r['res_salt']:+.1e}  {r['res_solid']:+.1e}  {rec['wall_s']:.1f}")
    assert all(abs(B[nb]["analytic"][k] - Aex[k]) <= 1e-12 * abs(Aex[k]) for nb in B for k in ("Da_star", "F", "delta_C", "C_top_star", "v"))
    out["B"] = {str(k): v for k, v in B.items()}
    r40 = B[40]
    say(f"  profile: max |C*_num - C*_top exp(-xi/delta_C)| / C*_top = {B[20]['profile_max_dev_over_Ctop']:.2e} / {r40['profile_max_dev_over_Ctop']:.2e} / {B[80]['profile_max_dev_over_Ctop']:.2e} (n_band 20/40/80)")
    say(f"  porosity range in the window (n_band 40): [{r40['eps_min']:.6f}, {r40['eps_max']:.6f}]")
    say(f"  Richardson (n_band 20/40/80, v): order {np.log((B[20]['v'] - B[40]['v']) / (B[40]['v'] - B[80]['v'])) / np.log(2):.2f}")
    # ---- Task B: transient from h0 = 0.3 mm, h(t) against h0 + v t
    say("\n## B'. transient reference from h0, 24 h, against the linear trajectory h0 + v_Prat t")
    lay = lr.Layering(d_sub=D_SUB, r=1.0, alpha=1.0, c_h=1.0)
    s = lr.ReferenceSolver(lay, U_EVAP, EPS_NEW, P, L=60 * Aex["delta_C"], n_band=40, stretch=1.03, h0=H0, **V0)
    dt = 0.5 * Aex["delta_C"] / Aex["v"]
    t0 = time.time(); hist = s.run(T_END, dt); wall = time.time() - t0
    t = np.array([h["t"] for h in hist]); h = np.array([h["h"] for h in hist]); v = np.array([h["v"] for h in hist])
    v_ss = r40["v"]
    h_lin = H0 + Aex["v"] * t
    settled = np.abs(v - v_ss) / v_ss < 1e-6
    t_trans = float(t[np.argmax(settled)]) if settled.any() else None
    cum = s.cumulative_balance()
    out["B_transient"] = dict(dt=dt, n_steps=len(hist), wall_s=wall, t=t.tolist(), h=h.tolist(), v=v.tolist(), h_linear=h_lin.tolist(),
                              v_steady_num=v_ss, t_transient_1e6=t_trans, h_end=float(h[-1]), h_end_linear=float(h_lin[-1]),
                              rel_err_h_growth_end=float(((h[-1] - H0) - (h_lin[-1] - H0)) / (h_lin[-1] - H0)),
                              rel_err_v_end=float((v[-1] - Aex["v"]) / Aex["v"]), h_offset_vs_steady_line=float(h[-1] - (H0 + v_ss * t[-1])),
                              cumulative=cum, picard_max=int(max(hh["picard"] for hh in hist)))
    say(f"  dt {dt:.1f} s, {len(hist)} steps, {wall:.1f} s; v settles to 1e-6 of the steady value at t = {t_trans:.0f} s ({t_trans / (Aex['delta_C'] / Aex['v']):.1f} delta_C/v)")
    say(f"  h(24 h) numeric {h[-1]*1e6:.4f} um, linear h0 + v_Prat t {h_lin[-1]*1e6:.4f} um: growth error {out['B_transient']['rel_err_h_growth_end']:+.2e}; "
        f"v(24 h) vs eq 29 {out['B_transient']['rel_err_v_end']:+.2e}; offset of h from the steady line h0 + v_ss t: {out['B_transient']['h_offset_vs_steady_line']*1e9:.3f} nm")
    say(f"  cumulative 24 h: salt in {cum['salt_in']:.6e} kg/m2, precipitated (front + internal) {cum['salt_out']:.6e}, acc {cum['salt_acc']:.3e}, res salt {cum['res_salt']:+.2e}, res solid {cum['res_solid']:+.2e}")
    # ---- Task C: balances of the steady solution (per unit area and time)
    say("\n## C. salt, solid and liquid-volume balances of the steady V0 solution (n_band 40)")
    rho, rho_cr = RHO, RHO_CR
    P_i, S_int, W_i, W_L = r40["P_i"], r40["S_int"], r40["W_i"], r40["W_L"]
    eps_i, eps_L = EPS_NEW, r40["eps_final"]
    v_n = r40["v"]
    salt_supplied = rho * W_L * C_EQ                        # advective, through xi = L (diffusion there ~ 0)
    C = dict(
        salt_supplied=salt_supplied, salt_front=P_i, salt_internal=S_int, salt_dissolved_acc=0.0,
        salt_residual=(salt_supplied - P_i - S_int) / salt_supplied,
        salt_delivered_eq=Aex["S_tot"], salt_front_over_total=P_i / (P_i + S_int), salt_internal_over_total=S_int / (P_i + S_int),
        F_analytic=Aex["F"], partition_vs_F=(P_i / (P_i + S_int) - Aex["F"]) / Aex["F"],
        solid_front_rate=rho_cr * (1 - eps_i) * v_n, solid_front_vs_P_i=(rho_cr * (1 - eps_i) * v_n - P_i) / P_i,
        pore_volume_lost_rate=(eps_i - eps_L) * v_n, solid_internal_volume_rate=S_int / rho_cr,
        pore_vs_solid_volume=((eps_i - eps_L) * v_n - S_int / rho_cr) / (S_int / rho_cr),
        brine_volume_consumed_rate=S_int / rho, liquid_flux_increase=W_L - W_i,
        liquid_vs_brine_consumed=((W_L - W_i) - S_int / rho) / (S_int / rho),
        new_pore_space_rate=eps_i * v_n, brine_filling_new_pores=rho * eps_i * v_n,
        evaporated_water=J_W, W_i_decomposition=dict(U_evap=U_EVAP, front_solid_term=(rho_cr / rho) * (1 - eps_i) * v_n, sum=U_EVAP + (rho_cr / rho) * (1 - eps_i) * v_n, W_i=W_i),
        water_at_interface=rho * W_i - P_i, water_vs_j=(rho * W_i - P_i - J_W) / J_W,
        eps_final=eps_L, eps_final_closure_F_num=1 - (1 - eps_i) / r40["F_num"], eps_final_closure_F_analytic=Aex["eps_final_closure"],
        eps_final_vs_closure_num=(eps_L - (1 - (1 - eps_i) / r40["F_num"])) / eps_L,
        eps_final_vs_closure_analytic=(eps_L - Aex["eps_final_closure"]) / eps_L,
        internal_share_of_evaporated_water=S_int / J_W,
        constant_density_overcount=dict(note="eq 3-4 remove brine VOLUME S/rho per precipitated salt mass S; only the salt leaves the liquid, "
                                             "so the water share (1-C_eq) S of that volume is an artefact of the constant-density closure",
                                        water_share_rate=(1 - C_EQ) * S_int, over_j=(1 - C_EQ) * S_int / J_W,
                                        over_U_evap_as_velocity=((1 - C_EQ) * S_int / rho) / U_EVAP))
    out["C"] = C
    say(f"  salt [kg/m2/s]: supplied rho W_L C_eq {salt_supplied:.6e}; front P_i {P_i:.6e}; internal S {S_int:.6e}; dissolved accumulation 0 (steady); residual {C['salt_residual']:+.2e}")
    say(f"      eq 29 delivery C_eq j/(1-C_eq) = {Aex['S_tot']:.6e} (supplied/delivered - 1 = {salt_supplied / Aex['S_tot'] - 1:+.2e})")
    say(f"      front share P_i/(P_i+S) = {C['salt_front_over_total']:.7f} vs F = {Aex['F']:.7f} ({C['partition_vs_F']:+.2e}); internal share {C['salt_internal_over_total']:.7f} vs 1-F {1 - Aex['F']:.7f}")
    say(f"  solid: front rho_cr (1-eps_new) v = {C['solid_front_rate']:.6e} = P_i ({C['solid_front_vs_P_i']:+.1e}); internal volume S/rho_cr {C['solid_internal_volume_rate']:.6e} m/s = pore volume lost (eps_new - eps_final) v {C['pore_volume_lost_rate']:.6e} ({C['pore_vs_solid_volume']:+.1e})")
    say(f"  liquid: W_L - W_i = {C['liquid_flux_increase']:.6e} m/s = brine volume consumed S/rho {C['brine_volume_consumed_rate']:.6e} ({C['liquid_vs_brine_consumed']:+.1e});")
    say(f"      W_i = U_evap + (rho_cr/rho)(1-eps_new) v = {U_EVAP:.4e} + {C['W_i_decomposition']['front_solid_term']:.4e} = {W_i:.6e}; water reaching the interface rho W_i - P_i = {C['water_at_interface']:.6e} = j ({C['water_vs_j']:+.1e})")
    say(f"      new pore space eps_new v = {C['new_pore_space_rate']:.4e} m/s, filled with brine {C['brine_filling_new_pores']:.4e} kg/m2/s (part of the liquid advected with the solid, not of W)")
    say(f"  porosity: eps_final {eps_L:.6f}; 1-(1-eps_new)/F_num {C['eps_final_closure_F_num']:.6f} ({C['eps_final_vs_closure_num']:+.1e}); 1-(1-eps_new)/F_analytic {Aex['eps_final_closure']:.6f} ({C['eps_final_vs_closure_analytic']:+.1e})")
    say(f"  constant-density closure: brine volume removed per precipitated salt S/rho; the water share (1-C_eq) S = {C['constant_density_overcount']['water_share_rate']:.3e} kg/m2/s "
        f"= {C['constant_density_overcount']['over_j']*100:.2f} % of j is liquid that physically stays (as diluted brine): an idealisation of eq 4, not a conservation error of the solver")
    say(f"\n[{time.time() - t_all:.0f} s]")
    RES.mkdir(exist_ok=True)
    (RES / "step2a_results.json").write_text(json.dumps(out, indent=1, default=float))
    with (RES / "step2a_table.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "analytic_a_v_36", "analytic_a_v_exact", "V0_n20", "V0_n40", "V0_n80"])
        for k in ("v", "v_over_U", "F", "C_top_star", "eps_final", "P_i", "S_int", "res_salt", "res_solid", "N"):
            row = [k, A36.get(k, A36.get("eps_final_closure") if k == "eps_final" else ""), Aex.get(k, Aex.get("eps_final_closure") if k == "eps_final" else "")]
            for nb in (20, 40, 80):
                row.append(B[nb].get(k, B[nb]["F_prat"] if k == "F" else ""))
            w.writerow(row)
    (RES / "step2a_summary.md").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
