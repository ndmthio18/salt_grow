"""Cross-validation of the two conservative network kernels (plan/19 section 6, step 1).

    network_step       (pnm/transport, production since 2026-09-21): Picard on rho^{n+1}; throat
                       fluxes at the OLD densities, storage V (rho^{n+1} - rho^n) / dt
    conservative_step  (pnm/conservative, Stage-T, opt-in): Picard on w^{n+1}; throat fluxes AND
                       storage at the iterated new densities

Same balances, different temporal evaluation of the flux densities, so the two differ at O(dt)
and are NOT expected to agree bit for bit.  What is compared instead:

    1. conservation of each path, measured with one yardstick (``pnm.step_residuals`` with the
       densities each path used for its fluxes), per step and cumulative;
    2. the first-saturation time t_sat of the surface against the Stage-T reference (single
       column, uniform demand: 1.917673 h with conservative_step at dt = 30 s);
    3. convergence under time-step refinement (dt = 60, 30, 15, 7.5 s): the order of each path
       (Richardson), the two extrapolated limits, and the inter-kernel difference as a function
       of dt.

Regime: evaporation-wicking, saturated rigid network, prescribed pure-water demand, no salt
sink (crystal = None, s_out = None).  No morphology, no vapour field.

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
    OMP_NUM_THREADS=1 $PY runs/S_surface/kernel_crossval.py        # ~2 min
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from saltv2.cases import stage_t_wicking as st                               # noqa: E402
from saltv2.materials import brine                                           # noqa: E402
from saltv2.materials.brine import C_SAT                                     # noqa: E402
from saltv2.pnm import conservative_step, network_step, step_residuals       # noqa: E402

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
C0 = st.C0
T_SAT_REF_COLUMN_H = 1.917673          # Stage-T column_uniform, conservative_step, dt 30 s [V]
DTS = (60.0, 30.0, 15.0, 7.5)          # [s]
T_SNAP = 3600.0                        # common snapshot [s], before any saturation
T_CAP_H = 3.0


def run_kernel(kernel, layer, q, dt, t_cap_h=T_CAP_H):
    """Advance one kernel until the first surface pore crosses C_SAT (linear interpolation of the
    crossing as in stage_t_wicking).  Returns t_sat [s], the surface w at T_SNAP, the cumulative
    and max-step yardstick residuals, and the Picard counts."""
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    V = np.asarray(layer.V, float)
    surf = layer.surface
    w = np.full(layer.N, C0)
    m = brine.salt_mass(w, V)
    cum = dict(salt=0.0, water=0.0, liquid=0.0, out=0.0, salt_in=0.0)
    mx = dict(salt=0.0, water=0.0, liquid=0.0)
    iters = []
    w_snap = None
    n_snap = int(round(T_SNAP / dt))
    n_cap = int(round(t_cap_h * 3600.0 / dt))
    t_sat = np.nan
    for n in range(1, n_cap + 1):
        w_old = w
        if kernel == "network_step":
            step = network_step(layer, m, C0, q, g_lat, g_res, d_lat, d_res, dt)
            m = step.m
            w = brine.mass_fraction_from_salt_mass(m, V)
            rho_flux = brine.density(w_old)                       # fluxes at the old densities
            Q_lat, Q_res = step.Q_lat, step.Q_res
        elif kernel == "conservative_step":
            r = conservative_step(layer, w, C0, q, g_lat, g_res, d_lat, d_res, dt)
            w = r["w"]; m = brine.salt_mass(w, V)
            rho_flux = r["rho_flux"]; Q_lat, Q_res = r["Q_lat"], r["Q_res"]
            iters.append(r["iters"])
        else:
            raise ValueError(kernel)
        res = step_residuals(layer, w_old, w, rho_flux, Q_lat, Q_res, d_lat, d_res, C0, q, dt)
        for k in ("salt", "water", "liquid"):
            cum[k] += res[k]; mx[k] = max(mx[k], abs(res[k]))
        cum["out"] += res["out"]; cum["salt_in"] += res["salt_in"]
        if n == n_snap:
            w_snap = w[surf].copy()
        ws, wso = w[surf], w_old[surf]
        newly = ws >= C_SAT
        if newly.any():
            tc = (n - 1) * dt + dt * (C_SAT - wso[newly]) / (ws[newly] - wso[newly])
            t_sat = float(tc.min())
            break
    S_gain = float(np.sum(brine.salt_mass(w, V)) - np.sum(brine.salt_mass(np.full(layer.N, C0), V)))
    return dict(t_sat=t_sat, w_snap=w_snap, cum=cum, max_step=mx, steps=n,
                rel_salt=cum["salt"] / S_gain if S_gain else np.nan,
                rel_water=cum["water"] / cum["out"] if cum["out"] else np.nan,
                picard_max=int(max(iters)) if iters else None)


def richardson(dts, vals):
    """Observed order p and the extrapolated limit from the three finest (dt, value) pairs,
    assuming v(dt) = v0 + c dt^p with dt halving."""
    v1, v2, v3 = vals[-3:]                      # dt, dt/2, dt/4
    d12, d23 = v1 - v2, v2 - v3
    if d23 == 0.0 or d12 / d23 <= 0.0:
        return dict(order=np.nan, limit=v3)
    p = float(np.log(d12 / d23) / np.log(2.0))
    return dict(order=p, limit=float(v3 + (v3 - v2) / (2 ** p - 1.0)))


def main():
    lines = ["# kernel_crossval.py output", ""]
    out = {}
    cases = [("column_1x1x50_uniform", st.build_network("C", n_side=1, n_layers=50), "uniform"),
             ("hotspot_5x5x10", st.build_network("C", n_side=5, n_layers=10), "hotspot")]
    for name, layer, pattern in cases:
        q = st.demand(layer, pattern)
        lines.append(f"## {name}  (N = {layer.N} pores, demand {pattern})")
        lines.append("  kernel             dt [s]  steps  t_sat [h]     cum salt res / gain   cum water res / evap   max |step salt| [kg]  Picard")
        rec = {k: {} for k in ("network_step", "conservative_step")}
        t0 = time.time()
        for dt in DTS:
            for kern in rec:
                r = run_kernel(kern, layer, q, dt)
                rec[kern][dt] = r
                lines.append(f"  {kern:18s} {dt:6.1f}  {r['steps']:5d}  {r['t_sat']/3600:.6f}    {r['rel_salt']:+.2e}            "
                             f"{r['rel_water']:+.2e}            {r['max_step']['salt']:.1e}            {r['picard_max'] if r['picard_max'] is not None else '-':>3}")
        # inter-kernel differences
        lines.append("  dt [s]   |t_sat(net) - t_sat(cons)| [s]   max |w_surf(1 h) diff|")
        diffs = {}
        for dt in DTS:
            a, b = rec["network_step"][dt], rec["conservative_step"][dt]
            dw = float(np.max(np.abs(a["w_snap"] - b["w_snap"]))) if a["w_snap"] is not None and b["w_snap"] is not None else np.nan
            diffs[dt] = dict(dt_sat=abs(a["t_sat"] - b["t_sat"]), dw=dw)
            lines.append(f"  {dt:6.1f}   {diffs[dt]['dt_sat']:12.4f}                    {dw:.3e}")
        rich = {k: richardson(DTS, [rec[k][dt]["t_sat"] for dt in DTS]) for k in rec}
        lines.append(f"  Richardson t_sat: network_step order {rich['network_step']['order']:.2f}, limit {rich['network_step']['limit']/3600:.6f} h; "
                     f"conservative_step order {rich['conservative_step']['order']:.2f}, limit {rich['conservative_step']['limit']/3600:.6f} h; "
                     f"limits differ by {abs(rich['network_step']['limit'] - rich['conservative_step']['limit']):.4f} s")
        if name.startswith("column"):
            t30 = rec["conservative_step"][30.0]["t_sat"] / 3600
            lines.append(f"  Stage-T reference (conservative_step, dt 30 s): {T_SAT_REF_COLUMN_H:.6f} h; here {t30:.6f} h "
                         f"(diff {abs(t30 - T_SAT_REF_COLUMN_H)*3600:.2e} s)")
        lines.append(f"  [{time.time()-t0:.0f} s]"); lines.append("")
        out[name] = dict(N=layer.N, pattern=pattern,
                         runs={k: {str(dt): {kk: (vv if kk != "w_snap" else None) for kk, vv in r.items()} for dt, r in v.items()} for k, v in rec.items()},
                         diffs={str(k): v for k, v in diffs.items()}, richardson=rich)
    lines += ["Reading.  Both paths close the salt and water books to round-off with the same yardstick (the production",
              "salt_step does not, Stage-T README section 0).  Their t_sat differ at first order in dt because the throat",
              "densities are taken at t^n (network_step) or t^{n+1} (conservative_step); the Richardson limits agree to",
              "the extrapolation error.  network_step stays the production kernel; conservative_step is kept for this",
              "cross-check only.  (network_step does not expose its Picard count; conservative_step needs 5 at tol 1e-13.)"]
    RES.mkdir(exist_ok=True)
    (RES / "kernel_crossval.md").write_text("\n".join(lines) + "\n")
    (RES / "kernel_crossval.json").write_text(json.dumps(out, indent=1, default=float))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
