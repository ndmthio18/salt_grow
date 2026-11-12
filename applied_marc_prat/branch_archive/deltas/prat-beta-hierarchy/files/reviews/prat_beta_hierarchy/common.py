"""Shared pieces of the beta / hierarchy review: parameters, the layered run driver on the extended
solver, analytic closure-layer formulas and the transition-duration measurement.

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; cd reviews/prat_beta_hierarchy
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from layered_reference_ext import (DAY, Layering, PratParams, ReferenceSolver, eps_new_for_final,
                                   homogeneous_steady, prat_analytic)

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; FIG = HERE / "figures"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
P0 = PratParams()
U_EVAP = 1e-2 / DAY
HYBRID = dict(evolve_eps=True, local_props=False, transport_eps_fixed=True)
CONSISTENT = dict(evolve_eps=True, local_props=True)


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


# ------------------------------------------------------------------ analytic closure layer
def d_c_crit(eps, p=P0):
    """Crystal size at which the one-scale Prat F equals the crust limit 1 - eps:
    sqrt(Da*) = (1-eps)/eps with Da* = beta^2 (1-eps)(1-C_eq) k_r d_c / (5 eps D_eff)  ->
    d_c,crit = 5 D_eff (1-eps) / (beta^2 eps (1-C_eq) k_r)."""
    return 5.0 * p.D_eff(eps) * (1 - eps) / (p.beta ** 2 * eps * (1 - p.C_eq) * p.k_r)


def closure_layer_analytic(d_sub, eps, r, p=P0):
    """First layer n with d_c(n) = 0.14 d_sub r^(n-1) < d_c,crit (inf for r = 1 and d_c(1) >= d_c,crit)."""
    dc1 = p.dc_over_dsub * d_sub; dcc = d_c_crit(eps, p)
    if dc1 < dcc:
        return 1
    if r >= 1.0:
        return np.inf
    return int(np.floor(np.log(dc1 / dcc) / np.log(1.0 / r))) + 2


def layer_bottom(n, d_sub, alpha, c_h):
    h1 = c_h * d_sub
    return h1 * (n - 1) if alpha == 1.0 else h1 * (alpha ** (n - 1) - 1.0) / (alpha - 1.0)


# ------------------------------------------------------------------ layered run on the extended solver
def layered_run(d_sub, eps_new, r=0.5, alpha=2.0, c_h=1.0, n_layers=5, p=P0, variant=None, n_band=40,
                L_over=60.0, stretch=1.03, dt_over_band=0.5, record_every=1, face_average="arithmetic",
                h_max=None, closure_check="iterate"):
    """Time-dependent layered run; returns the thinned history and the closure event."""
    variant = CONSISTENT if variant is None else variant
    lay = Layering(d_sub=d_sub, r=r, alpha=alpha, c_h=c_h)
    a1 = prat_analytic(lay.d_c(1, p), eps_new, U_EVAP, p)
    an = prat_analytic(lay.d_c(n_layers, p), eps_new, U_EVAP, p)
    s = ReferenceSolver(lay, U_EVAP, eps_new, p, L=L_over * a1["delta_C"], dxi0=an["delta_C"] / n_band,
                        stretch=stretch, face_average=face_average, closure_check=closure_check, **variant)
    h_end = layer_bottom(n_layers + 1, d_sub, alpha, c_h) if h_max is None else min(h_max, layer_bottom(n_layers + 1, d_sub, alpha, c_h))
    dt = dt_over_band * an["delta_C"] / max(an["v"], 1e-30)
    s.run(t_end=10 * h_end / max(0.2 * an["v"], 1e-30), dt=dt, h_end=h_end, record_every=record_every)
    H = s.history
    out = dict(d_sub=d_sub, eps_new=eps_new, r=r, alpha=alpha, c_h=c_h, n_layers=n_layers, beta=p.beta, N=int(s.xi.size),
               dt=dt, variant=("hybrid" if variant is HYBRID else "consistent"), face_average=face_average, closure_check=closure_check,
               t=np.array([q["t"] for q in H]), h=np.array([q["h"] for q in H]), v=np.array([q["v"] for q in H]),
               F_prat=np.array([q["F_prat"] for q in H]), F_share=np.array([q["F_num"] for q in H]),
               C_top_star=np.array([q["C_top_star"] for q in H]), eps_final=np.array([q["eps_final"] for q in H]),
               layer_top=np.array([q["layer_top"] for q in H]), picard=np.array([q["picard"] for q in H]),
               z_m_final=np.array([q["h"] for q in H]) - s.L, L=s.L, closure=s.closure_event,
               cumulative=s.cumulative_balance(), layer_bottoms=[layer_bottom(n, d_sub, alpha, c_h) for n in range(1, n_layers + 2)],
               delta_C=[prat_analytic(lay.d_c(n, p), eps_new, U_EVAP, p)["delta_C"] for n in range(1, n_layers + 2)])
    return out


def layer_plateaus(run, frac=0.8):
    """Quasi-steady F (Prat def.) and v per completed layer: mean over the last (1 - frac) of the layer."""
    out = {}
    for n in range(1, run["n_layers"] + 1):
        idx = np.where(run["layer_top"] == n)[0]
        if len(idx) < 5:
            continue
        tail = idx[int(frac * len(idx)):]
        out[n] = dict(F=float(np.mean(run["F_prat"][tail])), v=float(np.mean(run["v"][tail])),
                      C_top_star=float(np.mean(run["C_top_star"][tail])))
    return out


def transition_duration(run, n, tol=0.01):
    """Time after the interface enters layer n until v stays within tol * |jump| of the layer-n plateau;
    returned in seconds and in band-transit units delta_C(n)/v_n."""
    idx = np.where(run["layer_top"] == n)[0]
    prev = np.where(run["layer_top"] == n - 1)[0]
    if len(idx) < 10 or len(prev) < 5:
        return None
    t0 = run["t"][idx[0]]
    v_prev = run["v"][prev[-1]]; v_new = np.mean(run["v"][idx[int(0.8 * len(idx)):]])
    jump = abs(v_new - v_prev)
    if jump == 0:
        return None
    inside = np.abs(run["v"][idx] - v_new) <= tol * jump
    # first index after which it never leaves the band
    k = len(idx)
    for i in range(len(idx) - 1, -1, -1):
        if not inside[i]:
            k = i + 1
            break
        k = i
    tau = run["t"][idx[min(k, len(idx) - 1)]] - t0
    band = run["delta_C"][n - 1] / v_new
    return dict(tau_s=float(tau), tau_over_band=float(tau / band), v_prev=float(v_prev), v_new=float(v_new))
