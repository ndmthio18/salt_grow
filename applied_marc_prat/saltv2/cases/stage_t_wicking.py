"""Stage-T: pre-crystallisation water and salt transport in a reservoir-fed substrate.

Regime: evaporation-wicking.  A saturated, rigid pore network sits on a reservoir (fixed C0,
pressure datum); its surface pores lose pure water at a PRESCRIBED rate (the imposed demand,
identical in every network variant), the salt stays.  No vapour solve, no crystallisation, no
saturation variable: the run answers how the substrate's internal structure changes the surface
enrichment and the per-pore supply BEFORE the first crystal, and nothing after it.

The paired comparison (the point of the stage):

    C    full 3-D network: n x n surface pores, ``n_layers`` layers, lateral throats in every
         layer, vertical throats, reservoir throats under the bottom layer
    B2   C with its SUBSURFACE lateral throats removed (``pnm.drop_throats``) -- same pores,
         coordinates, volumes, initial state, vertical throats, reservoir throats, surface
         lateral throats and demand; the only difference is the subsurface lateral path
    B1   C with every lateral throat removed (independent columns), a control
    A    the production representation: the surface layer alone, each pore on one lumped
         reservoir throat as long as the substrate (``Lt_res`` = depth), surface lateral
         throats as in C

Physical depth, lateral size and numerical resolution are separate inputs and must not be
confused: ``n_layers`` at fixed ``DZ`` is the PHYSICAL depth (case E), ``n_side`` at fixed pitch
is the LATERAL domain (case D), ``dt`` is the only numerical resolution of a pore network
(adding layers is not grid refinement -- the pore is the medium's own scale).

Mass bookkeeping: ``mode="conservative"`` uses ``pnm.conservative_step`` (storage-resolving,
closes to round-off); ``mode="production"`` runs the production ``solve_pressure`` +
``salt_step`` pair unchanged, audited with the same ``pnm.step_residuals`` yardstick, so its
known missing storage terms are measured rather than assumed.

Values and provenance ([M] measured / stated by the source, [L] literature or library value,
[A] assumption of this project, [N] numerical):
"""
from __future__ import annotations

import time

import numpy as np

from ..materials import brine
from ..materials.brine import C_SAT, GAMMA
from ..pnm import (SUBSURFACE_LATERAL, VERTICAL, PoreLayer, conservative_step, drop_throats,
                   hydrostatic_head, net_inflow, salt_step, solve_pressure, step_residuals,
                   throat_kind)
from .vtp2014_column import C0, L_COLUMN

PORE_RADIUS = 270e-6        # inscribed radius of the QSPNM dgf pores, uniform [L] (the dgf's
                            # 270 +- 10 um spread is dropped so that the null control is exact)
PITCH = 0.7e-3              # lateral pore spacing of that dgf [L]
DZ = PITCH                  # vertical spacing = lateral: an isotropic lattice [A]
DEPTH = L_COLUMN            # 35 mm, the 2014 bead column = salt transport length [M]
N_LAYERS = int(round(DEPTH / DZ))       # 50 layers; reservoir node one DZ under the bottom one
L_S = 22e-3                 # 2014 height of the bead surface over the reservoir level [M]; the
                            # hydrostatic head of gravity assumption G2 (G1 takes H = DEPTH)
J_MEAN = 1.5854e-5          # mean evaporation flux per plane area [kg m^-2 s^-1, M]: the 2014
                            # Table 1 value (as saltv2.cases.pnm_column), applied as a demand
HOT_FACTOR = 4.0            # hotspot pore's demand / mean [A]: a test forcing of the order of the
                            # enhancement over a growing patch, not a calibrated value
STRIPE_AMPLITUDE = 0.5      # E ~ 1 + a cos(2 pi x / L_box) [A]
CV_HET = 0.25               # pore-radius coefficient of variation of the heterogeneity case [A]
SEED_HET = 1
DT = 30.0                   # base time step [s, N]; refinement in CASES
T_MAX_H = 24.0              # censoring time [h, N]
RECORD_S = 600.0            # time-series record interval [s, N]
SNAP_H = (0.25, 0.5, 1.0, 1.5)        # common snapshot times, before the first saturation [h, N]
COS_THETA = 1.0             # perfectly wetting brine on the substrate, capillary bounds [A]


def build_network(variant, n_side=5, n_layers=N_LAYERS, cv=0.0, seed=0):
    """The network of ``variant`` ('A', 'B1', 'B2', 'C').  B1 / B2 are C with throats dropped,
    so they share every pore and every remaining throat with C bit for bit."""
    if variant == "A":
        return PoreLayer(1, n_side=n_side, pore_radius=None if cv > 0 else PORE_RADIUS, cv=cv,
                         seed=seed, n_layers=1, Lt_res=n_layers * DZ)
    full = PoreLayer(1, n_side=n_side, pore_radius=None if cv > 0 else PORE_RADIUS, cv=cv,
                     seed=seed, n_layers=n_layers, dz=DZ,
                     Lt_res=DZ - 2.0 * PORE_RADIUS)   # reservoir = one more node, same throat
    kind = throat_kind(full)
    if variant == "C":
        return full
    if variant == "B2":
        return drop_throats(full, kind != SUBSURFACE_LATERAL)
    if variant == "B1":
        return drop_throats(full, kind == VERTICAL)
    raise ValueError(f"unknown variant {variant!r}")


def demand(layer, pattern):
    """Prescribed water loss per pore [kg/s]: surface pores only, mean J_MEAN * pitch^2 per pore
    whatever the pattern (the total is the same in every pattern), except ``hotspot_fixedbg``:
    the 5 x 5 hotspot's LOCAL field (hot pore 4 J_MEAN, every other pore (25 - 4) / 24 J_MEAN)
    kept on any box, so that a larger box changes only the spacing of the periodic images of
    the hotspot (case D); on 5 x 5 it is ``hotspot`` exactly."""
    s = layer.L_box / layer.n_side
    xs = layer.x[layer.surface]; ys = layer.y[layer.surface]
    n = layer.n_surface
    if pattern == "uniform":
        f = np.ones(n)
    elif pattern == "stripe":
        f = 1.0 + STRIPE_AMPLITUDE * np.cos(2.0 * np.pi * xs / layer.L_box)
    elif pattern == "hotspot":
        c = (layer.n_side // 2 + 0.5) * s
        hot = int(np.argmin((xs - c) ** 2 + (ys - c) ** 2))
        f = np.full(n, (n - HOT_FACTOR) / (n - 1)); f[hot] = HOT_FACTOR
    elif pattern == "hotspot_fixedbg":
        c = (layer.n_side // 2 + 0.5) * s
        hot = int(np.argmin((xs - c) ** 2 + (ys - c) ** 2))
        f = np.full(n, (25 - HOT_FACTOR) / 24); f[hot] = HOT_FACTOR
    else:
        raise ValueError(f"unknown pattern {pattern!r}")
    if pattern != "hotspot_fixedbg":
        f = f * n / f.sum()                           # exact total
    q = np.zeros(layer.N)
    q[layer.surface] = J_MEAN * s * s * f
    return q


def network_audit(layer):
    """Geometric depth, storage and resistances of one surface pore's column [SI]."""
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    kind = throat_kind(layer)
    s = layer.L_box / layer.n_side
    k = layer.N // layer.n_surface
    vert = kind == 2
    R_h = (np.sum(1.0 / g_lat[vert]) / layer.n_surface if vert.any() else 0.0) \
        + float(np.mean(1.0 / g_res[layer.bottom]))
    R_d = (np.sum(1.0 / d_lat[vert]) / layer.n_surface if vert.any() else 0.0) \
        + float(np.mean(1.0 / d_res[layer.bottom]))
    depth = (k - 1) * layer.dz + (layer.dz if k > 1 else float(layer.Lt_res[layer.bottom][0]))
    return dict(n_pores=int(layer.N), n_layers=int(k), n_side=int(layer.n_side),
                n_throats=int(len(layer.lat)),
                n_surface_lateral=int(np.sum(kind == 0)), n_subsurface_lateral=int(np.sum(kind == 1)),
                n_vertical=int(np.sum(vert)), depth_m=float(depth), pitch_m=float(s),
                pore_volume_per_column_m3=float(layer.V.sum() / layer.n_surface),
                porosity=float(layer.V.sum() / (layer.L_box ** 2 * depth)),
                hydraulic_resistance_per_column=float(R_h),    # [Pa s m^-3]
                diffusive_resistance_per_column=float(R_d))    # [s m^-3]


def capillary_bounds(layer):
    """Two stated upper bounds on the suction the menisci can hold [Pa], Young-Laplace with
    perfect wetting and the inscribed radius (square-duct corner flow and the Mayer-Stowe-Princen
    shape factor are ignored, an order-one correction): the surface pore opening (above it the
    meniscus would recede into the pore body -- desaturation, which the model cannot represent)
    and the narrowest throat (above it the network would be invaded)."""
    rt = np.array([t[2] for t in layer.lat] + list(layer.rt_res[layer.bottom]))
    return dict(P_open=2.0 * GAMMA * COS_THETA / layer.r[layer.surface],
                P_throat_min=float(2.0 * GAMMA * COS_THETA / rt.max()),
                P_throat_max=float(2.0 * GAMMA * COS_THETA / rt.min()))


def _surface_fluxes(layer, kind, bud):
    """Per surface pore: liquid / salt from below and net lateral [kg/s]."""
    a = np.array([t[0] for t in layer.lat], int); b = np.array([t[1] for t in layer.lat], int)
    surf = layer.surface
    out = {}
    for name in ("liquid", "salt"):
        f = bud[f"{name}_lat"]
        lateral = np.zeros(layer.N); vertical = np.zeros(layer.N)
        lat = kind != 2
        np.add.at(lateral, a[lat], -f[lat]); np.add.at(lateral, b[lat], f[lat])
        v = kind == 2
        np.add.at(vertical, a[v], -f[v]); np.add.at(vertical, b[v], f[v])
        vertical += bud[f"{name}_res"]
        out[f"{name}_from_below"] = vertical[surf]
        out[f"{name}_lateral_net"] = lateral[surf]
        out[f"{name}_reservoir"] = bud[f"{name}_res"][layer.bottom]
        out[f"{name}_subsurface_lateral_abs"] = float(np.sum(np.abs(f[kind == 1])))
        out[f"{name}_surface_lateral_abs"] = float(np.sum(np.abs(f[kind == 0])))
    return out


def run(variant="C", pattern="hotspot", n_side=5, n_layers=N_LAYERS, cv=0.0, seed=0, dt=DT,
        t_max_h=T_MAX_H, mode="conservative", record_s=RECORD_S, snap_h=SNAP_H, verbose=False):
    """One Stage-T run.  Returns (summary dict, arrays dict).  The run continues past the first
    surface saturation to t_max (or until every surface pore has crossed C_sat); everything
    after the FIRST crossing is a no-precipitation counterfactual and is labelled as such."""
    if mode not in ("conservative", "production"):
        raise ValueError("mode must be 'conservative' or 'production'")
    wall0 = time.time()
    layer = build_network(variant, n_side, n_layers, cv, seed)
    kind = throat_kind(layer)
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    q = demand(layer, pattern)
    rho_res = float(brine.density(C0))
    w = np.full(layer.N, C0)
    surf = layer.surface
    n_steps = int(round(t_max_h * 3600.0 / dt))
    rec_every = max(1, int(round(record_s / dt)))
    snap_steps = {int(round(h * 3600.0 / dt)): h for h in snap_h}
    S0 = float(np.sum(brine.salt_mass(w, layer.V))); L0 = float(np.sum(layer.V * brine.density(w)))
    cum = dict(evaporated=0.0, liquid_in=0.0, salt_in=0.0, res_salt=0.0, res_water=0.0,
               res_liquid=0.0)
    max_step_res = dict(salt=0.0, water=0.0, liquid=0.0)
    t_cross = np.full(layer.n_surface, np.nan)
    first = None
    rec = {k: [] for k in ("t", "w_surf", "var_surf", "evaporated", "liquid_in", "salt_in",
                           "res_salt", "res_water", "res_liquid", "liq_sub_abs", "salt_sub_abs",
                           "p_surf_min")}
    snaps = {}
    picard = []
    m = brine.salt_mass(w, layer.V)
    bounds = capillary_bounds(layer)
    depth = n_layers * DZ                 # surface pore centre to the reservoir node (G1 head)

    def snapshot(label, t, w_now, r, p):
        sf = _surface_fluxes(layer, kind, r["budget"])
        head_G1 = hydrostatic_head(layer, w_now, depth)
        head_G2 = hydrostatic_head(layer, w_now, min(L_S, depth))
        p_req_G1 = -p[surf] + head_G1; p_req_G2 = -p[surf] + head_G2
        snaps[label] = dict(t=t, w=w_now.copy(), p=p.copy(), demand=q[surf].copy(),
                            p_req_G1=p_req_G1, p_req_G2=p_req_G2, p_visc=-p[surf],
                            P_open=bounds["P_open"], **sf)

    for n in range(1, n_steps + 1):
        w_old = w
        if mode == "conservative":
            r = conservative_step(layer, w, C0, q, g_lat, g_res, d_lat, d_res, dt)
            w = r["w"]; p = r["p"]; rho_flux = r["rho_flux"]; Q_lat, Q_res = r["Q_lat"], r["Q_res"]
            picard.append(r["iters"])
        else:
            rho_flux = brine.density(w)
            p, Q_lat, Q_res = solve_pressure(layer, q, g_lat, g_res, rho_flux, rho_res)
            m = salt_step(layer, m, C0, Q_lat, Q_res, d_lat, d_res, dt)
            w = brine.mass_fraction_from_salt_mass(m, layer.V)
        res = step_residuals(layer, w_old, w, rho_flux, Q_lat, Q_res, d_lat, d_res, C0, q, dt)
        cum["evaporated"] += res["out"]; cum["liquid_in"] += res["liquid_in"]
        cum["salt_in"] += res["salt_in"]
        cum["res_salt"] += res["salt"]; cum["res_water"] += res["water"]
        cum["res_liquid"] += res["liquid"]
        for k in max_step_res:
            max_step_res[k] = max(max_step_res[k], abs(res[k]))
        t = n * dt
        ws, wso = w[surf], w_old[surf]
        newly = np.isnan(t_cross) & (ws >= C_SAT)
        if newly.any():
            t_cross[newly] = t - dt + dt * (C_SAT - wso[newly]) / (ws[newly] - wso[newly])
        if first is None and newly.any():
            first = dict(t_sat_s=float(np.nanmin(t_cross)), step=n,
                         pore=int(np.nanargmin(t_cross)))
            snapshot("first_saturation", t, w, res, p)
        if n in snap_steps:
            snapshot(f"t{snap_steps[n]:g}h", t, w, res, p)
        if n % rec_every == 0 or n == n_steps:
            bud = res["budget"]
            rec["t"].append(t); rec["w_surf"].append(ws.copy())
            rec["var_surf"].append(float(np.var(ws)))
            for k in ("evaporated", "liquid_in", "salt_in", "res_salt", "res_water", "res_liquid"):
                rec[k].append(cum[k])
            rec["liq_sub_abs"].append(float(np.sum(np.abs(bud["liquid_lat"][kind == 1]))))
            rec["salt_sub_abs"].append(float(np.sum(np.abs(bud["salt_lat"][kind == 1]))))
            rec["p_surf_min"].append(float(p[surf].min()))
            if verbose:
                print(f"t {t/3600:7.3f} h  w_surf max {ws.max():.6f}  var {np.var(ws):.3e}  "
                      f"res_salt {cum['res_salt']:+.2e}", flush=True)
        if np.all(np.isfinite(t_cross)):
            break
    S1 = float(np.sum(brine.salt_mass(w, layer.V))); L1 = float(np.sum(layer.V * brine.density(w)))
    t_end = n * dt
    if "end" not in snaps:
        snapshot("end", t_end, w, res, p)
    salt_gain = S1 - S0
    summary = dict(
        variant=variant, pattern=pattern, n_side=n_side, n_layers=n_layers, cv=cv, seed=seed,
        dt=dt, t_max_h=t_max_h, mode=mode, t_end_h=t_end / 3600.0, steps=n,
        wall_s=time.time() - wall0, audit=network_audit(layer),
        capillary=dict(P_open_min=float(bounds["P_open"].min()), **{k: v for k, v in bounds.items()
                                                                     if k != "P_open"}),
        censored=first is None,
        t_sat_first_h=None if first is None else first["t_sat_s"] / 3600.0,
        first_pore=None if first is None else first["pore"],
        t_cross_h=(t_cross / 3600.0).tolist(),
        w_surf_end=w[surf].tolist(),
        inventory=dict(salt_0=S0, salt_end=S1, liquid_0=L0, liquid_end=L1),
        cumulative=cum, max_step_residual=max_step_res,
        rel_res_salt=cum["res_salt"] / max(abs(salt_gain), 1e-300),
        rel_res_water=cum["res_water"] / cum["evaporated"],
        rel_res_salt_of_inventory=cum["res_salt"] / S0,
        picard_iters=dict(mean=float(np.mean(picard)), max=int(np.max(picard))) if picard else None,
    )
    arrays = {f"rec_{k}": np.asarray(v) for k, v in rec.items()}
    arrays["surface_x"] = layer.x[surf]; arrays["surface_y"] = layer.y[surf]
    arrays["t_cross_s"] = t_cross; arrays["demand"] = q[surf]
    for label, s in snaps.items():
        for k, v in s.items():
            arrays[f"snap_{label}_{k}"] = np.asarray(v)
    arrays["snap_labels"] = np.array(list(snaps))
    summary["snapshots"] = {label: _snap_stats(s) for label, s in snaps.items()}
    return summary, arrays


def _cv(x):
    x = np.asarray(x, float)
    return float(np.std(x) / abs(np.mean(x))) if np.mean(x) != 0 else float("nan")


def _snap_stats(s):
    """Scalar statistics of one snapshot for the JSON summary (the full fields go to the npz)."""
    ws = s["w"][:len(s["demand"])]
    total_E = float(np.sum(s["demand"]))
    hot = int(np.argmax(s["demand"]))
    return dict(
        t_h=s["t"] / 3600.0,
        w_surf_mean=float(np.mean(ws)), w_surf_max=float(np.max(ws)), w_surf_min=float(np.min(ws)),
        w_surf_var=float(np.var(ws)), w_surf_std=float(np.std(ws)),
        liquid_from_below_cv=_cv(s["liquid_from_below"]),
        salt_from_below_cv=_cv(s["salt_from_below"]),
        reservoir_liquid_cv=_cv(s["liquid_reservoir"]),
        reservoir_salt_cv=_cv(s["salt_reservoir"]),
        hot_pore=hot, hot_demand=float(s["demand"][hot]),
        hot_liquid_from_below=float(s["liquid_from_below"][hot]),
        hot_liquid_lateral=float(s["liquid_lateral_net"][hot]),
        hot_salt_from_below=float(s["salt_from_below"][hot]),
        hot_salt_lateral=float(s["salt_lateral_net"][hot]),
        subsurface_lateral_liquid_over_E=s["liquid_subsurface_lateral_abs"] / total_E,
        surface_lateral_liquid_over_E=s["liquid_surface_lateral_abs"] / total_E,
        subsurface_lateral_salt_abs=s["salt_subsurface_lateral_abs"],
        p_visc_max=float(np.max(s["p_visc"])),
        p_req_G1_max=float(np.max(s["p_req_G1"])), p_req_G2_max=float(np.max(s["p_req_G2"])),
        margin_open_G1=float(np.min(s["P_open"] / s["p_req_G1"])),     # bound / demand, per pore
        margin_open_G2=float(np.min(s["P_open"] / s["p_req_G2"])),
    )


# ---- the case table ---------------------------------------------------------------------------
def _case(**kw):
    base = dict(variant="C", pattern="hotspot", n_side=5, n_layers=N_LAYERS, cv=0.0, seed=0,
                dt=DT, t_max_h=T_MAX_H, mode="conservative")
    base.update(kw)
    return base


CASES = {}
for _pat in ("uniform", "stripe", "hotspot"):                        # core B2 / C pairs
    for _v in ("B2", "C"):
        CASES[f"{_v}_{_pat}"] = _case(variant=_v, pattern=_pat)
for _v in ("B2", "C"):                                                 # network heterogeneity
    CASES[f"{_v}_cv25"] = _case(variant=_v, pattern="uniform", cv=CV_HET, seed=SEED_HET)
for _v in ("A", "B1"):                                                 # controls
    for _pat in ("uniform", "hotspot"):
        CASES[f"{_v}_{_pat}"] = _case(variant=_v, pattern=_pat)
for _v in ("B2", "C"):
    for _dt in (60.0, 15.0, 7.5):                                      # time-step refinement
        CASES[f"{_v}_hotspot_dt{_dt:g}"] = _case(variant=_v, dt=_dt)
    CASES[f"{_v}_hotspot_production"] = _case(variant=_v, mode="production")   # legacy audit
    for _nl in (10, 25):                                               # E: physical depth
        CASES[f"{_v}_hotspot_L{_nl}"] = _case(variant=_v, n_layers=_nl)
    for _ns in (10, 20):                                               # D: lateral domain
        CASES[f"{_v}_hotspot_n{_ns}"] = _case(variant=_v, n_side=_ns)
CASES["C_uniform_production"] = _case(variant="C", pattern="uniform", mode="production")
CASES["column_uniform"] = _case(variant="C", pattern="uniform", n_side=1)   # 1-D null reference
# second batch (added after the first one had run; the entries above are unchanged)
for _v in ("B2", "C"):
    for _ns in (10, 20):                       # D without the background-demand confound
        CASES[f"{_v}_hotspotbg_n{_ns}"] = _case(variant=_v, pattern="hotspot_fixedbg", n_side=_ns)
    for _nl in (5, 100):                       # E: 3.5 mm and 70 mm
        CASES[f"{_v}_hotspot_L{_nl}"] = _case(variant=_v, n_layers=_nl)
    CASES[f"{_v}_cv25_dt15"] = _case(variant=_v, pattern="uniform", cv=CV_HET, seed=SEED_HET,
                                     dt=15.0)
    CASES[f"{_v}_stripe_dt15"] = _case(variant=_v, pattern="stripe", dt=15.0)
BATCH2 = [k for k in CASES if "hotspotbg" in k or k.endswith(("_L5", "_L100")) or k.endswith(
    ("cv25_dt15", "stripe_dt15"))]
# third batch: is the sign of the cv25 pair a property of one realisation?
for _v in ("B2", "C"):
    for _seed in (2, 3, 4, 5):
        CASES[f"{_v}_cv25_s{_seed}"] = _case(variant=_v, pattern="uniform", cv=CV_HET, seed=_seed)
BATCH3 = [k for k in CASES if "_cv25_s" in k]
