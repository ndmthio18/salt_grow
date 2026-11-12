"""S1 of plan/19: ONE surface pore, its column, the vapour field above it and ONE single-scale
Prat structure, with the dissolved, solid, liquid and water books closed globally.

Regime: evaporation-wicking.  The column is the Stage-T / QSPNM lattice (cubic pores of
inscribed radius 270 um, pitch 0.7 mm), 50 layers = 35 mm (the 2014 bead column, D-50) or the
10-layer prototype (D-10) whose lumped reservoir throat has the series conductance of the
unresolved 41 throats; the hydraulic head stays the physical 22 mm (G2, plan/19 section 2.3).
Evaporation comes from the existing vapour solver (``growth.vapor.solve_vapor``) with the
opening at a_w(C_pore) and the structure at a_w(C_sat) -- no prescribed demand (``feedback=
False`` keeps the Stage-T prescribed flux for the regression anchor only).  No hierarchy, no
porosity evolution, no sealing, no supply gate (S2); the capillary bound is reported only.

Structure (single scale, plan/19 section 3.2): a column over the opening of birth porosity
EPS_NEW (an input, a hypothesis [A]); the pore under it is held at C_sat (network_step
crystal / w_eq with the liquid sink), the salt X it loses is the structure's ONLY salt (no seed
injection); X = X_crystal + C_sat q_fill dt where q_fill is the brine filling the new pore space;
front   rho_cr (1 - eps_new) A dh/dt = F X_crystal / dt,   internal (1 - F) X_crystal / dt,
(the brine filling the new pores, rho eps_new A dh/dt, is drawn one step later: its water share is
in q_mass, its salt in X; the structure's liquid inventory is what was actually drawn),
F = F(Da*) of prat2026.core at (eps_new, d_sub).  With the reservoir at C_sat this reproduces
Prat eq 29 exactly; with C0 = 0.25 the column delivers less (reported as X / E).

    from saltv2.cases.surface_s1 import run, plane_check
"""
from __future__ import annotations

import time

import numpy as np

from ..coupling.interface import PoreOpenings
from ..growth.grid import GAS, SALT_WET, Grid
from ..growth.prat import prat
from ..growth.vapor import VaporProperties, solve_vapor
from ..materials import brine
from ..materials.air import moist_air
from ..materials.brine import C_SAT
from ..materials.salt import EPS_S, RHO_CR
from ..pnm import PoreLayer, hydrostatic_head, network_step, step_budget
from .stage_t_wicking import DZ, J_MEAN, PITCH, PORE_RADIUS, capillary_bounds
from .vtp2014_column import DELTA_EFF_MM, RH, L_COLUMN

# ---- case values ---------------------------------------------------------------------------
C0 = 0.25                       # reservoir and initial mass fraction [-, M]
T = brine.T_REF                 # 21 C [M]
DELTA_EFF = DELTA_EFF_MM * 1e-3     # far-field distance [m, B] (vtp2014_column)
N_LAYERS_REF = int(round(L_COLUMN / DZ))     # 50: the whole 35 mm transport length resolved
N_LAYERS_PROTO = 10             # the D-10 prototype: 7 mm resolved (Stage-T: t_sat depth-independent beyond)
H_HEAD = 22e-3                  # hydraulic head, bead surface over the reservoir level [m, M] (G2)
D_SUB = 4.0 * PORE_RADIUS       # grain size the salt closures are written in: the lattice rule [A]
EPS_NEW = EPS_S                 # birth porosity of the structure [-, A] (0.43, Dai; a hypothesis)
N_LAT = 14                      # lateral voxels across the 0.7 mm box: a = 50 um [N]
N_Z = 40                        # box height 2 mm [N]; the far field sits DELTA_EFF beyond z = 0
DT = 30.0                       # [s, N]


def build_column(n_layers):
    """n_layers pores over one lumped throat whose conductance is that of the remaining
    (N_LAYERS_REF - n_layers + 1) throats in series, so every column has the transport
    resistance of 35 mm; only the pore storage below the resolved depth is dropped."""
    n_lumped = N_LAYERS_REF - n_layers + 1
    return PoreLayer(1, n_side=1, pore_radius=PORE_RADIUS, n_layers=n_layers, dz=DZ,
                     Lt_res=n_lumped * (DZ - 2.0 * PORE_RADIUS))


def gas_side(layer):
    grid = Grid((N_LAT, N_LAT, N_Z), layer.L_box / N_LAT)
    openings = PoreOpenings(grid, layer.x[layer.surface], layer.y[layer.surface],
                            np.full(layer.n_surface, PORE_RADIUS), layer.L_box)
    air = moist_air(T, RH)
    props = VaporProperties(rho_D=air.rho_D, X_wet=air.X_salt, X_top=RH * air.W_sat)
    L_ext = DELTA_EFF - N_Z * grid.a
    if L_ext < 0.0:
        raise ValueError("the box is taller than the far-field distance")
    return grid, openings, air, props, L_ext


def plane_check():
    """Full-plane wet bottom at a_w(C_sat): the vapour solve must give the 1-D flux
    rho_D (X_salt - X_top) / delta_eff exactly (plan/19 section 4.1, check 3a)."""
    layer = build_column(1)
    grid, _, air, props, L_ext = gas_side(layer)
    state = np.full(grid.shape, GAS, int)
    mask = np.ones(grid.shape[:2], bool)
    fld = solve_vapor(grid, state, props, bottom=(mask, np.full(grid.shape[:2], air.X_salt)),
                      top_external_length=L_ext)
    j_1d = air.rho_D * (air.X_salt - props.X_top) / DELTA_EFF
    j_num = fld.M_base / fld.projected_area
    return dict(j_1d=j_1d, j_num=j_num, rel_err=(j_num - j_1d) / j_1d, imbalance=fld.imbalance,
                j_mean_2014=J_MEAN)


def run(n_layers=N_LAYERS_PROTO, c0=C0, feedback=True, growth=True, t_end_h=24.0, dt=DT,
        record_s=600.0, verbose=False):
    """One S1 run.  Returns (summary, series).  ``feedback=False``: the Stage-T prescribed water
    loss J_MEAN * pitch^2 instead of the vapour solve (regression anchor, no growth possible)."""
    layer = build_column(n_layers)
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    V = np.asarray(layer.V, float); surf = int(layer.surface[0])
    grid, openings, air, props, L_ext = gas_side(layer)
    foot = openings.foot[0]
    A_struct = float(openings.voxel_area()[0])                 # what the vapour solve sees
    F = float(prat.F(prat.Da_star(EPS_NEW, D_SUB)))
    rho_sat = float(brine.density(C_SAT)); rho_res = float(brine.density(c0))
    P_open = float(np.asarray(capillary_bounds(layer)["P_open"]).ravel()[0])
    w = np.full(layer.N, c0); m = brine.salt_mass(w, V)
    held = np.zeros(layer.N, bool)
    M_front = M_int = h = v = q_fill = M_liq = 0.0; m_seed = None; t_nuc = None
    cum = dict(E=0.0, liquid_in=0.0, salt_in=0.0, X=0.0)
    S0 = float(m.sum()); L0 = float((V * brine.density(w)).sum())
    series = {k: [] for k in ("t", "E", "E_base", "E_struct", "w_surf", "h", "v", "v_prat", "X_rate",
                              "M_front", "M_int", "M_liq", "p_req_G2", "p_req_G1", "res_salt",
                              "res_mass", "res_water")}
    n_steps = int(round(t_end_h * 3600.0 / dt)); rec_every = max(1, int(round(record_s / dt)))
    q = np.zeros(layer.N)
    state = np.full(grid.shape, GAS, int)
    wall0 = time.time(); stop = "t_end"
    for n in range(1, n_steps + 1):
        # ---- gas side: the structure's wet voxels, the opening at a_w(C_pore) ------------------
        k = int(h // grid.a)
        if k >= N_Z - 2:
            stop = "structure reached the box top"; break
        state[:] = GAS; state[foot, :k] = SALT_WET
        if feedback:
            X_open = np.array([float(brine.a_w(min(w[surf], C_SAT))) * air.W_sat])
            bmask, bval = openings.bottom_condition(X_open, blocked=state[:, :, 0] == SALT_WET)
            fld = solve_vapor(grid, state, props, bottom=(bmask, bval), top_external_length=L_ext)
            E_base, E_struct = fld.M_base, fld.M_salt
        else:
            E_base, E_struct = J_MEAN * PITCH ** 2, 0.0
        E = E_base + E_struct
        # ---- liquid side: water out, plus the water share of the brine filling new pores --------
        q[surf] = E + (1.0 - C_SAT) * q_fill
        w_old = w; rho_old = brine.density(w_old)
        step = network_step(layer, m, c0, q, g_lat, g_res, d_lat, d_res, dt, crystal=held,
                            w_eq=C_SAT, crystal_liquid_sink=True)
        m = step.m; w = brine.mass_fraction_from_salt_mass(m, V)
        X = float(step.crystallised[surf])
        bud = step_budget(layer, w, rho_old, step.Q_lat, step.Q_res, d_lat, d_res, c0, rho_res)
        cum["liquid_in"] += float(bud["liquid_res"].sum()) * dt
        cum["salt_in"] += float(bud["salt_res"].sum()) * dt
        cum["E"] += E * dt; cum["X"] += X
        # ---- structure: X is its only salt; split fill brine / crystal, then front / internal ---
        X_crystal = X - C_SAT * q_fill * dt                 # the fill brine of THIS step carried C_sat q_fill dt
        M_liq += q_fill * dt                                 # the structure's brine: what the pore delivered
        if growth and held[surf]:
            if X_crystal < 0.0:
                raise RuntimeError("the pore delivered less salt than the fill brine carries; reduce dt")
            if m_seed is None:
                m_seed = X                                   # the first crystal mass: sourced
            M_front += F * X_crystal; M_int += (1.0 - F) * X_crystal
            v = F * X_crystal / dt / (RHO_CR * (1.0 - EPS_NEW) * A_struct)
            h = M_front / (RHO_CR * (1.0 - EPS_NEW) * A_struct)
            q_fill = rho_sat * EPS_NEW * A_struct * v         # fills the new pores one step later (lag dt)
        t = n * dt
        if growth and not held[surf] and w[surf] >= C_SAT:
            held[surf] = True; t_nuc = t                      # nucleation at C_sat (S_NUC = 1)
        # ---- books (global): dissolved + solid = salt in; liquid + solid = liquid in - E ---------
        S = float(m.sum()) + C_SAT * M_liq; L = float((V * brine.density(w)).sum()) + M_liq
        solid = M_front + M_int
        res_salt = (S - S0) + solid - cum["salt_in"]
        res_mass = (L - L0) + solid - (cum["liquid_in"] - cum["E"])
        res_water = ((L - S) - (L0 - S0)) - (cum["liquid_in"] - cum["salt_in"] - cum["E"])
        if n % rec_every == 0 or n == n_steps:
            j_struct = (E_struct / A_struct) if k > 0 else (E / A_struct)
            head2 = float(hydrostatic_head(layer, w, H_HEAD)[0])
            head1 = float(hydrostatic_head(layer, w, N_LAYERS_REF * DZ)[0])
            for key, val in (("t", t), ("E", E), ("E_base", E_base), ("E_struct", E_struct),
                             ("w_surf", float(w[surf])), ("h", h), ("v", v),
                             ("v_prat", float(prat.growth_rate(j_struct, EPS_NEW, D_SUB))),
                             ("X_rate", X / dt), ("M_front", M_front), ("M_int", M_int), ("M_liq", M_liq),
                             ("p_req_G2", -float(step.p[surf]) + head2), ("p_req_G1", -float(step.p[surf]) + head1),
                             ("res_salt", res_salt), ("res_mass", res_mass), ("res_water", res_water)):
                series[key].append(val)
            if verbose:
                print(f"t {t/3600:7.3f} h  E {E:.3e}  w_s {w[surf]:.5f}  h {h*1e6:7.2f} um  "
                      f"res salt {res_salt:+.1e} mass {res_mass:+.1e} water {res_water:+.1e}", flush=True)
        if not growth and w[surf] >= C_SAT:
            stop = "first saturation"; break
    series = {k: np.array(val) for k, val in series.items()}
    summary = dict(n_layers=n_layers, c0=c0, feedback=feedback, growth=growth, dt=dt, t_end=n * dt,
                   stop=stop, t_nuc=t_nuc, m_seed=m_seed, F=F, A_struct=A_struct, A_open=(2 * PORE_RADIUS) ** 2,
                   h_end=h, M_front=M_front, M_int=M_int, cum=cum, P_open=P_open, L_box=layer.L_box,
                   res_salt=res_salt, res_mass=res_mass, res_water=res_water,
                   S_scale=float(m.sum()), L_scale=float((V * brine.density(w)).sum()),
                   wall_s=time.time() - wall0)
    return summary, series


def t_sat_of(series, w_cross=C_SAT):
    """First-saturation time [s] by linear interpolation of the recorded surface fraction."""
    t, ws = series["t"], series["w_surf"]
    k = int(np.argmax(ws >= w_cross))
    if ws[k] < w_cross or k == 0:
        return np.nan
    return float(t[k - 1] + (t[k] - t[k - 1]) * (w_cross - ws[k - 1]) / (ws[k] - ws[k - 1]))
