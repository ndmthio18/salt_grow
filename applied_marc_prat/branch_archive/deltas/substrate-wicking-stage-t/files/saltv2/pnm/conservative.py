"""OPT-IN mass-conservative transport step of the saturated pore network (Stage-T).

Regime: evaporation-wicking.  Every pore stays liquid-filled and rigid (liquid volume V_i fixed),
the reservoir holds a fixed mass fraction and the pressure datum.  Saturation is not a variable:
nothing here can drain a pore, move a meniscus or stop evaporation when supply fails.

Why a second path.  The production pair ``flow.solve_pressure`` + ``transport.salt_step`` is

    liquid:  sum rho^up Q = q                        (quasi-steady: no V d rho / dt)
    salt:    V rho(w^n) (w^{n+1} - w^n) / dt = fluxes  (non-conservative, lagged density)

so the part V w d rho of the pore salt and V d rho of the pore liquid have no flux behind them
(bounded, -0.21 % of the reservoir salt on the DuMuX column, see the execution record of the
2026-09-19 mass audit).  That is harmless for the closed loop it was built for, but a study of
how fast the surface enriches before crystallisation must close the salt and water books to
round-off, otherwise the "lateral redistribution" it looks for is of the size of the bookkeeping
error.  The production functions are not touched; this module wraps ``solve_pressure`` and
replaces ``salt_step`` only for callers that ask for it.

Equations, pore i (backward Euler, rho = ``brine.density``):

    liquid  V_i [rho(w_i^{n+1}) - rho(w_i^n)] / dt = sum_j rho^up_ij Q_ji + rho^up_res Q_res,i - q_i
    salt    V_i [rho(w_i^{n+1}) w_i^{n+1} - rho(w_i^n) w_i^n] / dt
                = sum_j [rho^up Q w^up + d rho_face (w_j - w_i)]_ji + (reservoir throat) - s_i

with the same upwind rule and face density as ``salt_step`` (mass-averaged Fick, D A / L).  The
diffusive flux carries no net liquid mass (mass-averaged frame), so it is absent from the liquid
balance.  The two are coupled through rho^{n+1}; a Picard iteration on rho^{n+1} solves them: the
storage V (rho* - rho^n) / dt enters ``solve_pressure`` as an extra liquid demand, and the salt
balance is linear in w^{n+1} once rho* is fixed.  At convergence both balances hold to the Picard
tolerance, which ``step_budget`` measures instead of assuming.

``q_out`` is the liquid MASS a pore loses [kg/s]: evaporated pure water at the surface (the salt
stays).  ``s_out`` [kg/s], the salt a pore loses to a structure, is accepted for completeness and
is zero before crystallisation.  Gravity is not in the network (as in ``solve_pressure``): p is
the viscous pressure with the reservoir at 0; the hydrostatic part of the capillary demand is
added by the caller (``hydrostatic_head``), under an assumption it must state.
"""
from __future__ import annotations

import copy

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from ..materials import brine
from ..materials.constants import G
from .flow import solve_pressure

SURFACE_LATERAL, SUBSURFACE_LATERAL, VERTICAL = 0, 1, 2
PICARD_TOL = 1e-13          # [-] max |w^{k+1} - w^k| at which the rho^{n+1} iteration stops
PICARD_MAX = 100


def throat_kind(layer):
    """Per throat of ``layer.lat``: SURFACE_LATERAL (both ends at z = 0), SUBSURFACE_LATERAL
    (both ends at the same depth below it) or VERTICAL (ends at different depths)."""
    z = np.asarray(layer.z, float)
    a = np.array([t[0] for t in layer.lat], int)
    b = np.array([t[1] for t in layer.lat], int)
    if a.size == 0:
        return np.zeros(0, int)
    same = z[a] == z[b]
    return np.where(same & (z[a] == 0.0), SURFACE_LATERAL,
                    np.where(same, SUBSURFACE_LATERAL, VERTICAL))


def drop_throats(layer, keep):
    """A copy of ``layer`` with only the throats where ``keep`` is True.  Pores, volumes,
    coordinates and reservoir throats are shared, not rebuilt, so the two networks differ in
    the dropped throats and in nothing else."""
    keep = np.asarray(keep, bool)
    if keep.shape != (len(layer.lat),):
        raise ValueError("keep must have one entry per throat")
    out = copy.copy(layer)
    out.lat = [t for t, k in zip(layer.lat, keep) if k]
    return out


def _ends(layer):
    a = np.array([t[0] for t in layer.lat], int).reshape(-1)
    b = np.array([t[1] for t in layer.lat], int).reshape(-1)
    return a, b


def step_budget(layer, w_new, rho_flux, Q_lat, Q_res, d_lat, d_res, C_res, rho_res):
    """Mass fluxes [kg/s] of one step, evaluated exactly as the step evaluated them.

    ``w_new`` is the end-of-step mass fraction (the implicit state), ``rho_flux`` the densities
    the step used for its upwind and face values (rho^{n+1} for ``conservative_step``, rho^n for
    the production ``salt_step``).  Returns per throat a -> b the liquid and salt flux, and per
    pore the liquid and salt that enter through the reservoir throat.  Summing them against the
    inventories is the conservation check; it is the same function for both paths, so the two
    are audited with one yardstick."""
    a, b = _ends(layer)
    rho = np.asarray(rho_flux, float)
    w = np.asarray(w_new, float)
    Q_lat = np.asarray(Q_lat, float); Q_res = np.asarray(Q_res, float)
    rho_up = np.where(Q_lat > 0.0, rho[a], rho[b])
    w_up = np.where(Q_lat > 0.0, w[a], w[b])
    liq_lat = rho_up * Q_lat
    rho_face = 0.5 * (rho[a] + rho[b])
    salt_lat = liq_lat * w_up - np.asarray(d_lat, float).reshape(-1) * rho_face * (w[b] - w[a])
    rho_up_res = np.where(Q_res > 0.0, rho_res, rho)
    liq_res = rho_up_res * Q_res
    w_up_res = np.where(Q_res > 0.0, C_res, w)
    salt_res = liq_res * w_up_res + np.asarray(d_res, float) * 0.5 * (rho + rho_res) * (C_res - w)
    return dict(liquid_lat=liq_lat, salt_lat=salt_lat, liquid_res=liq_res, salt_res=salt_res)


def net_inflow(layer, per_throat, per_reservoir):
    """Net inflow per pore [kg/s] from a per-throat a -> b flux and a per-pore reservoir flux."""
    a, b = _ends(layer)
    net = np.array(per_reservoir, float, copy=True)
    np.add.at(net, a, -per_throat)
    np.add.at(net, b, per_throat)
    return net


def _salt_matrix(layer, rho, Q_lat, Q_res, d_lat, d_res, rho_res, dt):
    """dt x (outflow + diffusion) operator of the salt balance at fixed densities ``rho``."""
    N = layer.N
    a, b = _ends(layer)
    qa = np.maximum(Q_lat, 0.0); qb = np.maximum(-Q_lat, 0.0)
    rf = 0.5 * (rho[a] + rho[b])
    ca = dt * (qa * rho[a] + d_lat * rf)            # what a loses per unit w_a (gained by b)
    cb = dt * (qb * rho[b] + d_lat * rf)            # what b loses per unit w_b (gained by a)
    diag = dt * (np.maximum(-Q_res, 0.0) * rho + d_res * 0.5 * (rho + rho_res))
    np.add.at(diag, a, ca)
    np.add.at(diag, b, cb)
    rows = np.concatenate([a, b]); cols = np.concatenate([b, a]); vals = np.concatenate([-cb, -ca])
    return diag, rows, cols, vals


def conservative_step(layer, w, C_res, q_out, g_lat, g_res, d_lat, d_res, dt, s_out=None,
                      tol=PICARD_TOL, max_iter=PICARD_MAX):
    """One backward-Euler step of the coupled liquid + salt balance (module docstring).

    ``w`` mass fraction per pore at t^n, ``q_out`` liquid mass lost per pore [kg/s].  Returns a
    dict with the new mass fraction ``w``, the pressure ``p`` [Pa, viscous, reservoir 0], the
    volumetric throat flows ``Q_lat`` (a -> b) and ``Q_res`` (reservoir -> pore), the densities
    ``rho_flux`` the fluxes were evaluated with, and ``iters``.  Raises if the Picard iteration
    does not reach ``tol``: a half-converged step would hide the very residual it is for."""
    V = np.asarray(layer.V, float)
    w_n = np.asarray(w, float)
    rho_n = brine.density(w_n)
    m_n = V * rho_n * w_n
    rho_res = float(brine.density(C_res))
    q_out = np.asarray(q_out, float)
    s = np.zeros(layer.N) if s_out is None else np.asarray(s_out, float)
    d_lat = np.asarray(d_lat, float); d_res = np.asarray(d_res, float)
    w_it = w_n.copy()
    for k in range(1, max_iter + 1):
        rho_it = brine.density(w_it)
        storage = V * (rho_it - rho_n) / dt                       # liquid mass the pore gains
        p, Q_lat, Q_res = solve_pressure(layer, q_out + storage, g_lat, g_res, rho_it, rho_res)
        diag, rows, cols, vals = _salt_matrix(layer, rho_it, Q_lat, Q_res, d_lat, d_res,
                                              rho_res, dt)
        diag = diag + V * rho_it
        rhs = (m_n - dt * s
               + dt * (np.maximum(Q_res, 0.0) * rho_res + d_res * 0.5 * (rho_it + rho_res)) * C_res)
        idx = np.arange(layer.N)
        A = coo_matrix((np.concatenate([vals, diag]),
                        (np.concatenate([rows, idx]), np.concatenate([cols, idx]))),
                       shape=(layer.N, layer.N)).tocsr()
        w_new = spsolve(A, rhs)
        change = float(np.max(np.abs(w_new - w_it)))
        w_it = w_new
        if change <= tol:
            break
    else:
        raise RuntimeError(f"conservative_step: rho^(n+1) Picard iteration did not settle in "
                           f"{max_iter} iterations (last change {change:.3e} > {tol:.0e})")
    if np.any(w_it < 0.0):
        raise RuntimeError("conservative_step: negative mass fraction; reduce dt")
    return dict(w=w_it, p=p, Q_lat=Q_lat, Q_res=Q_res, rho_flux=rho_it, iters=k)


def step_residuals(layer, w_n, w_new, rho_flux, Q_lat, Q_res, d_lat, d_res, C_res, q_out, dt,
                   s_out=None):
    """Global liquid / salt / water residuals of one step [kg], inventory change minus what the
    fluxes delivered, for EITHER path (pass the densities that path used for its fluxes).

    Inventories are the true ones, V rho(w) and V rho(w) w at both ends of the step; the
    production path's missing storage terms therefore show up here, not in its own algebra."""
    V = np.asarray(layer.V, float)
    rho_res = float(brine.density(C_res))
    bud = step_budget(layer, w_new, rho_flux, Q_lat, Q_res, d_lat, d_res, C_res, rho_res)
    s = 0.0 if s_out is None else float(np.sum(s_out))
    L0 = float(np.sum(V * brine.density(w_n))); L1 = float(np.sum(V * brine.density(w_new)))
    S0 = float(np.sum(brine.salt_mass(w_n, V))); S1 = float(np.sum(brine.salt_mass(w_new, V)))
    liq_in = float(np.sum(bud["liquid_res"])) * dt
    salt_in = float(np.sum(bud["salt_res"])) * dt
    out = float(np.sum(q_out)) * dt
    return dict(liquid=(L1 - L0) - (liq_in - out - s * dt),
                salt=(S1 - S0) - (salt_in - s * dt),
                water=((L1 - S1) - (L0 - S0)) - ((liq_in - salt_in) - out),
                liquid_in=liq_in, salt_in=salt_in, out=out, budget=bud)


def hydrostatic_head(layer, w, H):
    """Hydrostatic pressure [Pa] of the liquid between each SURFACE pore and a reservoir free
    surface a height ``H`` [m] below it, g * integral_0^H rho dz along the pore's own column.

    The column is the stack of pores under the surface pore (same x, y), layer j owning the
    depth interval [j dz, (j + 1) dz); the deepest pore also owns everything below the resolved
    column (for a one-layer network: the whole lumped reservoir throat, at the SURFACE pore's
    density, which overstates the head).  Gravity is not in the flow solve, so this is an
    a-posteriori demand, not a coupled term: it does not redistribute the flow (the gravity-driven
    solutal instability of a dense brine over a lighter one is absent)."""
    rho = brine.density(np.asarray(w, float))
    n_s = layer.n_surface
    k = layer.N // n_s
    col = rho.reshape(k, n_s)                       # layer by layer, surface first
    dz = float(layer.dz) if k > 1 else 0.0
    top = np.arange(k) * dz
    bottom = np.append(top[1:], np.inf)
    span = np.clip(np.minimum(bottom, H) - top, 0.0, None)
    return G * (span[:, None] * col).sum(axis=0)
