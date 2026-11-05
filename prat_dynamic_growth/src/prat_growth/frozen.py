"""Instantaneous quasi-steady concentration controls; structure is a coefficient.

Freezing a snapshot does NOT set epsilon_t=0. The solid-volume equation still
defines dE/dt, and L_left-L_right = rho*dE/dt + R*volume. Subtracting C times
this liquid balance from salt gives the nonlinear sink (1-C)R. Only dC/dt in
ALE coordinates is set to zero. Q0/Q1 omit concentration advection, not liquid
mass conservation; they are reduced operators, not stationary solid states.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import root
from .analytic import analytical_profile
from .grid import top_refined_grid


@dataclass
class FrozenResult:
    case: str
    concentration: np.ndarray
    c_top: float
    velocity: float
    reaction: np.ndarray
    epsilon: np.ndarray
    d_eff: np.ndarray
    a_v: np.ndarray
    liquid_flux: np.ndarray
    pore_rate: np.ndarray
    residual_inf: float
    interface_residual: float


def solve_frozen(p, n, height, epsilon, case="Q6", *, linear_reaction=None, linear_boundary=None):
    if case not in [f"Q{i}" for i in range(7)]:
        raise ValueError(case)
    f, x, dx = top_refined_grid(n.n_cells, n.grid_power)
    eps = np.full_like(epsilon, p.epsilon_ref) if case in ("Q0", "Q1", "Q2") else epsilon.copy()
    av = 36*(1-eps)/p.d_sub if case in ("Q5", "Q6") else np.full_like(eps, p.a_v_ref)
    deff = np.sqrt(eps)*p.D_molecular if case in ("Q4", "Q6") else np.full_like(eps, p.D_eff_ref)
    rho, rhos, ce, j = p.rho_solution, p.rho_crystal, p.C_eq, p.evaporation_mass_flux
    vol = height*dx
    k = rho*eps*deff
    conduct = np.empty(n.n_cells+1)
    conduct[0] = k[0]/(height*x[0])
    conduct[-1] = k[-1]/(height*(1-x[-1]))
    conduct[1:-1] = 1/(height*((f[1:-1]-x[:-1])/k[:-1]+(x[1:]-f[1:-1])/k[1:]))
    scale_c = 1e-5
    scale = (conduct[:-1]+conduct[1:]+rho*av*p.k_reaction*vol)*scale_c
    top_scale = max(j*ce, 1e-30)
    linear_reaction = case == "Q0" if linear_reaction is None else linear_reaction
    linear_boundary = case == "Q0" if linear_boundary is None else linear_boundary

    def evaluate(u, details=False):
        s, st = u[:-1]*scale_c, u[-1]*scale_c
        c, ct = ce+s, ce+st
        v = p.beta*rho/rhos*p.k_reaction*st
        jf = rhos*(1-p.epsilon_birth)*v
        reaction = rho*av*p.k_reaction*s
        diff = -conduct*np.diff(np.r_[0., s, st])
        F = np.empty(n.n_cells+1)
        F[-1] = j+(rho*p.epsilon_birth+rhos*(1-p.epsilon_birth))*v
        F[:-1] = F[-1]+np.cumsum(((1-rho/rhos)*reaction*vol)[::-1])[::-1]
        em = np.r_[0., -f[1:-1]*v*np.where(v >= 0, eps[1:], eps[:-1]), -v*p.epsilon_birth]
        L = F+rho*em
        erate = em[:-1]-em[1:]-reaction*vol/rhos
        sink = (1-ce if linear_reaction else 1-c)*reaction*vol
        residual = diff[:-1]-diff[1:]-sink
        if case not in ("Q0", "Q1"):
            sf = np.r_[0., np.where(L[1:-1] >= 0, s[:-1], s[1:]), st]
            residual += L[:-1]*(sf[:-1]-s)-L[1:]*(sf[1:]-s)
        boundary_c = ce if linear_boundary else ct
        boundary = diff[-1] - ((1-boundary_c)*jf-boundary_c*j)
        scaled = np.r_[residual/scale, boundary/top_scale]
        if details:
            return c, ct, v, reaction, F, erate, scaled
        return scaled

    cguess, ctguess = analytical_profile(p, height, x)
    guess = np.r_[cguess-ce, ctguess-ce]/scale_c
    solution = root(evaluate, guess, method="hybr", options={"xtol": 1e-10})
    c, ct, v, reaction, F, erate, residual = evaluate(solution.x, True)
    if np.max(np.abs(residual)) > 1e-8 or not np.all((c >= 0) & (c <= 1)) or not 0 <= ct <= 1:
        raise RuntimeError(f"{case}: {solution.message}; residual={max(abs(residual))}")
    return FrozenResult(case, c, ct, v, reaction, eps, deff, av, F, erate,
                        float(max(abs(residual[:-1]))), float(abs(residual[-1])))
