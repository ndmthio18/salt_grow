from __future__ import annotations

from dataclasses import dataclass
from math import sqrt, tanh

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from .grid import top_refined_grid
from .parameters import PhysicalParameters


@dataclass(frozen=True)
class DimensionlessGroups:
    Pe: float
    Da: float
    Da_top: float
    Da_star: float
    delta_C: float


def dimensionless_groups(p: PhysicalParameters, height: float, epsilon: float | None = None,
                         a_v: float | None = None, D_eff: float | None = None) -> DimensionlessGroups:
    eps = p.epsilon_ref if epsilon is None else epsilon
    av = p.a_v_ref if a_v is None else a_v
    deff = p.D_eff_ref if D_eff is None else D_eff
    one_minus_c = 1.0 - p.C_eq
    pe = p.U_evap * height / (eps * deff)
    da = height**2 * one_minus_c * av * p.k_reaction / (eps * deff)
    da_top = p.beta * (1.0 - eps) * one_minus_c * height * p.k_reaction / (eps * deff)
    da_star = (
        p.beta**2 * (1.0 - eps) ** 2 * one_minus_c * p.k_reaction
        / (eps * deff * av)
    )
    delta = sqrt(eps * deff / (av * p.k_reaction * one_minus_c))
    return DimensionlessGroups(pe, da, da_top, da_star, delta)


def general_growth_velocity(p: PhysicalParameters, height: float) -> float:
    g = dimensionless_groups(p, height)
    factor = sqrt(g.Da_star) * tanh(sqrt(g.Da))
    prefactor = (
        p.C_eq * p.evaporation_mass_flux
        / (p.rho_crystal * (1.0 - p.epsilon_ref) * (1.0 - p.C_eq))
    )
    return prefactor * factor / (1.0 + factor)


def large_da_growth_velocity(p: PhysicalParameters) -> float:
    g = dimensionless_groups(p, p.h0)
    factor = sqrt(g.Da_star)
    prefactor = (
        p.C_eq * p.evaporation_mass_flux
        / (p.rho_crystal * (1.0 - p.epsilon_ref) * (1.0 - p.C_eq))
    )
    return prefactor * factor / (1.0 + factor)


def analytical_profile(p: PhysicalParameters, height: float, x: np.ndarray) -> tuple[np.ndarray, float]:
    """Paper Eq. (19) and Robin condition, returned as mass fraction C."""
    g = dimensionless_groups(p, height)
    lam = sqrt(g.Da)
    alpha = lam / tanh(lam)
    cstar_top = g.Pe * p.C_eq / (g.Da_top + alpha)
    # Stable sinh(lam*x)/sinh(lam), including x=0 and very large lam.
    ratio = np.exp(lam * (x - 1.0)) * (-np.expm1(-2.0 * lam * x)) / (-np.expm1(-2.0 * lam))
    return p.C_eq + cstar_top * ratio, p.C_eq + cstar_top


@dataclass(frozen=True)
class StageAResult:
    n_cells: int
    x: np.ndarray
    concentration: np.ndarray
    c_top: float
    velocity: float
    velocity_exact: float
    velocity_large_da: float
    relative_velocity_error: float
    max_profile_error: float
    groups: DimensionlessGroups


def solve_stage_a_finite_volume(p: PhysicalParameters, height: float, n_cells: int,
                                grid_power: float = 2.0) -> StageAResult:
    """Independent FV solve of paper Eq. (18) with Eqs. (14), (22)."""
    faces, centers, widths = top_refined_grid(n_cells, grid_power)
    g = dimensionless_groups(p, height)
    matrix = lil_matrix((n_cells, n_cells), dtype=float)
    rhs = np.zeros(n_cells)

    for i in range(n_cells):
        # Right-face gradient contribution.
        if i < n_cells - 1:
            distance = centers[i + 1] - centers[i]
            matrix[i, i] += -1.0 / distance
            matrix[i, i + 1] += 1.0 / distance
        else:
            distance = faces[-1] - centers[i]
            denom = 1.0 + distance * g.Da_top
            matrix[i, i] += -g.Da_top / denom
            rhs[i] += -g.Pe * p.C_eq / denom

        # Minus left-face gradient contribution.
        if i > 0:
            distance = centers[i] - centers[i - 1]
            matrix[i, i] += -1.0 / distance
            matrix[i, i - 1] += 1.0 / distance
        else:
            distance = centers[i] - faces[0]
            matrix[i, i] += -1.0 / distance

        matrix[i, i] += -g.Da * widths[i]

    cstar = spsolve(matrix.tocsr(), rhs)
    top_distance = faces[-1] - centers[-1]
    top_gradient = (g.Pe * p.C_eq - g.Da_top * cstar[-1]) / (1.0 + top_distance * g.Da_top)
    cstar_top = cstar[-1] + top_distance * top_gradient
    c_top = p.C_eq + cstar_top
    velocity = p.beta * p.rho_solution / p.rho_crystal * p.k_reaction * cstar_top
    exact_velocity = general_growth_velocity(p, height)
    exact_profile, _ = analytical_profile(p, height, centers)
    profile_error = float(np.max(np.abs((p.C_eq + cstar) - exact_profile)))
    relative_error = abs(velocity - exact_velocity) / abs(exact_velocity)
    return StageAResult(
        n_cells=n_cells,
        x=centers,
        concentration=p.C_eq + cstar,
        c_top=c_top,
        velocity=velocity,
        velocity_exact=exact_velocity,
        velocity_large_da=large_da_growth_velocity(p),
        relative_velocity_error=relative_error,
        max_profile_error=profile_error,
        groups=g,
    )
