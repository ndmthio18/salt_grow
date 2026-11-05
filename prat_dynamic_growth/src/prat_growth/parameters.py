from __future__ import annotations

from dataclasses import dataclass, replace
from math import sqrt


@dataclass(frozen=True)
class PhysicalParameters:
    C_eq: float = 0.264
    rho_solution: float = 1200.0
    rho_crystal: float = 2170.0
    epsilon_ref: float = 0.43
    epsilon_birth: float = 0.43
    D_molecular: float = 1.3e-9
    k_reaction: float = 2.3e-3
    beta: float = 2.0
    d_sub: float = 25e-6
    h0: float = 1e-3
    U_evap: float = 1e-2 / 86400.0

    @property
    def evaporation_mass_flux(self) -> float:
        return self.rho_solution * self.U_evap

    @property
    def a_v_ref(self) -> float:
        return 36.0 * (1.0 - self.epsilon_ref) / self.d_sub

    @property
    def D_eff_ref(self) -> float:
        return sqrt(self.epsilon_ref) * self.D_molecular

    @property
    def crystal_scale(self) -> float:
        """Paper estimate d_c/d_sub approximately 0.14."""
        return 0.14 * self.d_sub

    def with_updates(self, **kwargs: float) -> "PhysicalParameters":
        return replace(self, **kwargs)


@dataclass(frozen=True)
class NumericalParameters:
    n_cells: int = 160
    grid_power: float = 2.0
    t_end: float = 6.0 * 3600.0
    output_interval: float = 300.0
    max_step: float = 60.0
    rtol: float = 1e-7
    atol_field: float = 1e-10
    atol_height: float = 1e-12
    atol_inventory: float = 1e-10
    porosity_stop_min: float = 0.02
    porosity_stop_max: float = 0.95
    concentration_stop_min: float = 0.0
    concentration_stop_max: float = 0.95

    def with_updates(self, **kwargs: float) -> "NumericalParameters":
        return replace(self, **kwargs)
