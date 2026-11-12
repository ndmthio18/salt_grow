"""Moist air over the brine: the gas side of the evaporation.

Two ways to get the three numbers the vapour problem needs (D_v, W_sat, rho_g):

* ``moist_air(T, rh)`` -- functions of the temperature the CASE declares, ported expression by
  expression from the DuMuX material library the reference case runs on (file and line at each
  function; tree: dumux_theresa/dumux/dumux/dumux/material/).
* the constants ``D_V``, ``RHO_G``, ``W_SAT`` (``LEGACY``), taken from QSPNM's shared_params and
  used by every run record to 2026-09-20.  They ARE those functions at one state, which the
  tests pin: T = 293.15 K, p = 1e5 Pa, and for the density the mean of the gas at a PURE-WATER
  surface and at RH = 0.00293 / 0.01468 (QSPNM's default environment) -- the averaging rule of
  the reference case's boundary condition.  So they are mutually consistent, at 20 C; what they
  are not consistent with is the liquid side (``brine.T_REF`` = 21 C, Pitzer parameters 25 C)
  and the 22 C of the 2014 experiments.  (Until 2026-09-21 this docstring said, from a
  back-calculation against other correlations, that W_SAT was a 20.2 C and RHO_G a 21.6 C
  value; that was wrong.)

The Dirichlet value at a wet surface stays ``a_w W_sat``, linear in the activity, as the vapour
solve has always taken it; the exact mass fraction at mole fraction a_w x_sat differs from it
by about (1 - a_w)(1 - M_w / M_air) x_sat relative, 0.2 % at saturation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .brine import C_SAT, P_REF, a_w
from .constants import M_H2O

M_AIR = 0.02896             # molar mass of air [kg/mol, L] DuMuX components/air.hh:50-51
R_IDEAL_GAS = 8.314472      # [J mol^-1 K^-1, L] DuMuX constants.hh:32, the value ITS ideal gas
                            # uses; materials.constants.R_GAS = 8.314 is the project's own and
                            # would move the density in the fifth digit

D_V = 2.4189e-5             # free-air vapour diffusivity [m^2/s, L] (QSPNM shared_params)
RHO_G = 1.18186             # moist-air density [kg/m^3, L] (QSPNM)
W_SAT = 0.01468             # saturated vapour mass fraction [-, L] (QSPNM)

# derived, computed once
RHO_D = RHO_G * D_V         # [kg m^-1 s^-1]: vapour flux = RHO_D grad X
X_SALT = a_w(C_SAT) * W_SAT  # vapour mass fraction at a wet salt surface (Dirichlet value)

T_LEGACY = 293.15           # [K] the state the three constants above belong to (see docstring)
RH_LEGACY = 0.00293 / 0.01468


def saturation_pressure(T):
    """Vapour pressure of pure water [Pa] at T [K] -- IAPWS-IF97 region 4, verbatim from DuMuX
    components/iapws/region4.hh:52-73 (``Region4::saturationPressure``)."""
    n = (0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
         0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
         -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
         0.65017534844798e3)
    T = np.asarray(T, float)
    sigma = T + n[8] / (T - n[9])
    A = (sigma + n[0]) * sigma + n[1]
    B = (n[2] * sigma + n[3]) * sigma + n[4]
    C = (n[5] * sigma + n[6]) * sigma + n[7]
    tmp = 2.0 * C / (np.sqrt(B * B - 4.0 * A * C) - B)
    return 1e6 * tmp ** 4


def vapour_diffusivity(T, p=P_REF):
    """Binary diffusion coefficient of water vapour in air [m^2/s] -- DuMuX
    binarycoefficients/h2o_air.hh:56-68 (``H2O_Air::gasDiffCoeff``, after Walker et al. 1981):
    2.13e-5 (1e5 / p) (T / 273.15)^1.8."""
    return 2.13e-5 * (1.0e5 / p) * (np.asarray(T, float) / 273.15) ** 1.8


def vapour_mass_fraction(x):
    """Mass fraction of water vapour in air at its MOLE fraction x."""
    x = np.asarray(x, float)
    return x * M_H2O / (x * M_H2O + (1.0 - x) * M_AIR)


def gas_density(T, x_vapour, p=P_REF):
    """Density of moist air [kg/m^3] as an ideal mixture at vapour mole fraction ``x_vapour``:
    p (x M_w + (1 - x) M_air) / (R T) -- DuMuX idealgas.hh:40,60 with components/air.hh:73-77,
    the form its H2OAir fluid system takes with ideal-gas components."""
    x = np.asarray(x_vapour, float)
    return p * (x * M_H2O + (1.0 - x) * M_AIR) / (R_IDEAL_GAS * np.asarray(T, float))


@dataclass(frozen=True)
class MoistAir:
    """The gas-side numbers of one run; ``T`` is None for the legacy constants."""
    T: float | None
    D_v: float              # [m^2/s]
    W_sat: float            # saturated vapour mass fraction over pure water [-]
    rho_g: float            # [kg/m^3]

    @property
    def rho_D(self):
        """[kg m^-1 s^-1]: vapour flux = rho_D grad X."""
        return self.rho_g * self.D_v

    @property
    def X_salt(self):
        """Vapour mass fraction at a wet salt surface (the Dirichlet value there)."""
        return a_w(C_SAT) * self.W_sat


LEGACY = MoistAir(T=None, D_v=D_V, W_sat=W_SAT, rho_g=RHO_G)


def moist_air(T, rh, a_w_surface=None, p=P_REF):
    """The gas side at the case's temperature ``T`` [K] and far-field humidity ``rh``.

    D_v and W_sat are the DuMuX expressions at (T, p).  rho_g is one number for the whole gas
    domain, as the vapour solve needs it (rho_g D_v is constant there): the mean of the moist-air
    density at the evaporating surface and in the far field, which is how the reference case's
    boundary condition averages it (pnm2pncmin2 problem.hh, ``diffWaterFlux``: gasDensityAvg =
    0.5 (pore gas + free-flow reference)).  The surface state is a wet salt surface,
    ``a_w_surface`` = a_w(C_sat), unless the caller states another [A]: bare openings below
    saturation are fresher, which moves rho_g by less than 0.1 %.
    """
    a_surface = float(a_w(C_SAT)) if a_w_surface is None else float(a_w_surface)
    x_sat = float(saturation_pressure(T)) / p
    rho_g = 0.5 * (float(gas_density(T, a_surface * x_sat, p))
                   + float(gas_density(T, rh * x_sat, p)))
    return MoistAir(T=float(T), D_v=float(vapour_diffusivity(T, p)),
                    W_sat=float(vapour_mass_fraction(x_sat)), rho_g=rho_g)
