"""NaCl brine properties, ported from DuMuX so that both codes use the same liquid.

Source: dumux/material/fluidsystems/brine.hh of the tree at
/home/zhixin/Dumux/dumux_theresa/dumux (the pnm2pncmin2 reference case).  Provenance tags:
[L] literature / library value, [A] assumption of this project.

    density      Batzle & Wang (1992), DOI 10.1190/1.1443207, exactly the expression DuMuX
                 uses (``Brine::density``), with the pure-water density held at the value
                 IAPWS gives for the project's conditions [L + A, see RHO_W]
    diffusion    D(NaCl in brine) = 1.54e-9 m^2/s, the constant DuMuX's Brine returns [L]
    molar masses 58.44 g/mol (NaCl) and 18.01518 g/mol (H2O, IAPWS) [L]
    activity     simplified Pitzer model for aqueous NaCl [L], from the reference case's own
                 header ``pnm2pncmin2/dumux/material/pitzernacl.hh`` (see ``water_activity``)

The mass fraction w and the mole fraction x of NaCl convert as

    x = (w / M_s) / (w / M_s + (1 - w) / M_w),      M(x) = x M_s + (1 - x) M_w

and the molar density of the brine is rho(w) / M(x).  ``saltv2.pnm.transport.salt_step`` solves
the MASS-averaged form rho(w) D A dw / L, which is the one the DuMuX reference case uses as well
(measured from its own steady profile: the mass balance returns their stated D to five digits,
the molar form would need 1.37e-9); the mole-fraction helpers are here for reading and writing
their output.

The constants below the DuMuX block (MU, GAMMA, RHO_L, C_SAT) and ``a_w`` are the project's own
liquid values; they lived with the pore network until the 2026-09-18 modularisation.
"""
from __future__ import annotations

import numpy as np

from .constants import M_H2O, M_NACL

D_NACL = 1.54e-9            # [m^2/s, L] DuMuX Brine::binaryDiffusionCoefficient
T_REF = 294.15              # [K, M] 21 C, the project's temperature
P_REF = 1.0e5               # [Pa, A] the gas pressure of the experiments
RHO_W = 998.0               # [kg/m^3, L] pure water at T_REF; DuMuX calls IAPWS H2O here,
                            # which varies by < 0.2 % over the temperatures of this project [A]

MU = 1.9e-3                 # viscosity of the saturated brine [Pa s, L] (QSPNM shared_params)
GAMMA = 0.075               # brine surface tension [N/m, L] (Wijnhorst 2025)
RHO_L = 1200.0              # brine density of the WATER book-keeping (volume <-> mass of the
                            # evaporation and supply fluxes) [kg/m^3, A]; density(C_SAT) = 1200.97
                            # at 21 C, so it stands for a saturated brine to 0.08 %.  The salt
                            # TRANSPORT uses density(w) point by point.
C_SAT = 0.264               # NaCl saturation mass fraction [-, M] (Prat 2026, Materials)


def mole_fraction(w):
    """NaCl mole fraction from its mass fraction."""
    w = np.asarray(w, float)
    n_s = w / M_NACL
    return n_s / (n_s + (1.0 - w) / M_H2O)


def mass_fraction(x):
    """NaCl mass fraction from its mole fraction."""
    x = np.asarray(x, float)
    return x * M_NACL / (x * M_NACL + (1.0 - x) * M_H2O)


def molar_mass(x):
    """Mean molar mass of the brine [kg/mol] at NaCl mole fraction x."""
    x = np.asarray(x, float)
    return x * M_NACL + (1.0 - x) * M_H2O


def density(w, T=T_REF, p=P_REF, rho_w=RHO_W):
    """Brine mass density [kg/m^3] at NaCl MASS fraction w -- Batzle & Wang (1992).

    Verbatim from DuMuX ``Brine::density``; ``rho_w`` stands in for its call to IAPWS water.
    """
    w = np.maximum(np.asarray(w, float), 0.0)
    TempC = T - 273.15
    pMPa = p / 1.0e6
    return rho_w + 1000.0 * w * (0.668 + 0.44 * w + 1.0e-6 * (
        300.0 * pMPa - 2400.0 * pMPa * w
        + TempC * (80.0 + 3.0 * TempC - 3300.0 * w - 13.0 * pMPa + 47.0 * pMPa * w)))


def molar_density(w, **kw):
    """Brine molar density [mol/m^3] at NaCl mass fraction w."""
    return density(w, **kw) / molar_mass(mole_fraction(w))


# ---------------------------------------------------------------- water activity (Pitzer)
# Simplified Pitzer model for aqueous NaCl, verbatim from the reference case's own header
# pnm2pncmin2/dumux/material/pitzernacl.hh (struct PitzerNaCl).  The parameters are hard-coded
# there for T = 298.15 K and the header carries a TODO for the missing temperature dependence;
# this project runs at 294.15 K, so using them is an ASSUMPTION [A] on top of the literature
# values [L].  They reproduce the header's own check: phi(6.14) = 1.283, a_w(6.14) = 0.753.
PITZER_A = 0.39159          # Debye-Hueckel constant [(kg/mol)^0.5, L]
PITZER_B = 2.22718          # [(kg/mol)^0.5, L]
PITZER_BETA_MX = 0.05383    # single-B Pitzer parameter [-, L]
PITZER_C_MX = 0.00134       # third-order Pitzer parameter [-, L]
PITZER_M_MAX = 6.5          # [mol/kg] their defensive clamp, see osmotic_coefficient
NU_NACL = 2.0               # ions per formula unit (1:1 electrolyte) [-, L]


def molality(w):
    """NaCl molality [mol per kg of WATER] from its mass fraction in the brine."""
    w = np.asarray(w, float)
    return w / (M_NACL * (1.0 - w))


def mass_fraction_at_molality(m):
    """Invert ``molality``: the NaCl mass fraction of a brine of molality m [mol/kg water]."""
    m = np.asarray(m, float)
    return m * M_NACL / (1.0 + m * M_NACL)


def osmotic_coefficient(m):
    """Osmotic coefficient phi [-] of aqueous NaCl at molality m -- ``PitzerNaCl::osmoticCoefficient``.

    Gibbs-Duhem counterpart of their mean ionic activity coefficient, so the two are consistent
    by construction.  The clamp to [0, 6.5] mol/kg is theirs: outside it the parameter set is not
    calibrated (NaCl saturates at 6.14), and their Newton must not be handed a NaN.  The drivers
    here never call it above saturation, so the clamp is inert on the production path.
    """
    m = np.clip(np.asarray(m, float), 0.0, PITZER_M_MAX)
    s = np.sqrt(m)
    return (1.0 - PITZER_A * s / (1.0 + PITZER_B * s)
            + m * PITZER_BETA_MX + 2.0 * m * m * PITZER_C_MX)


def water_activity(w):
    """Water activity a_w [-] of a brine at NaCl mass fraction w -- ``PitzerNaCl::waterActivity``.

    a_w = exp(-phi nu m M_w), the Gibbs-Duhem-consistent replacement for the ideal Raoult value.
    Saturation (w = 0.264, m = 6.14) gives 0.7530 against the ideal 0.900 and against the 0.75
    the project used as a linear anchor until 2026-09-17.
    """
    m = np.clip(molality(w), 0.0, PITZER_M_MAX)
    return np.exp(-osmotic_coefficient(m) * M_H2O * NU_NACL * m)


def salt_mass(w, V, **kw):
    """Dissolved NaCl [kg] held by a pore of liquid volume V at mass fraction w."""
    return density(w, **kw) * np.asarray(V, float) * np.asarray(w, float)


def mass_fraction_from_salt_mass(m, V, **kw):
    """Invert ``salt_mass``: the mass fraction whose brine holds m kg of NaCl in volume V.

    rho(w) w is strictly increasing over the physical range, so a short Newton iteration from
    the constant-density guess converges; it is solved to machine precision, not tabulated.
    """
    m = np.asarray(m, float); V = np.asarray(V, float)
    w = m / (RHO_W * V)                                   # first guess: pure-water density
    for _ in range(50):
        f = density(w, **kw) * V * w - m
        dw = 1e-8 + 1e-8 * np.abs(w)
        df = (density(w + dw, **kw) * V * (w + dw) - density(w, **kw) * V * w) / dw
        step = f / df
        w = np.maximum(w - step, 0.0)
        if np.all(np.abs(step) <= 1e-15 + 1e-14 * np.abs(w)):
            break
    else:
        raise RuntimeError("brine mass fraction did not converge")
    return w


def a_w(C):
    """Water activity of the brine at NaCl mass fraction C [L, Pitzer], clipped at pure water.

    Until 2026-09-17 this was linear, 1 at pure water to 0.75 at saturation [A].  It is now the
    simplified Pitzer model of the DuMuX reference case (``water_activity``), so the vapour
    boundary condition and the liquid it sits on come from the same property library.
    The two anchors barely move (0.7500 -> 0.7530 at C_sat, so the evaporation driver
    a_w - RH_env rises by 0.55 % and the delta_eff = 18 mm calibration with it), but the SHAPE
    in between does: a_w is up to 3.4 % higher at C = 0.6 C_sat, which is where the still-bare
    pore openings sit.  Callers clip at C_sat, so the supersaturated branch is never used.
    """
    return water_activity(np.clip(C, 0.0, None))
