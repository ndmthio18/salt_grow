"""Analytical growth law of Prat (2026) -- the single source of truth.

Every function is pure, stateless, SI-unit, and numpy-broadcasting, so the
figure scripts can sweep parameters without redefining any physics. Each
docstring names the paper equation it implements.

Reference
---------
M. Prat, "The Evaporation-Driven Growth Law of Porous Salt Structures",
Langmuir 2026, DOI 10.1021/acs.langmuir.6c02211.

Physical picture
----------------
A porous salt structure sits on a substrate that feeds it saturated brine by
capillarity. Evaporation at its top surface drives a Darcy flow that carries
dissolved salt upward; the salt precipitates in a thin reaction-diffusion
boundary layer at the moving surface, and the structure advances. The whole
paper reduces this to: the interface normal velocity is proportional to the
local evaporation flux, with a prefactor set by the height-independent
Damkoehler number Da*.
"""

import numpy as np

from params import (
    BETA,
    C_EQ,
    DS_T21,
    EA,
    EPS_S,
    KR_T21,
    R_GAS,
    RHO,
    RHO_CR,
)

__all__ = [
    "specific_surface",
    "Ds_eff",
    "Da",
    "Pe",
    "Da_top",
    "Da_star",
    "alpha",
    "F",
    "conc_profile",
    "boundary_layer_thickness",
    "C_top_star",
    "growth_rate",
    "growth_ratio",
    "U_over_Uevap",
    "dMs_dt",
    "kr_at_temperature",
]


# =============================================================================
# Stereology and transport properties
# =============================================================================


def specific_surface(eps_s=EPS_S, d_sub=None):
    """Pore wall area per unit volume a_v [1/m] -- eq (17).

    Eq (16) models the salt structure as a packing of cubic crystals of size
    d_c, giving a_v = 5(1-eps_s)/d_c. With the observed ratio
    d_c/d_sub ~= 0.14 this becomes a_v ~= 36(1-eps_s)/d_sub, i.e. the salt
    structure's specific surface is inherited from the substrate grain size.
    """
    return 36.0 * (1.0 - eps_s) / d_sub


def Ds_eff(eps_s=EPS_S, Ds=DS_T21):
    """Effective salt diffusivity in the pore liquid [m^2/s].

    D_s,eff = D_s / tau with the packing tortuosity tau ~= eps_s^-0.5, hence
    D_s,eff = eps_s^0.5 * D_s. Note the structure porosity appears *separately*
    in eq (1)'s diffusive term, so eps_s is deliberately not folded in here.
    """
    return np.sqrt(eps_s) * Ds


# =============================================================================
# Dimensionless groups
# =============================================================================


def Da(h, eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21):
    """Damkoehler number Da = h^2 (1-C_eq) a_v k_r / (eps_s D_s,eff) -- eq (13).

    Precipitation reaction vs diffusion, non-dimensionalised on the *whole*
    structure height, so Da grows as h^2. It exceeds 1 once the structure is
    more than a few microns thick, which is why the paper can take Da >> 1.
    """
    av = specific_surface(eps_s, d_sub)
    return h**2 * (1.0 - C_EQ) * av * kr / (eps_s * Ds_eff(eps_s, Ds))


def Pe(U_evap, h, eps_s=EPS_S, Ds=DS_T21):
    """Peclet number Pe = U_evap h / (eps_s D_s,eff) -- eq (13).

    Evaporation-induced convection vs diffusion. Typically 0-4 for the
    conditions of interest, hence Pe << Da and the convective term of eq (13)
    is dropped to reach eq (18).
    """
    return U_evap * h / (eps_s * Ds_eff(eps_s, Ds))


def Da_top(h, eps_s=EPS_S, kr=KR_T21, Ds=DS_T21, beta=BETA):
    """Top-surface Damkoehler number -- defined just below eq (22).

    Da_top = beta (1-eps_s)(1-C_eq) h k_r / (eps_s D_s,eff). It carries the
    interfacial growth kinetics of eq (20) into the top boundary condition.
    """
    return beta * (1.0 - eps_s) * (1.0 - C_EQ) * h * kr / (eps_s * Ds_eff(eps_s, Ds))


def Da_star(eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA, av=None):
    """Height-independent Damkoehler number Da* -- eq (27).

    Da* = (Da_top/sqrt(Da))^2 = beta^2 (1-eps_s)^2 (1-C_eq) k_r
                                / (eps_s D_s,eff a_v)

    Both h-dependences cancel, leaving a purely local group built from the
    structure's stereology, its diffusivity and the reaction coefficient. This
    is the single number that sets what fraction of the evaporation-supplied
    salt actually converts into front advance.

    Substituting eq (17) gives the equivalent eq (34) form
    Da* = beta^2 (1-eps_s)(1-C_eq) k_r d_sub / (36 eps_s D_s,eff),
    which makes the linear growth with substrate grain size explicit.

    Pass av to supply the salt structure's specific surface [1/m] directly
    instead of inheriting it from a substrate grain size -- an external
    efflorescence module has its own a_v,e and must not be forced to express
    it through d_sub.

    KNOWN DISCREPANCY: this reproduces the paper's *formula* but sits 8-12%
    above the three Da* values quoted in the paper's text (169.6 vs 156,
    62.5 vs 57, 102.9 vs 92). The algebra has been re-derived along two
    independent routes and checked against eq (34). Because F(Da*) is very flat
    here, no growth-rate result is affected. See the top-level README.
    """
    if av is None:
        av = specific_surface(eps_s, d_sub)
    return (
        beta**2
        * (1.0 - eps_s) ** 2
        * (1.0 - C_EQ)
        * kr
        / (eps_s * Ds_eff(eps_s, Ds) * av)
    )


def kr_at_temperature(T_K, kr_ref=KR_T21, T_ref_K=294.15, Ea=EA):
    """Arrhenius rescaling of the precipitation coefficient -- paper, Fig 9 text.

    k_r(T) = k_r(T_ref) exp[-(Ea/R)(1/T - 1/T_ref)].
    Used to take k_r from 21 C to the 45 C of the surface-heave experiment,
    where the paper quotes 10.5e-3 m/s.
    """
    return kr_ref * np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_ref_K))


# =============================================================================
# Concentration field inside the salt structure
# =============================================================================


def alpha(Da_value):
    """alpha = sqrt(Da) coth(sqrt(Da)) -- the profile's surface slope, eq (23).

    Obtained by differentiating eq (19) at z* = 1. For Da >> 1 it tends to
    sqrt(Da); for Da -> 0 it tends to 1 + Da/3. numpy's tanh saturates cleanly
    at large arguments, so only the small-Da limit needs a guard.
    """
    s = np.sqrt(np.asarray(Da_value, dtype=float))
    small = s < 1e-6
    return np.where(small, 1.0 + s**2 / 3.0, s / np.tanh(np.where(small, 1.0, s)))


def conc_profile(z_star, Da_value):
    """Normalised supersaturation profile C*(z*)/C*_top -- eq (19).

    Eq (19) reads sinh(sqrt(Da) z*)/sinh(sqrt(Da)), which OVERFLOWS for the
    paper's own parameters: at h = 1 mm and d_sub = 10 um, Da ~ 9.5e6 so
    sqrt(Da) ~ 3078 and sinh(3078) = inf, turning Figure 5b into all-NaN.

    The algebraically identical form used here,

        [exp(-s(1-z*)) - exp(-s(1+z*))] / [1 - exp(-2s)],   s = sqrt(Da)

    has every exponent <= 0 for 0 <= z* <= 1, so it never overflows. It
    degrades gracefully: for s -> 0 it returns the linear profile z*, and for
    s >> 1 it becomes exp(-s(1-z*)), the exponential boundary layer of
    thickness 1/sqrt(Da) that is the paper's central structural result.
    """
    z = np.asarray(z_star, dtype=float)
    s = np.sqrt(np.asarray(Da_value, dtype=float))
    tiny = s < 1e-8
    s_safe = np.where(tiny, 1.0, s)
    num = np.exp(-s_safe * (1.0 - z)) - np.exp(-s_safe * (1.0 + z))
    den = 1.0 - np.exp(-2.0 * s_safe)
    return np.where(tiny, z, num / den)


def boundary_layer_thickness(eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, av=None):
    """Supersaturation layer thickness delta_C [m] -- paper, after eq (19).

    delta_C = h/sqrt(Da) = sqrt(eps_s D_s,eff / ((1-C_eq) a_v k_r)).

    The two h's cancel: the layer is an intrinsic reaction-diffusion length, so
    a taller structure does not get a thicker gradient zone, it just adds more
    material sitting essentially at equilibrium. This is why the paper can skip
    discretising the structure interior -- and why anyone who *does* discretise
    it needs a mesh far finer than the structure height.

    Pass av to supply the specific surface directly (see Da_star).
    """
    if av is None:
        av = specific_surface(eps_s, d_sub)
    return np.sqrt(eps_s * Ds_eff(eps_s, Ds) / ((1.0 - C_EQ) * av * kr))


def C_top_star(U_evap, eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA):
    """Supersaturation at the moving surface, C_top - C_eq -- eq (31).

    C*_top = U_evap C_eq
             / [beta(1-eps_s)(1-C_eq)k_r + sqrt((1-C_eq) a_v k_r eps_s D_s,eff)]

    Independent of structure height, for the same reason delta_C is. For NaCl
    k_r is large, so a supersaturation of order 1e-5 already sustains the
    observed growth -- the structure is essentially at equilibrium everywhere
    except in delta_C.
    """
    av = specific_surface(eps_s, d_sub)
    De = Ds_eff(eps_s, Ds)
    denom = beta * (1.0 - eps_s) * (1.0 - C_EQ) * kr + np.sqrt(
        (1.0 - C_EQ) * av * kr * eps_s * De
    )
    return U_evap * C_EQ / denom


# =============================================================================
# Growth law
# =============================================================================


def F(Da_star_value):
    """Kinetic efficiency F = sqrt(Da*)/(1 + sqrt(Da*)) -- eq (29).

    Bounded in (0, 1): the fraction of the evaporation-supplied salt that the
    interface kinetics can actually convert into advance. F -> sqrt(Da*) when
    reaction-limited, F -> 1 when supply-limited (the mass-conservation ceiling).
    The paper stresses F cannot be taken as 1 in general -- for d_sub of a few
    hundred microns Da* = O(100) and F ~ 0.9.
    """
    s = np.sqrt(Da_star_value)
    return s / (1.0 + s)


def growth_ratio(eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA, h=None):
    """(dh/dt)/U_evap -- eq (32), the dimensionless form of the growth law.

    (dh/dt)/U_evap = rho C_eq / (rho_cr (1-eps_s)(1-C_eq)) * F(Da*)

    Pass h to use the exact eq (28), which retains tanh(sqrt(Da)); leave it None
    for the Da >> 1 form of eq (29) that the paper uses for every figure. The
    switch exists so stage 2 can quantify the error of that last approximation.

    NOTE: depends on RHO, which is back-calibrated rather than given by the
    paper (see params.RHO).
    """
    ds = Da_star(eps_s, d_sub, kr, Ds, beta)
    s = np.sqrt(ds)
    if h is None:
        kinetic = s / (1.0 + s)  # eq (29): tanh(sqrt(Da)) -> 1
    else:
        th = np.tanh(np.sqrt(Da(h, eps_s, d_sub, kr, Ds)))  # eq (28)
        kinetic = s * th / (1.0 + s * th)
    return RHO * C_EQ / (RHO_CR * (1.0 - eps_s) * (1.0 - C_EQ)) * kinetic


def growth_rate(j, eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA, h=None):
    """Salt structure growth rate dh/dt [m/s] from evaporation flux j -- eq (29).

    dh/dt = C_eq j / (rho_cr (1-eps_s)(1-C_eq)) * F(Da*)

    This is the paper's punchline: growth is *linear* in the local evaporation
    flux. Its 3D generalisation, eq (36), simply reads this along the interface
    normal, which is what makes the law usable as a moving-boundary closure.
    """
    return growth_ratio(eps_s, d_sub, kr, Ds, beta, h) * j / RHO


def U_over_Uevap(eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA, h=None):
    """Darcy velocity in the structure, normalised by U_evap -- eq (33).

    U/U_evap = 1 + (rho eps_s + rho_cr(1-eps_s)) C_eq
                   / (rho_cr(1-eps_s)(1-C_eq)) * F(Da*)

    Always > 1: the liquid drawn up must not only replace the evaporated water
    but also fill the pore space of the freshly created structure and supply
    the salt that becomes crystal.
    """
    ratio = growth_ratio(eps_s, d_sub, kr, Ds, beta, h)
    return 1.0 + (eps_s + RHO_CR / RHO * (1.0 - eps_s)) * ratio


def dMs_dt(j, area, eta, eps_s=EPS_S, d_sub=None, kr=KR_T21, Ds=DS_T21, beta=BETA):
    """Precipitated salt mass rate [kg/s] -- eq (35).

    dM_s/dt = rho_cr(1-eps_s) eta A dh/dt = C_eq j eta A F(Da*) / (1-C_eq)

    eta is the fraction of the substrate surface actually covered by
    efflorescence (patchy efflorescence does not cover it all), and A the
    substrate cross-section. Note the rho_cr and (1-eps_s) cancel, so the mass
    rate is insensitive to the structure porosity -- unlike dh/dt.
    """
    ds = Da_star(eps_s, d_sub, kr, Ds, beta)
    return C_EQ * j * eta * area * F(ds) / (1.0 - C_EQ)
