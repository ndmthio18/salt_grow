"""Internal consistency checks that use no numbers from the paper.

These test the mathematics of the implementation itself -- limits, degenerate
cases, numerical stability -- so they stay valid even if the paper's quoted
values turn out to be wrong. They are the counterweight to test_anchors.py.
"""

import numpy as np
import pytest

import core
from params import DAY


# =============================================================================
# Limits of the kinetic efficiency F
# =============================================================================


def test_F_is_bounded_between_zero_and_one():
    ds = np.logspace(-6, 6, 200)
    f = core.F(ds)
    assert np.all(f > 0.0) and np.all(f < 1.0)
    assert np.all(np.diff(f) > 0)  # monotonically increasing in Da*


def test_F_reaction_limited_limit():
    """Da* << 1 -> F ~ sqrt(Da*): interface kinetics throttle the growth."""
    ds = 1e-8
    assert core.F(ds) == pytest.approx(np.sqrt(ds), rel=1e-3)


def test_F_supply_limited_limit():
    """Da* >> 1 -> F -> 1: growth reaches the mass-conservation ceiling."""
    assert core.F(1e12) == pytest.approx(1.0, rel=1e-5)


def test_growth_rate_ceiling_is_never_exceeded():
    """dh/dt must stay below C_eq j / (rho_cr(1-eps)(1-C_eq)) -- below eq (28)."""
    from params import C_EQ, RHO, RHO_CR

    j = RHO * 0.01 / DAY
    for d_sub in (1e-6, 1e-5, 1e-4, 1e-3):
        eps_s = 0.43
        ceiling = C_EQ * j / (RHO_CR * (1.0 - eps_s) * (1.0 - C_EQ))
        assert core.growth_rate(j, eps_s=eps_s, d_sub=d_sub) < ceiling


# =============================================================================
# eq (28) vs eq (29): the Da >> 1 approximation
# =============================================================================


def test_exact_and_simplified_growth_agree_for_thick_structures():
    """Passing h (eq 28) must converge to the h-free eq (29) once Da >> 1.

    Da > 10 already suffices, which is why the paper drops tanh(sqrt(Da)).
    """
    simplified = core.growth_ratio(eps_s=0.43, d_sub=10e-6)
    exact = core.growth_ratio(eps_s=0.43, d_sub=10e-6, h=1e-3)
    assert exact == pytest.approx(simplified, rel=1e-6)


def test_exact_growth_is_smaller_for_very_thin_structures():
    """For a structure only a few nm thick, Da < 1 and tanh(sqrt(Da)) < 1.

    Physically: the structure is thinner than the reaction-diffusion length, so
    the bottom boundary still drains supersaturation and growth is slower than
    the Da >> 1 law predicts. This is the regime eq (29) is not valid in.
    """
    simplified = core.growth_ratio(eps_s=0.43, d_sub=10e-6)
    exact = core.growth_ratio(eps_s=0.43, d_sub=10e-6, h=1e-9)
    assert exact < simplified


# =============================================================================
# Concentration profile: correctness and numerical stability
# =============================================================================


def test_profile_satisfies_its_boundary_conditions():
    """C*(0) = 0 (eq 14) and C*(1) = C*_top (eq 15), at every Da."""
    for Da_value in (1e-3, 1.0, 1e3, 1e7):
        assert core.conc_profile(0.0, Da_value) == pytest.approx(0.0, abs=1e-12)
        assert core.conc_profile(1.0, Da_value) == pytest.approx(1.0, rel=1e-12)


def test_profile_matches_naive_sinh_where_sinh_still_works():
    """Validate the stable rewrite against the literal eq (19) at moderate Da."""
    z = np.linspace(0.0, 1.0, 51)
    for Da_value in (0.5, 2.0, 25.0, 100.0):
        s = np.sqrt(Da_value)
        naive = np.sinh(s * z) / np.sinh(s)
        assert np.allclose(core.conc_profile(z, Da_value), naive, rtol=1e-10)


def test_profile_is_finite_where_the_naive_form_overflows():
    """The Figure 5b regime: naive sinh gives inf, the stable form must not.

    Da ~ 9.5e6 at h = 1 mm, d_sub = 10 um -- exactly the paper's own example.
    """
    Da_value = core.Da(h=1e-3, eps_s=0.43, d_sub=10e-6)
    assert Da_value > 1e6
    with np.errstate(over="ignore"):  # the overflow is the point of this line
        assert not np.isfinite(np.sinh(np.sqrt(Da_value)))

    z = np.linspace(0.999, 1.0, 500)
    profile = core.conc_profile(z, Da_value)
    assert np.all(np.isfinite(profile))
    assert np.all(np.diff(profile) >= 0)  # monotone increasing toward the surface

    # The layer is 1/sqrt(Da) ~ 3.2e-4 thick in z*, so the left edge of the
    # Figure 5b window sits ~3 layer thicknesses down and must have decayed by
    # the corresponding exp(-3.08) ~ 0.046. Ten layers down it is numerically nil.
    assert profile[0] == pytest.approx(np.exp(-np.sqrt(Da_value) * 1e-3), rel=1e-6)
    assert core.conc_profile(1.0 - 10.0 / np.sqrt(Da_value), Da_value) < 1e-4


def test_profile_degenerates_to_linear_without_reaction():
    """Da -> 0 means no precipitation sink, so pure diffusion gives C* ~ z*."""
    z = np.linspace(0.0, 1.0, 21)
    assert np.allclose(core.conc_profile(z, 1e-12), z, atol=1e-6)


def test_profile_becomes_a_pure_exponential_layer_at_large_Da():
    """C*/C*_top -> exp(-sqrt(Da)(1-z*)), the singular-perturbation inner solution."""
    Da_value = 1e6
    z = np.linspace(0.99, 1.0, 100)
    expected = np.exp(-np.sqrt(Da_value) * (1.0 - z))
    assert np.allclose(core.conc_profile(z, Da_value), expected, rtol=1e-8)


# =============================================================================
# alpha
# =============================================================================


def test_alpha_large_and_small_Da_limits():
    assert core.alpha(1e8) == pytest.approx(np.sqrt(1e8), rel=1e-9)
    assert core.alpha(1e-10) == pytest.approx(1.0, rel=1e-6)
    assert core.alpha(0.03) == pytest.approx(1.0 + 0.03 / 3.0, rel=1e-3)


def test_alpha_is_the_profile_slope_at_the_surface():
    """alpha must equal d/dz* of eq (19) at z* = 1 -- the link from (19) to (23)."""
    for Da_value in (0.5, 10.0, 500.0):
        h = 1e-7
        slope = (
            core.conc_profile(1.0, Da_value) - core.conc_profile(1.0 - h, Da_value)
        ) / h
        assert slope == pytest.approx(core.alpha(Da_value), rel=1e-4)


# =============================================================================
# Boundary layer thickness
# =============================================================================


def test_boundary_layer_is_independent_of_structure_height():
    """delta_C = h/sqrt(Da) must give the same length for any h.

    Da grows as h^2, so the two h-dependences cancel exactly. This is the
    paper's key structural claim and the reason a fixed fine mesh -- not one
    scaled to the structure -- is what a numerical solver would need.
    """
    intrinsic = core.boundary_layer_thickness(eps_s=0.43, d_sub=10e-6)
    for h in (1e-5, 1e-4, 1e-3, 1e-2):
        from_Da = h / np.sqrt(core.Da(h, eps_s=0.43, d_sub=10e-6))
        assert from_Da == pytest.approx(intrinsic, rel=1e-12)


def test_boundary_layer_grows_with_substrate_grain_size():
    """Coarser substrate -> fewer crystal faces per volume -> weaker sink."""
    fine = core.boundary_layer_thickness(eps_s=0.43, d_sub=1e-6)
    coarse = core.boundary_layer_thickness(eps_s=0.43, d_sub=1e-3)
    assert coarse > fine


# =============================================================================
# Cross-relations between the derived quantities
# =============================================================================


def test_Da_star_equals_the_eq34_substituted_form():
    """Da* from (Da_top/sqrt(Da))^2 must equal the eq (34) closed form."""
    from params import BETA, C_EQ, KR_T21

    for eps_s in (0.2, 0.34, 0.43, 0.5):
        for d_sub in (1e-6, 5e-5, 1e-3):
            direct = core.Da_star(eps_s=eps_s, d_sub=d_sub)
            eq34 = (
                BETA**2
                * (1.0 - eps_s)
                * (1.0 - C_EQ)
                * KR_T21
                * d_sub
                / (36.0 * eps_s * core.Ds_eff(eps_s))
            )
            assert direct == pytest.approx(eq34, rel=1e-12)


def test_Da_star_equals_Da_top_over_sqrt_Da_squared():
    """The defining relation, checked at an arbitrary height."""
    h, eps_s, d_sub = 3e-4, 0.43, 4e-5
    lhs = core.Da_star(eps_s=eps_s, d_sub=d_sub)
    rhs = (core.Da_top(h, eps_s=eps_s) / np.sqrt(core.Da(h, eps_s=eps_s, d_sub=d_sub))) ** 2
    assert lhs == pytest.approx(rhs, rel=1e-12)


def test_darcy_velocity_matches_the_interface_mass_balance():
    """eq (33) must be reproducible straight from eq (6).

    eq (6):  rho U = [rho eps_s + rho_cr(1-eps_s)] dh/dt + j   (rho_v neglected)
    so       U/U_evap = [eps_s + (rho_cr/rho)(1-eps_s)] (dh/dt)/U_evap + 1
    """
    from params import RHO, RHO_CR

    eps_s, d_sub = 0.43, 10e-6
    ratio = core.growth_ratio(eps_s=eps_s, d_sub=d_sub)
    expected = 1.0 + (eps_s + RHO_CR / RHO * (1.0 - eps_s)) * ratio
    assert core.U_over_Uevap(eps_s=eps_s, d_sub=d_sub) == pytest.approx(expected, rel=1e-12)


def test_darcy_velocity_always_exceeds_evaporation_velocity():
    """Liquid must also fill new pore space and deliver the crystallising salt."""
    for d_sub in (1e-6, 1e-5, 1e-4, 1e-3):
        assert core.U_over_Uevap(eps_s=0.43, d_sub=d_sub) > 1.0


def test_mass_rate_is_consistent_with_the_height_rate():
    """eq (35) must equal rho_cr(1-eps_s) eta A dh/dt exactly."""
    from params import RHO, RHO_CR

    eps_s, d_sub, eta, area = 0.43, 0.58e-3, 0.76, 1e-4
    j = RHO * 7.7e-3 / DAY
    from_rate = (
        RHO_CR * (1.0 - eps_s) * eta * area * core.growth_rate(j, eps_s=eps_s, d_sub=d_sub)
    )
    assert core.dMs_dt(j, area, eta, eps_s=eps_s, d_sub=d_sub) == pytest.approx(
        from_rate, rel=1e-12
    )


def test_mass_rate_is_insensitive_to_structure_porosity():
    """rho_cr(1-eps_s) cancels in eq (35) -- unlike dh/dt, dMs/dt barely moves.

    Only the weak eps_s dependence hidden inside Da* survives.
    """
    from params import RHO

    j = RHO * 7.7e-3 / DAY
    lo = core.dMs_dt(j, 1e-4, 0.76, eps_s=0.30, d_sub=0.58e-3)
    hi = core.dMs_dt(j, 1e-4, 0.76, eps_s=0.50, d_sub=0.58e-3)
    assert abs(hi - lo) / lo < 0.05


# =============================================================================
# The identity that isolates beta
# =============================================================================


def test_Da_star_over_Da_depends_only_on_beta():
    """Da*/Da = beta^2 d_sub^2 / (1296 h^2) -- everything else cancels.

    This is what makes the paper's Figure 2a (which plots Da) a decisive test
    for the Figure 4 discrepancy (which is in Da*): k_r, D_s, eps_s, C_eq and
    the a_v prefactor all drop out, so beta is the only quantity that can shift
    Da* without shifting Da.
    """
    from params import BETA

    for eps_s in (0.2, 0.34, 0.43, 0.5):
        for d_sub in (1e-6, 2e-5, 5e-4):
            for h in (1e-5, 1e-3, 1e-2):
                lhs = core.Da_star(eps_s=eps_s, d_sub=d_sub) / core.Da(
                    h, eps_s=eps_s, d_sub=d_sub
                )
                rhs = BETA**2 * d_sub**2 / (1296.0 * h**2)
                assert lhs == pytest.approx(rhs, rel=1e-12)


def test_a_change_in_beta_leaves_Da_untouched():
    """beta enters the top boundary condition only, never the bulk Da."""
    base = core.Da(1e-3, eps_s=0.43, d_sub=10e-6)
    for beta in (1.5, 1.8, 2.0, 2.5):
        assert core.Da_star(d_sub=10e-6, beta=beta) == pytest.approx(
            core.Da_star(d_sub=10e-6, beta=2.0) * (beta / 2.0) ** 2, rel=1e-12
        )
    assert core.Da(1e-3, eps_s=0.43, d_sub=10e-6) == base
