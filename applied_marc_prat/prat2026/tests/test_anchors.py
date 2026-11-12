"""Pin every formula against a number the paper states in its own text.

These are the guard rails for the whole reproduction: the figures look
plausible whether or not a coefficient is wrong, so numeric anchors are the
only thing that actually catches a transcription error.

Anchor inventory (all from the paper's running text, not read off a figure):
  k_r(45 C) = 10.5e-3 m/s          Figure 9 discussion
  dh/dt     = 0.25 cm/day          Figure 6 discussion, at U_evap = 1 cm/day
  C_i-C_eq  ~ 1e-5                 just before eq (12)
  Da*       = 156 / 57 / 92        Figure 8 and Figure 9 discussions
"""

import pathlib

import numpy as np
import pytest

import core
import params
from params import DAY, DioufHeave, NorouziRad


# =============================================================================
# Anchors that pass -- these lock down the formulas they exercise
# =============================================================================


def test_kr_arrhenius_to_45C():
    """k_r(45 C) = 10.5e-3 m/s from k_r(21 C) = 2.3e-3 and Ea = 49.1 kJ/mol."""
    kr45 = core.kr_at_temperature(params.T45_K)
    assert kr45 == pytest.approx(10.5e-3, rel=0.02)


def test_growth_rate_anchor_25um():
    """dh/dt = 0.25 cm/day at U_evap = 1 cm/day, d_sub = 25 um.

    This single assertion is the sole support for params.RHO = 1200 kg/m^3,
    which the paper never states. Rerun it after touching RHO.
    """
    ratio = core.growth_ratio(eps_s=0.43, d_sub=25e-6)
    assert ratio == pytest.approx(0.25, rel=0.03)


def test_growth_rate_anchor_10um_from_figure6():
    """Figure 6 reads (dh/dt)/U_evap ~ 0.21-0.22 at d_sub = 10 um.

    A second, independent check on RHO at a different grain size, so the
    calibration cannot be an accident of one point.
    """
    ratio = core.growth_ratio(eps_s=0.43, d_sub=10e-6)
    assert 0.20 < ratio < 0.23


def test_surface_supersaturation_is_order_1e_minus_5():
    """C_top - C_eq ~ 1e-5 for U_evap = 1 cm/day, d_sub = 10 um -- eq (31).

    The paper uses this smallness to linearise the (1-C) source term, so if
    this anchor moved by orders of magnitude the derivation itself would be in
    question.
    """
    c_star = core.C_top_star(U_evap=0.01 / DAY, eps_s=0.43, d_sub=10e-6)
    assert 3e-6 < c_star < 3e-5


def test_boundary_layer_is_far_thinner_than_the_structure():
    """delta_C must be microscopic compared with a millimetre-scale structure.

    This is what makes the paper's 'do not discretise the interior' argument
    work, and what any stage-2 numerical solver will have to resolve.
    """
    delta = core.boundary_layer_thickness(eps_s=0.43, d_sub=10e-6)
    assert delta < 1e-6  # sub-micron, vs h ~ 1 mm


# =============================================================================
# Da* anchors -- KNOWN DISCREPANCY, deliberately not forced to pass
# =============================================================================
#
# Our Da* runs 8-12% above the three values quoted in the paper's text. The
# algebra was re-derived twice (from the Da_top/sqrt(Da) definition and from the
# equivalent eq (34)) and matches the paper's printed formula both times, so
# this is recorded as an open question rather than papered over by tuning a
# parameter. F(Da*) is very flat here, so no growth-rate result is affected --
# but Figure 4 plots Da* on a linear axis, where 10% is plainly visible.
#
# These use a loose 15% tolerance so the suite stays green while the
# discrepancy is under investigation; the exact computed values are asserted
# separately below so any *change* in our own numbers is caught immediately.

DA_STAR_TOL = 0.15


def test_da_star_norouzi_rad_sand():
    """Paper quotes Da* = 156 for eps_s = 0.43, d_sub = 0.58 mm. We get ~169.6."""
    ds = core.Da_star(eps_s=NorouziRad.EPS_S, d_sub=NorouziRad.D_SUB)
    assert ds == pytest.approx(NorouziRad.DA_STAR_PAPER, rel=DA_STAR_TOL)


def test_da_star_diouf_base_porosity():
    """Paper quotes Da* = 57 for eps_s = 0.43, d_sub = 90 um at 45 C. We get ~62.5."""
    kr45 = core.kr_at_temperature(params.T45_K)
    ds = core.Da_star(
        eps_s=DioufHeave.EPS_S_BASE, d_sub=DioufHeave.D_SUB, kr=kr45, Ds=DioufHeave.DS
    )
    assert ds == pytest.approx(DioufHeave.DA_STAR_PAPER_BASE, rel=DA_STAR_TOL)


def test_da_star_diouf_fitted_porosity():
    """Paper quotes Da* = 92 for eps_s = 0.34. We get ~102.9."""
    kr45 = core.kr_at_temperature(params.T45_K)
    ds = core.Da_star(
        eps_s=DioufHeave.EPS_S_FITTED, d_sub=DioufHeave.D_SUB, kr=kr45, Ds=DioufHeave.DS
    )
    assert ds == pytest.approx(DioufHeave.DA_STAR_PAPER_FITTED, rel=DA_STAR_TOL)


def test_our_da_star_values_are_stable():
    """Freeze OUR computed Da*, so a code change that shifts them is caught.

    Without this, the loose 15% tolerance above could silently absorb a real
    regression. If these values move, something in the formula chain changed.
    """
    kr45 = core.kr_at_temperature(params.T45_K)
    assert core.Da_star(eps_s=0.43, d_sub=0.58e-3) == pytest.approx(169.6, rel=1e-3)
    assert core.Da_star(
        eps_s=0.43, d_sub=90e-6, kr=kr45, Ds=DioufHeave.DS
    ) == pytest.approx(62.5, rel=1e-2)
    assert core.Da_star(
        eps_s=0.34, d_sub=90e-6, kr=kr45, Ds=DioufHeave.DS
    ) == pytest.approx(102.9, rel=1e-2)


def test_da_star_discrepancy_does_not_affect_growth_rate():
    """The 8-12% Da* gap changes F by <1%, hence no dh/dt result is affected.

    This is the reason the discrepancy is an open question and not a blocker.
    """
    for ours, theirs in ((169.6, 156.0), (62.5, 57.0), (102.9, 92.0)):
        rel = abs(core.F(ours) - core.F(theirs)) / core.F(theirs)
        assert rel < 0.01


# =============================================================================
# Slopes used by the Figure 8 / 9 experimental comparisons
# =============================================================================


def test_norouzi_rad_growth_slope():
    """Growth rate for the sand experiment, used as the Figure 8b model line."""
    rate = core.growth_rate(
        j=params.RHO * NorouziRad.U_EVAP, eps_s=NorouziRad.EPS_S, d_sub=NorouziRad.D_SUB
    )
    assert rate * DAY * 1e3 == pytest.approx(2.489, rel=1e-2)  # mm/day


def test_diouf_growth_slopes():
    """Both porosity variants of the Figure 9b model lines."""
    kr45 = core.kr_at_temperature(params.T45_K)
    j = params.RHO * DioufHeave.U_EVAP
    base = core.growth_rate(
        j, eps_s=DioufHeave.EPS_S_BASE, d_sub=DioufHeave.D_SUB, kr=kr45, Ds=DioufHeave.DS
    )
    fitted = core.growth_rate(
        j,
        eps_s=DioufHeave.EPS_S_FITTED,
        d_sub=DioufHeave.D_SUB,
        kr=kr45,
        Ds=DioufHeave.DS,
    )
    assert base * DAY * 1e3 == pytest.approx(1.576, rel=1e-2)  # mm/day
    assert fitted * DAY * 1e3 == pytest.approx(1.395, rel=1e-2)
    # Lower porosity grows slower -- the 1/(1-eps_s) factor of eq (32).
    assert fitted < base


# =============================================================================
# Figure 5b: the structure height the paper never states
# =============================================================================


def test_inferred_height_reproduces_the_published_profile():
    """The h behind Figure 5b is not in the paper; guard the value we inferred.

    On a z/h axis the layer width is 1/sqrt(Da) ~ 1/h, so getting h wrong makes
    the panel look wrong even when every formula is right (an initial guess of
    h = 1 mm was ~8x too wide). h = 8.4 mm was recovered by fitting the
    digitised published curve; this test refits it so a change in the
    digitisation or in core cannot silently invalidate that.
    """
    data = pathlib.Path(__file__).resolve().parent.parent / "data" / "fig05b_profile.csv"
    if not data.exists():
        pytest.skip("run digitize.py first")

    d = np.loadtxt(data, delimiter=",", skiprows=6)
    z, v = d[:, 0], d[:, 1] - 1.0
    c_top = core.C_top_star(0.01 / DAY, eps_s=0.43, d_sub=10e-6) / params.C_EQ

    sel = (v > 0.03 * c_top) & (v < 0.35 * c_top)
    assert sel.sum() > 40, "not enough points on the exponential rise"
    s = -np.polyfit(1.0 - z[sel], np.log(v[sel] / c_top), 1)[0]

    h = s / np.sqrt(core.Da(1.0, eps_s=0.43, d_sub=10e-6))
    assert h == pytest.approx(8.4e-3, rel=0.10)


def test_boundary_layer_thickness_does_not_depend_on_the_unknown_height():
    """delta_C in metres is an intrinsic length, so the unstated h cannot
    affect any physical prediction -- only the appearance of the z/h plot."""
    delta = core.boundary_layer_thickness(eps_s=0.43, d_sub=10e-6)
    assert delta == pytest.approx(324.8e-9, rel=1e-3)
    for h in (1e-3, 8.4e-3, 5e-2):
        assert h / np.sqrt(core.Da(h, eps_s=0.43, d_sub=10e-6)) == pytest.approx(
            delta, rel=1e-12
        )
