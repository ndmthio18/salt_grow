"""NaCl brine properties (``saltv2.materials.brine``), ported from the DuMuX reference case.

Split like prat2026's tests: the ANCHORS pin numbers the source states, the CONSISTENCY checks
pin limits and monotonicity with no source number in them.
"""
import numpy as np
import pytest

from saltv2.materials.brine import C_SAT, a_w
from saltv2.materials import brine


# ------------------------------------------------------------------ anchors
def test_the_pitzer_header_s_own_numerical_check():
    """pitzernacl.hh: "at saturation (m = 6.14 mol/kg) this returns phi ~ 1.283, giving
    a_w ~ 0.753 as reported in the manuscript"."""
    assert float(brine.osmotic_coefficient(6.14)) == pytest.approx(1.283, abs=5e-4)
    assert float(brine.water_activity(brine.mass_fraction_at_molality(6.14))) == \
        pytest.approx(0.753, abs=5e-4)


def test_saturation_is_their_6_14_molal_solution():
    """C_sat = 0.264 mass fraction is the same brine as the reference case's 6.14 mol/kg."""
    assert float(brine.molality(C_SAT)) == pytest.approx(6.14, rel=4e-3)


def test_the_ideal_and_the_linear_anchor_bracket_the_pitzer_value():
    """The header contrasts a_w(m_sat) = 0.753 with the ideal Raoult value 0.900; the project's
    own pre-2026-09-17 closure was a straight line to 0.75.  Pitzer sits just above the line."""
    x_water = 1.0 - brine.mole_fraction(C_SAT)
    assert float(x_water) == pytest.approx(0.900, abs=5e-4)
    assert 0.75 < float(brine.water_activity(C_SAT)) < 0.7531


def test_density_matches_batzle_wang_at_the_two_anchors():
    """Batzle & Wang through DuMuX's Brine::density: pure water falls back on RHO_W, and a
    saturated brine at 21 C is 1200.97 kg/m^3 -- which is what makes the project's constant
    RHO_L = 1200, used for the WATER book-keeping, a saturated brine to 0.08 %."""
    assert float(brine.density(0.0)) == pytest.approx(brine.RHO_W, rel=1e-12)
    assert float(brine.density(C_SAT)) == pytest.approx(1200.97, abs=5e-3)


# -------------------------------------------------------------- consistency
def test_water_activity_is_one_at_pure_water_and_decreases():
    w = np.linspace(0.0, C_SAT, 50)
    aw = brine.water_activity(w)
    assert aw[0] == pytest.approx(1.0, rel=1e-14)
    assert np.all(np.diff(aw) < 0.0)
    assert np.all((aw > 0.0) & (aw <= 1.0))


def test_the_molality_clamp_is_flat_and_never_nan():
    """Their defensive clamp at 6.5 mol/kg: above it the curve stops rather than going wild."""
    w_clamp = brine.mass_fraction_at_molality(brine.PITZER_M_MAX)
    above = brine.water_activity(np.array([w_clamp, 0.5, 0.9]))
    assert np.all(np.isfinite(above))
    assert np.allclose(above, above[0])
    assert float(brine.water_activity(np.array(-0.1))) == pytest.approx(1.0, rel=1e-14)


def test_the_driver_s_a_w_is_the_brine_s():
    """runs/E1_pre wires the vapour boundary condition to the same property library as the
    liquid it stands on -- that is the whole point of the 2026-09-17 change."""
    C = np.linspace(0.0, C_SAT, 17)
    assert np.array_equal(a_w(C), brine.water_activity(C))


def test_mass_and_mole_fractions_invert_each_other():
    w = np.linspace(0.0, 0.35, 40)
    assert np.allclose(brine.mass_fraction(brine.mole_fraction(w)), w, rtol=1e-14, atol=1e-16)


def test_salt_mass_inverts_to_machine_precision():
    w = np.linspace(1e-6, 0.35, 40)
    V = np.full_like(w, 3e-11)
    assert np.allclose(brine.mass_fraction_from_salt_mass(brine.salt_mass(w, V), V), w,
                       rtol=1e-13, atol=0.0)
