"""``materials.air``: the DuMuX expressions for the gas side, and the legacy constants as one
state of them.

The port is checked against numbers that do NOT come from this project: IAPWS-IF97 table values
of the saturation pressure, and the closed form of the diffusion coefficient.  That the three
constants every record was made with are these functions at 20 C is a separate statement, pinned
to the precision the constants were written with.
"""
import numpy as np
import pytest

from saltv2.cases import dumux1d
from saltv2.coupling import driver3d
from saltv2.materials import air, brine


def test_saturation_pressure_is_iapws_if97_region_4():
    # IAPWS-IF97 (2007 revision), table 35: p_s(300 K), p_s(500 K), p_s(600 K) [MPa]
    for T, p_MPa in ((300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)):
        assert float(air.saturation_pressure(T)) == pytest.approx(p_MPa * 1e6, rel=1e-8)


def test_vapour_diffusivity_is_the_dumux_power_law():
    assert float(air.vapour_diffusivity(273.15, 1e5)) == pytest.approx(2.13e-5, rel=1e-15)
    assert float(air.vapour_diffusivity(300.0, 2e5)) == pytest.approx(
        2.13e-5 * 0.5 * (300.0 / 273.15) ** 1.8, rel=1e-15)


def test_gas_density_is_an_ideal_mixture():
    dry = float(air.gas_density(293.15, 0.0))
    assert dry == pytest.approx(1e5 * 0.02896 / (8.314472 * 293.15), rel=1e-15)
    # water vapour is lighter than air: moist air is lighter, linearly in the mole fraction
    x = 0.02
    assert float(air.gas_density(293.15, x)) == pytest.approx(
        dry * (1.0 - x * (1.0 - air.M_H2O / air.M_AIR)), rel=1e-14)


def test_the_legacy_constants_are_the_dumux_expressions_at_20_C():
    """D_V, W_SAT and RHO_G are one state: 293.15 K, 1e5 Pa, the density averaged between a
    pure-water surface and QSPNM's default environment (RH = 0.00293 / 0.01468), the rule of
    the reference case's boundary condition.  Each to the digits it was written with."""
    m = air.moist_air(air.T_LEGACY, air.RH_LEGACY, a_w_surface=1.0)
    assert round(m.D_v, 9) == air.D_V                      # 2.4189e-5: five significant digits
    assert round(m.W_sat, 5) == air.W_SAT                  # 0.01468
    assert round(m.rho_g, 5) == air.RHO_G                  # 1.18186
    assert air.LEGACY.rho_D == air.RHO_D and air.LEGACY.X_salt == air.X_SALT


def test_a_stated_temperature_reaches_every_closure_that_carries_one():
    kw = dict(dumux1d.RUN_3D, t_end=900.0, verbose=False)
    T = 295.15
    res = driver3d.run(**kw, T=T)
    expect = air.moist_air(T, kw["rh"])
    p = res["params"]
    assert (p["T"], p["D_v"], p["W_sat"], p["rho_g"]) == (T, expect.D_v, expect.W_sat, expect.rho_g)
    assert res["closure"]["T"] == T and res["closure"]["D_v"] == expect.D_v
    from saltv2.growth.prat import prat
    assert p["kr"] == pytest.approx(float(prat.kr_at_temperature(T)), rel=1e-15) and p["kr"] > prat.KR_T21
    # the vapour problem is linear in rho_g D_v W_sat, so the evaporation of the first step
    # (M_ref, what E_total is normalised by) scales with exactly that product
    legacy = driver3d.run(**kw)
    assert legacy["params"]["T"] is None and legacy["params"]["kr"] == prat.KR_T21
    assert legacy["closure"]["T"] == 294.15
    ratio = res["M_ref"] / legacy["M_ref"]
    assert ratio == pytest.approx(expect.rho_D * expect.W_sat / (air.RHO_D * air.W_SAT), rel=1e-12)
    assert 1.12 < ratio < 1.16                              # +14 %: the 20 -> 22 C mismatch
    # and the network's brine is weighed at T
    assert np.allclose(brine.salt_mass(res["pore_C"][0], res["layer"].V, T=T),
                       brine.salt_mass(kw["C0"], res["layer"].V, T=T), rtol=1e-12)


def test_temperature_is_refused_where_it_is_not_carried():
    with pytest.raises(ValueError, match="precipitation"):
        driver3d.run(**dict(dumux1d.RUN_3D, t_end=60.0, verbose=False), T=295.15,
                     pore_precipitation=True)
