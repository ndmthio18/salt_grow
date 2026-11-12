"""In-pore precipitation (``saltv2.pnm.precipitation``), the opt-in subflorescence closure.

Ported from the DuMuX mineralisation model, equilibrium branch: everything above C_sat becomes
solid at once, the solid takes liquid volume (V = (1 - phi_s) V_0), the throats narrow by the
same volumetric fraction, and the displaced liquid is released into the network.  The production
path does NOT use it (it caps the pores at C_sat instead), so the tests here are about the
closure itself and about the driver switch being off by default.
"""
import numpy as np
import pytest

from saltv2.materials.brine import C_SAT
from saltv2.materials.salt import RHO_CR
from saltv2.pnm import PoreLayer
from saltv2.materials import brine
from saltv2.pnm.precipitation import PoreClogged, PorePrecipitation


def _column(n_layers=4, r=2e-4, dz=6e-4):
    return PoreLayer(1, n_side=1, pore_radius=r, n_layers=n_layers, dz=dz, Lt_res=dz)


def test_nothing_happens_below_saturation():
    layer = _column()
    pp = PorePrecipitation(layer, C_SAT, RHO_CR)
    V0 = np.array(layer.V)
    m = brine.salt_mass(0.2, layer.V)
    out = pp.step(m, dt=1.0)
    assert np.array_equal(out, m)
    assert np.all(pp.m_solid == 0.0) and np.allclose(layer.V, V0)
    assert np.all(pp.last_volume_rate == 0.0)


def test_salt_is_conserved_and_the_pore_ends_up_saturated():
    layer = _column()
    pp = PorePrecipitation(layer, C_SAT, RHO_CR)
    m = brine.salt_mass(np.array([0.30, 0.28, 0.20, 0.10]), layer.V)
    total = m.sum()
    out = pp.step(m, dt=10.0)
    assert out.sum() + pp.m_solid.sum() == pytest.approx(total, rel=1e-14)
    w = brine.mass_fraction_from_salt_mass(out, layer.V)
    assert np.all(w <= C_SAT + 1e-12)                    # nothing stays supersaturated
    assert w[0] == pytest.approx(C_SAT, rel=1e-10)       # and the ones that precipitated sit on it
    assert w[1] == pytest.approx(C_SAT, rel=1e-10)
    assert w[2] == pytest.approx(0.20, rel=1e-10)        # the others are untouched


def test_the_case_temperature_is_carried_through_the_equilibrium_projection():
    layer = _column(n_layers=1)
    pp = PorePrecipitation(layer, C_SAT, RHO_CR, T=293.15)
    m = brine.salt_mass(1.1 * C_SAT, layer.V, T=293.15)
    out = pp.step(m, dt=1.0)
    w = brine.mass_fraction_from_salt_mass(out, layer.V, T=293.15)
    assert w[0] == pytest.approx(C_SAT, rel=1e-10)


def test_the_solid_takes_liquid_volume_and_releases_it_as_a_rate():
    layer = _column()
    V0 = np.array(layer.V)
    pp = PorePrecipitation(layer, C_SAT, RHO_CR)
    m = brine.salt_mass(np.array([0.32, 0.10, 0.10, 0.10]), layer.V)
    dt = 25.0
    pp.step(m, dt=dt)
    phi = pp.m_solid / RHO_CR / V0
    assert np.allclose(layer.V, (1.0 - phi) * V0)         # DuMuX's poreVolume()
    assert pp.last_volume_change[0] == pytest.approx(V0[0] - layer.V[0], rel=1e-12)
    assert pp.last_volume_rate[0] == pytest.approx(pp.last_volume_change[0] / dt, rel=1e-12)
    assert np.all(pp.last_volume_change[1:] == 0.0)


def test_throats_narrow_with_the_precipitate_and_lose_conductance():
    layer = _column()
    r0 = np.array([t[2] for t in layer.lat])
    g0, _ = layer.conductances(); d0, _ = layer.diffusive()
    pp = PorePrecipitation(layer, C_SAT, RHO_CR)
    pp.step(brine.salt_mass(np.full(layer.N, 0.40), layer.V), dt=10.0)
    r1 = np.array([t[2] for t in layer.lat])
    g1, _ = layer.conductances(); d1, _ = layer.diffusive()
    assert np.all(r1 < r0) and np.all(r1 > 0.0)
    assert np.all(g1 < g0) and np.all(d1 < d0)
    # square throats: r = 0.5 sqrt(A), A = A_0 exp(-(phi_a + phi_b)/2) -- their update multiplies
    # the CURRENT area, so the throat decays exponentially and never reaches zero
    assert np.allclose(pp._radius(pp.throat_area), r1)
    phi = pp.solid_volume_fraction
    a, b = pp.ends
    A0 = pp._area(pp.initial_r)
    assert np.allclose(pp.throat_area, A0 * np.exp(-0.5 * (phi[a] + phi[b])), rtol=1e-14)


def test_the_throat_law_is_the_one_dumux_own_output_shows():
    """Anchors from the 1-D saturated reference case, both temperatures (README section 13.6).
    A throat between a pore at phi_s and an untouched one keeps r_0 exp(-phi_s / 4)."""
    for phi_s, r_dumux in ((0.3476, 100.850e-6), (0.8882, 88.0951e-6)):
        assert 110e-6 * np.exp(-phi_s / 4.0) == pytest.approx(r_dumux, rel=1e-4)


def test_a_filling_pore_stops_the_run_instead_of_crawling():
    layer = _column(n_layers=2)
    pp = PorePrecipitation(layer, C_SAT, RHO_CR, max_solid_fraction=0.5)
    with pytest.raises(PoreClogged) as e:
        for _ in range(200):                             # keep pouring salt into pore 0
            m = brine.salt_mass(np.array([0.40, 0.10]), layer.V)
            pp.step(m, dt=1.0)
    assert e.value.fraction >= 0.5
    assert 0 in e.value.pores.tolist()


def test_the_driver_switch_is_off_by_default_and_ends_cleanly_when_on():
    """The 3-D loop: the production path caps the pores, the option precipitates in them."""
    from saltv2.coupling import driver3d as v
    from tests.conftest import CASE_3D

    off = v.run(**dict(CASE_3D, t_end=2 * 3600, verbose=False))
    assert off["precipitation"] is None
    assert off["ledger"]["salt_pnm_precipitate"] == 0.0

    on = v.run(**dict(CASE_3D, t_end=6 * 3600, pore_precipitation=True, verbose=False))
    assert on["precipitation"]["solid_mass"] > 0.0
    assert on["ledger"]["salt_pnm_precipitate"] == pytest.approx(
        on["precipitation"]["solid_mass"], rel=1e-9)
    assert on["ledger"]["salt_pnm_cap"] == 0.0            # nothing is discarded
    # Until 2026-09-19 this case ended on "pore clogged": the covered pores kept precipitating
    # inside.  A pore under salt now crystallises into the STRUCTURE (the crystal it touches), so
    # in-pore solid forms only before the opening is covered and the run reaches t_end; the
    # clogging stop itself stays covered by the PoreClogged test above.
    assert on["stop_reason"] == "t_end"
    assert on["ledger"]["salt_crystallised"] > 10 * on["precipitation"]["solid_mass"]
    assert np.abs(v.col(on["hist"], "res_salt")).max() < 1e-9
