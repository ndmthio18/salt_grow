"""The stacked pore column: geometry, and the steady salt profile against its analytic solution.

Regime: evaporation-wicking.  ``PoreLayer(n_layers=k)`` resolves the substrate depth as k pore
layers joined by vertical throats, with the reservoir (fixed C0, pressure datum, saturated)
under the BOTTOM layer only; k = 1 is the single surface layer every production run uses.

Reference solution.  Water leaves the top pore at rate Q and the salt stays in the network, so
at steady state the salt flux vanishes in every throat: advection up balances back-diffusion.
Per throat (advection upwind, ``salt_step``'s own discretisation, densities from saltv2.materials.brine)
that is

    Q rho_below w_below = d rho_face (w_above - w_below)
        ->  w_above / w_below = 1 + (Q / d) (rho_below / rho_face),

and in the continuum limit, with u = Q/A and the column diffusivity D = d dz / A,

    u C = D dC/dz  ->  C(z) = C_bottom exp(u z / D),   u dz / D = Q/d = Pe_cell,

i.e. the 2014 JFM section 4 profile (their Pe = J L / (rho_l A eps D_s)).  The discrete solution
is the first-order upwind form (1 + Pe)^n of exp(Pe n), so refining dz must make the two agree.
"""
import numpy as np
import pytest

from saltv2.materials import brine
from saltv2.materials.brine import C_SAT
from saltv2.pnm import PoreLayer, network_step, salt_step, solve_pressure


# water sink at the top pore [m^3/s].  The throats of this column carry d ~ 5e-16 m^3/s, so this
# is a cell Peclet number Q/d of about 0.3 at dz = 200 um -- enrichment that stays a brine.
Q_COLUMN = 1.5e-16


def _column(n_layers, dz=200e-6, pore_radius=10e-6, lam=1):
    """One pore per layer: a 1-D column.  r << dz keeps the throat length ~ dz, so that d * dz
    (the column diffusivity times area) barely moves when dz is refined."""
    return PoreLayer(lam, n_side=1, pore_radius=pore_radius, n_layers=n_layers, dz=dz,
                     Lt_res=dz)


def _steady_profile(layer, Q_evap, C_res, n_sweeps=60, dt=1e5):
    """Salt profile after the network has settled: backward Euler with a long step, no salt sink
    (the salt that advection carries up stays in the network)."""
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    m = brine.salt_mass(C_res, layer.V)
    q = np.zeros(layer.N); q[layer.surface] = Q_evap          # [m^3/s] out of the top pore
    rho = brine.density(brine.mass_fraction_from_salt_mass(m, layer.V))
    _, Q_lat, Q_res = solve_pressure(layer, q * rho, g_lat, g_res,
                                     rho, float(brine.density(C_res)))
    for _ in range(n_sweeps):
        m_new = salt_step(layer, m, C_res, Q_lat, Q_res, d_lat, d_res, dt)
        if np.max(np.abs(m_new - m)) <= 1e-14 * np.max(np.abs(m)):
            m = m_new
            break
        m = m_new
    return brine.mass_fraction_from_salt_mass(m, layer.V), Q_lat, d_lat


def test_one_layer_is_exactly_the_network_the_production_runs_use():
    a = PoreLayer(4, n_side=3)
    b = PoreLayer(4, n_side=3, n_layers=1)
    assert a.N == b.N == 9 and a.n_layers == 1
    assert np.array_equal(a.r, b.r) and a.lat == b.lat
    assert np.array_equal(a.rt_res, b.rt_res) and np.array_equal(a.surface, a.bottom)


def test_stack_geometry_and_where_the_reservoir_hangs():
    k, n_side, dz = 4, 3, 300e-6
    layer = PoreLayer(4, n_side=n_side, n_layers=k, dz=dz)
    per_layer = n_side * n_side
    assert layer.N == per_layer * k
    assert np.array_equal(layer.surface, np.arange(per_layer))
    assert np.array_equal(layer.bottom, np.arange(per_layer) + per_layer * (k - 1))
    # lateral throats in every layer (2 per pore on the periodic lattice) + one vertical per pore
    assert len(layer.lat) == 2 * per_layer * k + per_layer * (k - 1)
    # the reservoir touches the bottom layer only
    _, g_res = layer.conductances()
    assert np.all(g_res[layer.bottom] > 0.0)
    assert np.all(np.delete(g_res, layer.bottom) == 0.0)
    # the surface layer sits at z = 0 and the stack goes down by dz
    assert np.allclose(layer.z[layer.surface], 0.0)
    assert np.allclose(layer.z[layer.bottom], -(k - 1) * dz)
    # the pores of one column share a lateral position
    assert np.allclose(layer.x[layer.surface], layer.x[layer.bottom])


def test_steady_profile_is_the_zero_flux_balance_of_every_throat():
    layer = _column(6)
    C, Q_lat, d_lat = _steady_profile(layer, Q_evap=Q_COLUMN, C_res=0.05)
    assert C.max() < C_SAT                              # still a brine, no precipitation implied
    # vertical throats are the last (k - 1) entries and run surface -> down
    vertical = [(a, b, Q, d) for (a, b, _, _), Q, d in zip(layer.lat, Q_lat, d_lat)
                if abs(layer.z[a] - layer.z[b]) > 0.0]
    assert len(vertical) == layer.n_layers - 1 == len(layer.lat)   # one column: no lateral throat
    from saltv2.materials import brine
    rho = brine.density(C)
    for a, b, Q, d in vertical:                          # a above b, Q > 0 means a -> b
        ratio = C[a] / C[b]
        expected = 1.0 + (abs(Q) / d) * rho[b] / (0.5 * (rho[a] + rho[b]))
        assert ratio == pytest.approx(expected, rel=1e-8)
    assert np.all(np.diff(C) < 0.0)                      # enrichment towards the surface


def test_discrete_profile_converges_to_the_analytic_exponential():
    """C(z) = C_bottom exp(u z / D): the upwind solution (1 + Pe)^n must approach it as Pe -> 0."""
    H = 1.0e-3                                           # column height [m], fixed
    C_res, Q = 0.05, Q_COLUMN
    errors, pecl = [], []
    for k in (3, 5, 9, 17, 33):
        dz = H / (k - 1)
        layer = _column(k, dz=dz)
        C, Q_lat, d_lat = _steady_profile(layer, Q_evap=Q, C_res=C_res)
        vertical = [(Q, d) for (a, b, _, _), Q, d in zip(layer.lat, Q_lat, d_lat)
                    if abs(layer.z[a] - layer.z[b]) > 0.0]
        Pe = np.mean([abs(Q_) / d for Q_, d in vertical])     # cell Peclet number u dz / D
        depth = -layer.z                                      # 0 at the surface, H at the bottom
        analytic = C[layer.bottom[0]] * np.exp((depth[layer.bottom[0]] - depth) * (Pe / dz))
        errors.append(float(np.max(np.abs(C - analytic) / analytic)))
        pecl.append(float(Pe))
    assert pecl[0] > pecl[-1] > 0.0                      # refinement lowers the cell Peclet number
    assert errors[0] > errors[-1]                        # and the error with it
    # first-order upwind: the error falls roughly in proportion to Pe
    assert errors[-1] / errors[0] < 1.5 * pecl[-1] / pecl[0]
    assert errors[-1] < 0.02


def test_salt_is_conserved_while_the_column_settles():
    """What the reservoir hands over is what the network gained, to round-off (until 2026-09-21
    this test could only say that the column enriches: the scheme was not conservative, see
    tests/test_pnm_conservation.py)."""
    layer = _column(8)
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    C_res = 0.05; dt = 10.0
    rho_res = float(brine.density(C_res))
    m = brine.salt_mass(C_res, layer.V)
    q = np.zeros(layer.N); q[layer.surface] = Q_COLUMN
    m0 = m.sum(); net_in = 0.0
    for _ in range(20):
        rho = brine.density(brine.mass_fraction_from_salt_mass(m, layer.V))
        m_new, _, _, _, Q_res, storage = network_step(layer, m, C_res, q * rho, g_lat, g_res,
                                                      d_lat, d_res, dt)[:6]
        w_new = m_new / (layer.V * (rho + storage * dt / layer.V))
        rho_face = 0.5 * (rho + rho_res)                  # advection in minus back-diffusion out
        net_in += dt * float(((np.maximum(Q_res, 0.0) * rho_res + d_res * rho_face) * C_res
                              - (np.maximum(-Q_res, 0.0) * rho + d_res * rho_face) * w_new).sum())
        m = m_new
    C = brine.mass_fraction_from_salt_mass(m, layer.V)
    assert m.sum() > m0                                   # the column enriches
    assert np.all(C >= C_res - 1e-12)                     # nothing falls below the reservoir
    assert net_in > 0.0
    assert m.sum() - m0 == pytest.approx(net_in, rel=1e-10)


def test_the_driver_runs_on_a_column_and_keeps_its_balances():
    """The 3-D loop with a resolved substrate: the reservoir feeds the bottom layer, the openings
    are still the surface layer, and both registered balances stay tight."""
    from saltv2.coupling import driver3d as v
    from tests.conftest import CASE_3D

    res = v.run(**dict(CASE_3D, n_layers=3, t_end=4 * 3600, verbose=False))
    layer = res["layer"]
    assert layer.n_layers == 3 and layer.N == 3 * layer.n_surface
    h = res["hist"]
    assert np.abs(v.col(h, "res_water")).max() < 1e-9
    assert np.abs(v.col(h, "res_salt")).max() < 1e-9
    assert abs(res["ledger"]["relative_residual_salt_network"]) < 1e-9
    # the openings the vapour solve sees are the surface layer only
    assert res["params"]["n_layers"] == 3
    assert np.isfinite(res["t_nuc"]).sum() <= layer.n_surface


# ------------------------------------------------- the DuMuX top boundary condition (13.6)
def test_the_boundary_layer_flux_reproduces_the_dumux_reference_run():
    """``saltv2.pnm.boundary_layer_evaporation`` against the number the reference case
    itself reports, with nothing fitted.

    Measured from /home/zhixin/Dumux/dumux_theresa/dumux/pnm2pncmin2, 1-D dilute case rerun on
    2026-09-17 with SpatialParams.Temperature = 293.15: at the steady state its top pore sits at
    a NaCl mass fraction of 0.19969 and evaporation.csv reports 2.68462e-11 kg/s through the
    (2 x 200 um)^2 opening.  The residual is their averaging of rho_g and D_v between the pore
    and the free-flow reference, which are constants here.
    """
    from saltv2.cases.dumux1d import DELTA_BL, RH_ENV
    from saltv2.pnm import boundary_layer_evaporation
    w_top, area = 0.19969, 4.0 * (200e-6) ** 2
    assert float(boundary_layer_evaporation(w_top, area, DELTA_BL, RH_ENV)) == pytest.approx(2.68462e-11, rel=2e-3)


def test_the_boundary_layer_flux_is_a_difference_of_vapour_mass_fractions():
    """Its three linearities, and the sign it is NOT clamped at.

    Condensation cannot be reached by concentration alone: a_w only falls to 0.735 at the
    Pitzer clamp (6.5 mol/kg), far above the 20 % far field.  It takes a humid far field, which
    is why the function takes ``rh`` rather than reading the case's RH_ENV.
    """
    from saltv2.cases.dumux1d import DELTA_BL, RH_ENV
    from saltv2.materials.air import D_V, RHO_G, W_SAT
    from saltv2.materials.brine import a_w
    from saltv2.pnm import boundary_layer_evaporation

    def bl(C, area, delta=DELTA_BL, rh=RH_ENV):          # the kernel has no case defaults
        return boundary_layer_evaporation(C, area, delta, rh)

    area = 1.6e-7
    assert float(bl(0.2, 2 * area)) == pytest.approx(2 * float(bl(0.2, area)), rel=1e-14)
    assert float(bl(0.2, area, delta=2 * DELTA_BL)) == pytest.approx(
        0.5 * float(bl(0.2, area)), rel=1e-14)
    assert float(bl(0.0, area)) == pytest.approx(
        RHO_G * D_V * area * W_SAT * (1.0 - RH_ENV) / DELTA_BL, rel=1e-14)
    w_sat_brine = 0.264
    assert float(bl(w_sat_brine, area, rh=float(a_w(w_sat_brine)))) == pytest.approx(0.0, abs=1e-30)
    assert float(bl(w_sat_brine, area, rh=0.9)) < 0.0          # condensation, not clamped away
    assert float(a_w(0.35)) > RH_ENV                           # never negative from salt alone


def test_the_column_can_be_driven_by_its_own_water_activity():
    """``run_column(top_bc='boundary_layer')``: the draw is no longer imposed, so it falls as the
    surface pore enriches -- which is exactly what the prescribed-flux case cannot show.

    The two deltas are not interchangeable: 1.6 mm is the pore-scale layer of the DuMuX case,
    18 mm is the calibrated far field of the 2014 experiment (README section 13.6).
    """
    from saltv2.cases import pnm_column as cc
    j2014 = 1.5925e-5
    fast = cc.run_column(11, j2014, top_bc="boundary_layer")
    slow = cc.run_column(11, j2014, top_bc="boundary_layer", delta_bl=18e-3)
    fixed = cc.run_column(11, j2014)
    assert fixed["top_bc"] == "prescribed" and fast["top_bc"] == "boundary_layer"
    assert fast["M_top"] / slow["M_top"] == pytest.approx(18e-3 / 1.6e-3, rel=0.02)
    assert fast["t_c"] < fixed["t_c"] < slow["t_c"]        # 0.08 h < 0.53 h < 1.40 h
    with pytest.raises(ValueError, match="top_bc"):
        cc.run_column(11, j2014, top_bc="nonsense")


# ------------------------------------------------ the mass-conserving network (README 13.7)
def test_the_network_balances_liquid_mass_not_volume():
    """``solve_pressure`` conserves the liquid MASS the reservoir hands over, which is what a
    brine whose density runs 998 -> 1233 kg/m^3 requires (and what DuMuX does)."""
    layer = _column(6)
    g_lat, g_res = layer.conductances()
    C_res = 0.05
    # pore 0 is the surface (the one that loses liquid), the last one sits on the reservoir
    w = np.linspace(0.26, C_res, layer.N)                 # a real spread of densities
    rho = brine.density(w); rho_res = float(brine.density(C_res))
    top = int(layer.surface[0])
    q_mass = np.zeros(layer.N); q_mass[top] = 3e-12                # [kg/s] leaving the top pore
    p, Q_lat, Q_res = solve_pressure(layer, q_mass, g_lat, g_res, rho, rho_res)
    rho_up_res = np.where(Q_res > 0.0, rho_res, rho)
    assert float((rho_up_res * Q_res).sum()) == pytest.approx(float(q_mass.sum()), rel=1e-10)
    # the old VOLUME balance would have delivered q_mass / rho_top instead: 16 % more here
    old_volumetric = float(q_mass.sum()) / float(rho[top])
    assert float(Q_res.sum()) == pytest.approx(float(q_mass.sum()) / rho_res, rel=1e-10)
    assert old_volumetric / float(Q_res.sum()) == pytest.approx(rho_res / rho[top], rel=1e-10)
    assert abs(rho_res / rho[top] - 1.0) > 0.1


def test_a_uniform_density_reduces_it_to_the_old_volumetric_solve():
    """The invariant that makes the change auditable: with rho constant the density cancels, so
    the throat volume fluxes depend only on q_mass / rho and not on rho itself."""
    layer = _column(6)
    g_lat, g_res = layer.conductances()
    q_vol = np.zeros(layer.N); q_vol[layer.surface] = 2.5e-15      # [m^3/s]
    out = []
    for rho0 in (998.0, 1200.0, 1232.65):
        p, Q_lat, Q_res = solve_pressure(layer, q_vol * rho0, g_lat, g_res,
                                         np.full(layer.N, rho0), rho0)
        out.append((p, Q_lat, Q_res))
    for k in range(3):
        assert np.allclose(out[0][k], out[1][k], rtol=1e-12, atol=0.0)
        assert np.allclose(out[0][k], out[2][k], rtol=1e-12, atol=0.0)
    assert float(out[0][2].sum()) == pytest.approx(float(q_vol.sum()), rel=1e-10)
