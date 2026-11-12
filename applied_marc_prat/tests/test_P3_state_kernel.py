"""P3 acceptance tests of plan/17: the WET / SEALED semantics and the per-voxel
property kernel.

Each state must have its gas-side and liquid-side meaning pinned separately
(plan/17 section 3.2): a SEALED voxel is a wall to the vapour solver unless the
driver hands over its conductance through the ``porous`` argument, and a new
voxel's properties must come from an explicit flip mask -- never from a
porosity test.
"""
import numpy as np
import pytest

from saltv2 import (SALT_SEALED, SALT_WET, BirthClosure,
                    CrustField, Grid, GeometryEngine, effective_vapour_diffusivity)
from saltv2.growth.crust import HIERARCHICAL, LEGACY_MIXED, SINGLE_SCALE, SaltStructure
from saltv2.growth.vapor import VaporProperties, solve_vapor

A = 1e-4
PROPS = VaporProperties(rho_D=3.0e-5, X_wet=0.02, X_top=0.002)
D_SUB = 68e-6


def engine_with_layer(shape=(6, 6, 20), h=3):
    g = Grid(shape, A)
    eng = GeometryEngine(g)
    eng.state[..., :h] = SALT_WET
    return g, eng


def crust(g, eng, **kw):
    return CrustField(g, eng, BirthClosure(d_sub=D_SUB, D_v=2.4189e-5, **kw))


# ---------------------------------------------------------------- gas side
def test_sealed_cap_is_a_wall_by_default_and_conducts_through_the_porous_argument():
    g, eng = engine_with_layer()
    eng.state[:, :, 3:8] = SALT_SEALED       # sealed cap over the wet layer
    cf = crust(g, eng, kinetic_efficiency=0.7, reaction_band=0.6e-6)   # births every voxel
    walls = solve_vapor(g, eng.state, PROPS)                 # no porous argument: a wall
    assert walls.M_salt == 0.0 and walls.M_base == 0.0
    porous = cf.porous_argument(0.02)
    assert porous is not None and np.array_equal(porous[0], eng.state == SALT_SEALED)
    open_field = solve_vapor(g, eng.state, PROPS, porous=porous)
    assert open_field.M_salt > 0.0                           # the wet layer evaporates through it
    assert open_field.M_top == pytest.approx(open_field.M_salt, rel=1e-10)
    # the gas above the wall is still solved, but it has no source: its top flux is round-off
    assert abs(walls.M_top) < 1e-12 * open_field.M_top


def test_porous_cap_resistance_is_the_series_law_with_graded_factors():
    """A cap whose voxels have different conductances must give the harmonic (series)
    resistance sum_k a / (f_k D), not an average -- the law the sealed crust relies on."""
    g = Grid((4, 4, 24), A)
    eng = GeometryEngine(g)
    eng.state[:, :, :2] = SALT_WET
    eng.state[:, :, 2:6] = SALT_SEALED
    factors = np.array([0.5, 0.2, 0.1, 0.05])
    factor = np.ones(g.shape); factor[:, :, 2:6] = factors
    field = solve_vapor(g, eng.state, PROPS, porous=(eng.state == SALT_SEALED, factor))
    R = (g.n_z - 6) * g.a / 1.0 + g.a * np.sum(1.0 / factors)   # gas + cap, /rho_D
    expected = PROPS.rho_D * (PROPS.X_wet - PROPS.X_top) / R
    assert field.face_flux(g.z_axis, -1)[:, :, 2] == pytest.approx(expected, rel=1e-10)
    assert abs(field.imbalance) < 1e-10


def test_a_fully_sealed_cavity_carries_no_flux():
    """Wet salt walled off by SEALED voxels must not evaporate: it is not
    connected to the far field (plan/17 sections 4.1 and 4.3)."""
    g = Grid((9, 9, 14), A)
    eng = GeometryEngine(g)
    eng.state[3:6, 3:6, 2:4] = SALT_WET               # a wet block ...
    shell = np.zeros(g.shape, bool)
    shell[2:7, 2:7, 1:5] = True
    shell[3:6, 3:6, 2:4] = False
    eng.state[shell] = SALT_SEALED                     # ... completely enclosed
    cf = crust(g, eng)
    field = solve_vapor(g, eng.state, PROPS)
    assert field.M_salt == 0.0
    assert np.all(np.isnan(field.X[3:6, 3:6, 2:4]) | (field.X[3:6, 3:6, 2:4] == PROPS.X_wet))
    assert not field.active_gas[4, 4, 2]
    del cf


def test_evaporation_falls_monotonically_with_porous_cap_thickness():
    g = Grid((5, 5, 30), A)
    E = []
    for cap in range(0, 9):
        eng = GeometryEngine(g)
        eng.state[:, :, :2] = SALT_WET
        porous = None
        if cap:
            eng.state[:, :, 2:2 + cap] = SALT_SEALED
            porous = (eng.state == SALT_SEALED, 0.05)
        E.append(solve_vapor(g, eng.state, PROPS, porous=porous).M_top)
    E = np.array(E)
    assert np.all(np.diff(E) < 0.0)
    assert E[-1] < 0.5 * E[0]


# ---------------------------------------------------------------- liquid side
def test_state_machine_is_irreversible_and_records_its_times():
    g, eng = engine_with_layer(shape=(4, 4, 8), h=2)
    cf = crust(g, eng)
    # sealing takes any salt voxel whose porosity has been filled
    cf.eps_s[:, :, 1] = 0.01
    sealed = cf.seal(0.02, 300.0)
    assert sealed.sum() == 16 and np.all(eng.state[:, :, 1] == SALT_SEALED)
    assert np.all(cf.t_sealed[:, :, 1] == 300.0)
    assert cf.seal(0.02, 400.0).sum() == 0            # irreversible, not re-fired


def test_sealed_voxels_are_not_growth_candidates():
    g, eng = engine_with_layer(shape=(6, 6, 10), h=1)
    cf = crust(g, eng)
    eng.state[:, :, 0] = SALT_SEALED
    assert not eng.salt_mask().any()
    rate = eng.volume_rate(1e-8)
    assert rate.max() == 0.0, "only SALT_WET faces may drive growth"
    assert cf.volumes()["wet"] == 0.0


# ---------------------------------------------------------------- properties
def test_birth_properties_come_from_the_flip_mask_only():
    g = Grid((6, 6, 12), A)
    eng = GeometryEngine(g)
    eng.state[:, :, 0] = SALT_WET
    cf = crust(g, eng)
    assert np.all(np.isfinite(cf.t_birth[:, :, 0]))           # seeds got theirs
    assert np.all(np.isnan(cf.t_birth[:, :, 1:]))
    rate = np.zeros(g.shape); rate[2, 2, 1] = 1.0 * g.cell_volume
    rep = eng.advance(rate, 1.0)
    assert rep.flipped.sum() == 1 and rep.flipped[2, 2, 1]
    cf.birth(rep.flipped, 500.0)
    assert cf.t_birth[2, 2, 1] == 500.0
    assert cf.eps_s[2, 2, 1] == pytest.approx(cf.closure.eps_s0)
    with pytest.raises(ValueError, match="twice"):
        cf.birth(rep.flipped, 600.0)


def test_hierarchy_gives_coarse_permeability_inside_and_fine_threshold_outside():
    """The parameter rule of CLAUDE.md: k_s and P_cth are independent inputs.
    Deeper (older, lower) voxels keep the coarse k_s; the newest layer carries
    the fine P_cth."""
    c = BirthClosure(d_sub=500e-6)
    z = np.array([0.0, 0.4e-3, 0.9e-3, 2.0e-3, 5.0e-3])
    n = c.generation(z)
    assert list(n) == [1, 1, 2, 3, 4]
    k = c.permeability(n); Pc = c.capillary_threshold(n)
    assert np.all(np.diff(k) <= 0.0) and np.all(np.diff(Pc) >= 0.0)
    # k_s(1) for a 240 um sand should land near the Piotrowski 2020 measurement
    # 4.8e-12 m^2 to within the O(1) prefactor c_k (knowledge base anchor)
    k1 = float(BirthClosure(d_sub=240e-6, c_k=1.0).permeability(np.array(1.0)))
    assert 1e-12 < k1 < 1e-10
    # each layer halves the pore, so P_cth doubles: exactly the 2x jumps that
    # make need - P_cth only piecewise monotonic (first_failure_height)
    assert Pc[3] / Pc[2] == pytest.approx(2.0)


@pytest.mark.parametrize("structure", [LEGACY_MIXED, HIERARCHICAL, SINGLE_SCALE],
                         ids=lambda s: s.name)
def test_closure_and_supply_check_describe_the_same_stack(structure):
    """crust.BirthClosure and growth.supply.supply_need must describe the same stack for EVERY
    salt structure, or the properties a voxel is born with contradict the supply test that
    decides whether it can be fed."""
    from saltv2.growth.supply import supply_need
    from saltv2.materials.brine import GAMMA, MU, RHO_L
    from saltv2.materials.constants import G

    c = BirthClosure(d_sub=D_SUB, c_k=3e-3, gamma=GAMMA, structure=structure)
    U = 1e-6
    for h in (1e-5, 1e-4, 5e-4, 2e-3):
        n = int(c.generation(np.array(h)))
        need, Pc = supply_need(h, U, D_SUB, 0.0, 3e-3, structure)
        assert float(c.capillary_threshold(np.array(float(n)))) == pytest.approx(Pc, rel=1e-12)
        # the viscous drop is the closure's own k_s summed over the layers below h
        layers = np.arange(1, n + 1, dtype=float)
        top = D_SUB * (2.0 ** layers - 1.0)                       # h_1 = d_sub, doubling
        thick = np.minimum(top, h) - np.concatenate([[0.0], top[:-1]])
        dp = float(np.sum(MU * U * thick / c.permeability(layers)))
        assert need == pytest.approx(RHO_L * G * h + dp, rel=1e-12)


def test_a_single_scale_structure_has_no_height_in_its_properties():
    c = BirthClosure(d_sub=D_SUB, c_k=1.0, structure=SINGLE_SCALE)
    n = np.arange(1.0, 9.0)
    for prop in (c.pore_radius(n), c.permeability(n), c.capillary_threshold(n)):
        assert np.all(prop == prop[0])
    assert c.eps_s0 == 0.43 and BirthClosure(structure=HIERARCHICAL).eps_s0 == 0.20
    # the first layer is the same slit stack in all three structures
    for other in (HIERARCHICAL, LEGACY_MIXED):
        o = BirthClosure(d_sub=D_SUB, c_k=1.0, structure=other)
        assert float(o.permeability(np.array(1.0))) == float(c.permeability(np.array(1.0)))


def test_a_salt_structure_cannot_be_given_by_half():
    """Pore-size ratio and porosity are one record: the closure has no porosity of its own."""
    with pytest.raises(TypeError):
        BirthClosure(d_sub=D_SUB, eps_s0=0.20)
    with pytest.raises(ValueError):
        SaltStructure("no_such_salt", 0.0, 0.43)
    with pytest.raises(ValueError):
        SaltStructure("no_such_salt", 0.5, 1.0)


def test_hierarchy_matches_the_supply_check_of_the_pre_experiment():
    from saltv2.materials.brine import GAMMA, MU

    c = BirthClosure(d_sub=D_SUB, c_k=3e-3, gamma=GAMMA)
    # viscous slope mu U / k_n of the two implementations must agree as well
    U = 1e-6
    n = 3.0
    k_n = float(c.permeability(np.array(n)))
    a3 = 0.14 * D_SUB / 4 * 0.5 ** 2
    assert k_n == pytest.approx(3e-3 * 2 ** 3 * a3 ** 3 / (12 * 0.14 * D_SUB), rel=1e-12)
    assert MU * U / k_n > 0.0


def test_effective_diffusivity_is_knudsen_limited_only_for_nanopores():
    D_v = 2.4189e-5
    coarse = effective_vapour_diffusivity(0.43, 1e-5, D_v, 294.15)
    nano = effective_vapour_diffusivity(0.25, 3.5e-9, D_v, 294.15)
    from saltv2.growth.crust import knudsen_diffusivity
    assert coarse == pytest.approx(0.43 ** 1.5 * D_v, rel=0.05)   # bulk limit
    # Knudsen limit: D_K(3.5 nm) = 6.9e-7 m^2/s << D_v, so D_eff -> eps^1.5 D_K
    assert nano == pytest.approx(0.25 ** 1.5 * knudsen_diffusivity(3.5e-9, 294.15), rel=0.05)
    assert nano < 1e-2 * D_v and nano > 0.0


def test_per_voxel_d_sub_is_used_instead_of_one_global_mean():
    """plan/17 P9 anti-pattern: with a pore-size distribution, every voxel must be born from the
    grain size under its own column, not from the box average."""
    g = Grid((4, 4, 10), A)
    eng = GeometryEngine(g)
    eng.state[:, :, 0] = SALT_WET
    d_map = np.full((4, 4), 100e-6)
    d_map[:2] = 25e-6                                   # half the box sits on a finer substrate
    c = BirthClosure(d_sub=float(d_map.mean()))
    cf = CrustField(g, eng, c, d_sub_field=d_map)
    # fine side: smaller pores -> lower k_s, higher P_cth
    assert cf.k_s[0, 0, 0] < cf.k_s[3, 3, 0]
    assert cf.P_cth[0, 0, 0] > cf.P_cth[3, 3, 0]
    assert cf.k_s[0, 0, 0] == pytest.approx(float(c.permeability(np.array(1.0), 25e-6)), rel=1e-12)
    assert cf.k_s[3, 3, 0] == pytest.approx(float(c.permeability(np.array(1.0), 100e-6)), rel=1e-12)
    # without the field every voxel would get the mean, which is neither
    flat = CrustField(g, GeometryEngine(g, eng.state), c)
    assert flat.k_s[0, 0, 0] == flat.k_s[3, 3, 0]
    assert flat.k_s[0, 0, 0] / cf.k_s[0, 0, 0] > 10.0   # k_s ~ d_sub^2: 25 um vs 62 um
    with pytest.raises(ValueError, match="lateral or the full grid shape"):
        CrustField(g, GeometryEngine(g, eng.state), c, d_sub_field=np.ones(3))
