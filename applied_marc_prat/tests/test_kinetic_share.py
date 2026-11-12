"""The (1 - F) share of the Prat growth law (ideas/计算结果系统评估_2026-09-06.md section 3).

A face evaporating j hands over ``beta * j * A`` of salt, ``beta = C_sat/(1 - C_sat)``.  The
growth law converts the fraction ``F(Da*)`` of it into interface advance, and only where the face
may advance.  The rest has to precipitate into the pores of the voxel behind the face, or the
model silently destroys salt -- which is what happened until 2026-09-06 (the driver used F in all
three places: growing faces, frozen faces and internal wet/dry fronts).

These tests read the shares off the drivers' own ledgers, so they fail if any of the three places
regresses.  They say nothing about whether putting the leftover into the pores is the RIGHT
physical closure; that is an assumption [A] and the paper's own open question.
"""
import numpy as np
import pytest

from saltv2.coupling import driver2d as fb
from saltv2.coupling import driver3d as v
from saltv2.growth.prat import prat
from saltv2.materials.brine import D_NACL
from saltv2.materials.salt import EPS_S
from tests.conftest import SMALL_2D, SMALL_3D

E1_SMALL = dict(SMALL_3D, ell_over_dsub=0.25)
E2_SMALL = dict(SMALL_2D, ell_over_dsub=0.25)


def _placed(led):
    return (led["salt_grown"] + led["salt_front_external"] + led["salt_crust_front"]
            + led["salt_base"])


def test_3d_places_every_kilogram_evaporation_hands_over():
    """The balance the driver registers: what evaporation released must be what was placed.

    The two sides are measured independently -- the left from the growth law, the right from the
    vapour field -- so a share that does not sum to one shows up here.
    """
    res = v.run(**E1_SMALL, t_end=6 * 3600.0)
    led = res["ledger"]
    assert led["salt_released_by_evaporation"] > 0.0
    assert led["relative_residual_salt_precipitation"] == pytest.approx(0.0, abs=1e-9)
    assert _placed(led) == pytest.approx(
        led["salt_released_by_evaporation"] + led["salt_released_by_base"], rel=1e-9)


def test_3d_growing_faces_place_all_of_their_salt_through_the_birth_porosity():
    """Task B: an advancing face does not hand its (1 - F) share to the voxel behind it any
    more.  The share is carried by the porosity the new material is BORN with,

        eps_birth = 1 - (1 - eps_s0) / F,

    so rho_cr (1 - eps_birth) V_n is exactly what evaporation delivered and the pore channels
    stay empty.  Neither expression contains a grid length, which is the point of task B."""
    res = v.run(**E1_SMALL, t_end=6 * 3600.0, c_k=1.0)     # c_k = 1: supply never fails
    led = res["ledger"]
    assert led["salt_front_external"] == 0.0
    assert led["salt_grown"] == pytest.approx(led["salt_released_by_evaporation"], rel=1e-9)
    # eps_birth is per voxel: F comes from the pore that feeds the column (plan/17 P9), and even
    # at cv = 0 the dgf pores are 270 +- 10 um.  Every voxel must therefore sit inside the
    # bracket the pore spread allows, and the mean must be the mean-F value.
    F = np.asarray(prat.F(prat.Da_star(eps_s=EPS_S, d_sub=4 * res["layer"].r, Ds=D_NACL)), float)
    lo = 1.0 - (1.0 - EPS_S) / F.min()
    hi = 1.0 - (1.0 - EPS_S) / F.max()
    born = res["crust"].eps_birth[res["eng"].solid_mask_full()]
    assert lo - 1e-12 <= born.min() and born.max() <= hi + 1e-12
    # The population mean must sit NEAR the mean-F value, not at an edge of that bracket: a mean
    # pushed to an edge would mean the growth had concentrated on the widest or narrowest pore.
    # Measured as a fraction of the bracket width (the only scale this statement has): 0.068 with
    # the linear a_w that this project used until 2026-09-17, 0.096 with the Pitzer one -- the
    # fresher, more humid bare openings of the Pitzer curve skew the growth slightly more.  The
    # bound below is twice that; it is a claim about the spread, not a numerical tolerance.
    assert abs(born.mean() - (1.0 - (1.0 - EPS_S) / F.mean())) < 0.2 * (hi - lo)
    assert hi < EPS_S                                               # denser than self-replication


def test_birth_porosity_and_reaction_band_contain_no_grid_length():
    """The two closures task B introduces are functions of F and delta_C only."""
    from saltv2 import BirthClosure, CrustField, GeometryEngine, Grid, SALT_WET
    d_sub = 34e-6
    F = float(prat.F(prat.Da_star(eps_s=0.43, d_sub=d_sub)))
    dC = float(prat.boundary_layer_thickness(eps_s=0.43, d_sub=d_sub))
    closure = BirthClosure(d_sub=d_sub, kinetic_efficiency=F, reaction_band=dC)
    caps, eps = [], []
    for a in (5e-6, 2.5e-6, 1.25e-6):
        g = Grid((4, 4, 6), a)
        eng = GeometryEngine(g); eng.state[..., :2] = SALT_WET
        cf = CrustField(g, eng, closure)
        caps.append(cf.skin_capacity(0.02)[eng.solid_mask_full()][0])
        eps.append(cf.eps_s[0, 0, 0])
    assert caps[0] == pytest.approx(caps[1]) == pytest.approx(caps[2])
    assert eps[0] == pytest.approx(eps[1]) == pytest.approx(eps[2])
    assert eps[0] == pytest.approx(1.0 - (1.0 - 0.43) / F)
    # and F <= 1 - eps_s0 is the crust limit, refused rather than clipped
    with pytest.raises(ValueError, match="crust limit"):
        BirthClosure(kinetic_efficiency=0.5).birth_porosity(np.array(1.0))


def test_3d_water_ledger_separates_the_bare_substrate_share():
    """Salt released by evaporation is counted on WET SALT interfaces only: the bare openings
    evaporate too, but their salt stays in the pore layer and drives nucleation instead.

    The water side of that statement is its own channel, ``water_evaporated_on_salt``, measured
    on the same faces and with the same ``max(j, 0)`` as the salt, so the two are one identity.
    It is NOT the difference of the two signed totals: while a pore opening is still bare its
    brine is fresher than the salt's (higher ``a_w``), and its plume condenses on the salt faces
    around it.  The growth closure has no dissolution branch and drops those faces, so the
    signed difference is the SMALLER number -- by 1.0e-7 relative in this case (30 steps of 273,
    ending exactly when the last opening closes; |j_neg| / j_max = 1.5e-5).  That share is
    bounded here as a physical claim, not used as a tolerance.
    """
    res = v.run(**E1_SMALL, t_end=6 * 3600.0)
    led = res["ledger"]
    assert 0.0 < led["water_evaporated_base"] < led["water_evaporated"]
    beta = v.C_SAT / (1 - v.C_SAT)
    on_salt = led["water_evaporated_on_salt"]
    assert on_salt > 0.0
    # Until 2026-09-19 this was the identity released = beta(C_sat) * on_salt: every face handed
    # over the saturated value whatever the network had delivered.  Now beta(C_sat) -- Prat (2026)
    # eq (29) -- is the upper bound, reached when the pores' crystallised salt covers it; in this
    # run the four seeds hold more salt than their pores had given, and the first hours of growth
    # pay that back (5 % of the total).  The identity that replaces it is the interface balance.
    assert 0.9 * beta * on_salt < led["salt_released_by_evaporation"] <= beta * on_salt * (1 + 1e-12)
    assert led["relative_residual_salt_interface"] == pytest.approx(0.0, abs=1e-10)
    signed = led["water_evaporated"] - led["water_evaporated_base"]
    assert 0.0 < signed <= on_salt                       # condensation only ever removes from it
    assert on_salt / signed - 1.0 < 1e-6                 # and it stays a rounding-level share


def test_2d_places_every_kilogram_evaporation_hands_over():
    res = fb.run(**E2_SMALL, c_k=3e-3)
    led = res["ledger"]
    assert led["salt_released_by_evaporation"] > 0.0
    assert led["relative_residual_salt_precipitation"] == pytest.approx(0.0, abs=1e-9)


def test_2d_a_frozen_face_leaves_all_of_its_salt_in_a_band_of_physical_thickness():
    """A frozen face cannot advance, so its whole delivery goes into the reaction band behind
    it.  The band closes at (eps_birth - eps_min) rho_cr delta_C, a mass per unit area: that is
    what seals the voxel, and its closing time carries no grid length (task B).  Once sealed,
    the voxel is the crust front: the faces evaporating into it keep filling its BULK until it
    is down to eps_min before the front moves on (front tracking, 2026-09-08), so a sealed voxel
    ends at eps_min without ever going below it."""
    res = fb.run(**E2_SMALL, c_k=1e-7)                                # supply fails early
    led = res["ledger"]
    assert led["salt_front"] > 0.0                                    # something froze
    assert led["relative_residual_salt_precipitation"] == pytest.approx(0.0, abs=1e-9)
    sealed = res["state"] == fb.SEALED
    assert sealed.any()
    cap = res["skin_capacity"]
    assert (res["m_skin"][sealed] >= cap[sealed] * (1 - 1e-9)).all()   # every band is closed
    eps = res["eps_s"][sealed]
    assert eps.min() >= -1e-14                                        # never overfilled
    assert (eps <= 0.02 + 1e-12).mean() > 0.5                         # the front finished most of them


def test_2d_sealing_time_and_extent_do_not_depend_on_the_grid():
    """The acceptance test of task B.  With the physics held fixed, refining the mesh must not
    change when the structure seals nor how much of it seals.  Before task B the same scan gave
    'never / 0.50 h / 1.13 h' for vpp 4 / 8 / 16, because the band was one voxel thick.

    Supply fails at once here (c_k = 1e-7), so the whole low structure seals and the crust
    front then fills it to the substrate: the sealed region is an AREA (count x a^2) and its
    height is that of the structure, both of which must be mesh-independent.  (Under the wall
    model this replaced, exactly one row of cells sealed and the sealed LENGTH count x a was the
    invariant; see checkpoint 95c1acb.)  Measured 2026-09-14: t_seal 1.07 / 0.96 / 0.99 h,
    sealed area 1457 / 1270 / 1228 um^2, height 16.5 / 16.5 / 15.2 um for vpp 4 / 8 / 16."""
    t_seal, area, height, a_all = [], [], [], []
    for vpp in (4, 8, 16):
        r = fb.run(d_sub=34e-6, n_seeds=4, cells_per_pitch=vpp, z_box_mm=1.0, t_end_h=60.0,
                   c_k=1e-7, eps_min=0.02, verbose=False, ell_over_dsub=0.25,
                   snap_hours=None)
        h = r["hist"]; se = h[:, 6] > 0
        assert se.any()
        t_seal.append(h[se, 0][0])
        a = (0.7e-3 * 34e-6 / 1080e-6) / vpp
        sealed = r["state"] == fb.SEALED
        area.append(int(sealed.sum()) * a * a)                         # sealed area [m^2]
        height.append((int(np.where(sealed.any(axis=0))[0].max()) + 1) * a)
        a_all.append(a)
    assert max(t_seal) / min(t_seal) < 1.25                           # was unbounded before
    assert max(area) / min(area) < 1.25                               # same physical extent
    assert max(height) - min(height) <= max(a_all)                    # to within the coarsest voxel
