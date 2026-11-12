"""Vapour conduction through the sealed crust (2026-09-07).

Regime: evaporation-wicking.  Before this change a SALT_SEALED voxel was a zero-flux wall, so
the resolved evaporation fell from its full value to exactly zero in one step and the crust leak
had to be supplied by a calibrated leak.  These tests pin the replacement: the closed
reaction band in series with the voxel's remaining bulk pore space, and a driver that keeps
evaporating and precipitating underneath the crust once it has formed.
"""
import numpy as np
import pytest

from saltv2 import BirthClosure, CrustField, Grid, GeometryEngine
from saltv2.growth.crust import T_K, effective_vapour_diffusivity
from saltv2.growth.grid import SALT_SEALED

EPS_MIN = 0.02


def _field(a=4e-6, delta_C=0.6e-6):
    g = Grid((3, 6), a)
    eng = GeometryEngine(g)
    closure = BirthClosure(d_sub=34e-6, c_k=1e-5,
                           kinetic_efficiency=0.759, reaction_band=delta_C)
    crust = CrustField(g, eng, closure, rho_cr=2160.0)
    seeds = np.zeros(g.shape, bool)
    seeds[:, 0] = True
    eng.state[seeds] = 1                              # SALT_WET, so birth() has voxels to fill
    crust.birth(seeds, 0.0)
    return g, eng, crust, seeds


def test_sealed_conductance_is_the_series_of_band_and_bulk():
    """a/D_cell = delta_C/D(eps_min) + (a - delta_C)/D(eps_bulk), with r_p/r_birth = sqrt(eps/eps_birth)."""
    a, delta_C = 4e-6, 0.6e-6
    g, eng, crust, seeds = _field(a, delta_C)
    c = crust.closure
    cond = crust.sealed_conductance(EPS_MIN)[seeds]

    eps_b = crust.eps_birth[seeds]
    r_band = crust.r_birth[seeds] * np.sqrt(EPS_MIN / eps_b)
    D_band = effective_vapour_diffusivity(EPS_MIN, r_band, c.D_v, T_K)
    D_bulk = crust.D_eff[seeds]
    f = delta_C / a
    expect = 1.0 / ((f / D_band + (1.0 - f) / D_bulk) * c.D_v)
    assert np.allclose(cond, expect, rtol=1e-12)
    # the band dominates: it is 15 % of the voxel but two orders less conductive
    assert np.all(cond < 0.1 * D_bulk / c.D_v)
    assert np.all(cond > D_band / c.D_v)


def test_band_thicker_than_the_voxel_leaves_no_bulk_term():
    """delta_C >= a: the band fills the voxel, so the conductance is the band's alone."""
    a = 0.4e-6
    g, eng, crust, seeds = _field(a, delta_C=0.6e-6)
    c = crust.closure
    cond = crust.sealed_conductance(EPS_MIN)[seeds]
    r_band = crust.r_birth[seeds] * np.sqrt(EPS_MIN / crust.eps_birth[seeds])
    expect = effective_vapour_diffusivity(EPS_MIN, r_band, c.D_v, T_K) / c.D_v
    assert np.allclose(cond, expect, rtol=1e-12)


def test_porous_argument_is_the_sealed_voxels_with_their_conductance():
    g, eng, crust, seeds = _field()
    assert crust.porous_argument(EPS_MIN) is None               # nothing sealed: no porous cell
    crust.state[seeds] = SALT_SEALED
    mask, factor = crust.porous_argument(EPS_MIN)
    assert np.array_equal(mask, seeds)
    assert np.all(factor[mask] > 0.0) and np.all(factor[mask] < 1.0)
    assert np.allclose(factor[mask], crust.sealed_conductance(EPS_MIN)[mask])


def test_sealed_conductance_needs_the_reaction_band():
    g = Grid((3, 6), 4e-6); eng = GeometryEngine(g)
    closure = BirthClosure(d_sub=34e-6, c_k=1e-5, kinetic_efficiency=0.759)
    crust = CrustField(g, eng, closure, rho_cr=2160.0)
    with pytest.raises(ValueError, match="reaction_band"):
        crust.sealed_conductance(EPS_MIN)


# --------------------------------------------------------------------- driver
def test_conducting_crust_declines_without_a_cliff_and_conserves_salt():
    """The resolved crust turns sealing into a decline, with the salt balance intact.

    The wall model this replaced (checkpoint 95c1acb) dropped E_vap from ~1.0 to exactly 0 in a
    single step as the last gas path closed; with the crust conducting, no single step may take
    away more than a fifth of the evaporation, E must still be falling rather than pinned, and
    it never reaches a hard zero.
    """
    from saltv2.coupling import driver2d as fb
    kw = dict(d_sub=34e-6, n_seeds=4, cells_per_pitch=6, z_box_mm=0.6, rh=0.07, c_k=1e-5,
              eps_min=EPS_MIN, t_end_h=22.0, verbose=False, ell_over_dsub=0.25, snap_hours=0.0)
    res = fb.run(**kw)
    E_open = res["hist"][:, 1]
    assert (res["state"] == fb.SEALED).any()                     # sealing did happen
    assert np.min(np.diff(E_open)) > -0.2                        # no cliff
    assert E_open[-1] < 0.9 * E_open[0]                          # but it does come down
    assert E_open.min() > 0.0                                    # and never reaches a hard zero
    led = res["ledger"]
    assert abs(led["relative_residual_salt_precipitation"]) < 1e-9
    assert abs(led["relative_residual_salt_structure"]) < 1e-10


# ------------------------------------------------------------- base (substrate) salt
def test_bare_substrate_salt_is_booked_and_placed_under_the_crust():
    """The brine feeding the bare substrate carries salt; none of it may vanish unbooked.

    The covered columns precipitate it into the voxel standing on them (crust base); the rest
    falls on columns the structure never covered, which a model with prescribed seeds cannot
    nucleate, and is booked as ``salt_base_uncovered`` so the precipitation balance closes.
    """
    from saltv2.coupling import driver2d as fb
    # c_k = 1e-6 keeps the structure thin (h_max ~ 0.1 mm), so the sealing front sweeps ALL of it
    # inside the window and the z = 0 row really does become conducting crust.  Since 2026-09-08
    # the front pays a whole voxel's pore space per voxel instead of a delta_C-thick layer, so at
    # the reference c_k = 1e-5 it does not reach the substrate even in 150 h -- that is the point
    # of the change, not a reason to test the base-salt path on a structure it never reaches.
    led = fb.run(d_sub=34e-6, n_seeds=4, cells_per_pitch=6, z_box_mm=0.3, rh=0.07, c_k=1e-6,
                 eps_min=EPS_MIN, t_end_h=60.0, verbose=False, ell_over_dsub=0.25,
                 snap_hours=0.0)["ledger"]
    assert led["salt_released_by_base"] > 0.0
    assert abs(led["relative_residual_salt_precipitation"]) < 1e-9
    # most of the release reaches the crust base (the rest falls on columns the structure never
    # covered, which a model with prescribed seeds cannot nucleate)
    placed = 1.0 - led["salt_base_uncovered"] / led["salt_released_by_base"]
    assert 0.5 < placed < 1.0, f"placed fraction {placed:.3f}"


def test_a_band_thicker_than_the_voxel_can_hold_is_refused_by_name():
    """delta_C / a > eps_birth / (eps_birth - eps_min): closing the band would need more salt than
    the voxel's pore space holds.  Before 2026-09-07 this surfaced as an unexplained porosity
    overfill deep inside deposit_skin; it must name the condition instead."""
    a = 0.2e-6                      # far below delta_C = 0.6 um
    g, eng, crust, seeds = _field(a, delta_C=0.6e-6)
    with pytest.raises(ValueError, match="reaction band too thick for the mesh"):
        crust.skin_capacity(EPS_MIN)


def test_the_band_model_is_valid_up_to_the_stated_ratio():
    """At delta_C = a the band exactly empties the voxel's pore space down to eps_min, which is
    the continuum limit the sealing rule stands for, and it is admissible."""
    a = 0.6e-6
    g, eng, crust, seeds = _field(a, delta_C=0.6e-6)
    cap = crust.skin_capacity(EPS_MIN)[seeds]
    cost = cap / (crust.rho_cr * a)
    assert np.allclose(cost, crust.eps_birth[seeds] - EPS_MIN, rtol=1e-12)
    assert np.all(cost < crust.eps_birth[seeds])
