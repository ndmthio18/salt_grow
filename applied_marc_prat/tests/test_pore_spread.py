"""Pore-to-pore grain-size spread in the 2-D driver (plan/17 P9 ported, 2026-09-07).

Regime: evaporation-wicking.  Before this, the whole 2-D box shared one d_sub, so every pore had
the same F, the same reaction band and the same hydraulic reach, and they all sealed inside one
time step -- which is what made the modelled evaporation collapse as a step instead of a decline.
The pores now come from ``saltv2.pnm.PoreLayer``, the generator the 3-D driver uses, and each
column inherits the grain size of the pore it stands on.
"""
import numpy as np

from saltv2.coupling import driver2d as fb

BASE = dict(d_sub=34e-6, n_seeds=8, cells_per_pitch=6, z_box_mm=0.6, rh=0.07, c_k=1e-5,
            eps_min=0.02, t_end_h=24.0, verbose=False, ell_over_dsub=0.25,
            snap_hours=0.0)


def test_even_cv_zero_keeps_the_dgf_spread():
    """The dgf pores are 270 +- 10 um before the lambda scaling, so cv = 0 is NOT a uniform box.

    Deriving every voxel from the box mean while the growth used the local value is exactly the
    anti-pattern plan/17 P9 forbids; the columns must differ already at cv = 0.
    """
    res = fb.run(cv=0.0, **BASE)
    d = res["d_sub_col"]
    assert d.min() < d.max()
    assert 0.005 < (d.std() / d.mean()) < 0.10          # the dgf's own few-per-cent spread
    # F and delta_C follow it
    assert res["F_col"].min() < res["F_col"].max()
    assert res["band_col"].min() < res["band_col"].max()


def test_cv_widens_the_spread_and_still_conserves_salt():
    wide = fb.run(cv=0.4, pore_seed=3, **BASE)
    narrow = fb.run(cv=0.0, **BASE)
    assert wide["d_sub_col"].std() / wide["d_sub_col"].mean() > \
        3 * narrow["d_sub_col"].std() / narrow["d_sub_col"].mean()
    for res in (wide, narrow):
        led = res["ledger"]
        assert abs(led["relative_residual_salt_structure"]) < 1e-10
        assert abs(led["relative_residual_salt_precipitation"]) < 1e-9


def test_birth_porosity_follows_the_column_not_the_box_mean():
    """eps_birth = 1 - (1 - eps_s)/F with the column's own F, so a spread in d_sub must show up
    as a spread in the porosity the material is born with."""
    res = fb.run(cv=0.4, pore_seed=3, **BASE)
    eps_b = res["eps_b_col"]
    assert eps_b.min() < eps_b.max()
    # monotone in the grain size: bigger grains -> larger Da* -> larger F -> more porous birth
    order = np.argsort(res["d_sub_col"])
    assert np.all(np.diff(eps_b[order]) >= -1e-12)


def test_the_pores_no_longer_seal_in_a_single_step():
    """With one d_sub for the whole box every column sealed together.  With the spread the sealed
    fraction has to climb through intermediate values instead of jumping straight to 1."""
    res = fb.run(cv=0.4, pore_seed=3, **BASE)
    f_sealed = res["hist"][:, 5]
    assert f_sealed.max() > 0.0
    distinct = np.unique(np.round(f_sealed[f_sealed > 0.0], 6))
    assert distinct.size >= 3, f"sealed fraction took only {distinct.size} distinct values"
