"""Inputs of ``driver3d.run`` that belong to the case: the substrate grain size ``d_sub``, the
recording interval ``record_hours`` and the salt structure.  None may change a run that does not
use it.
"""
import numpy as np
import pytest

from saltv2.cases import dumux1d
from saltv2.coupling import driver3d
from saltv2.coupling.driver3d import col
from saltv2.growth.crust import HIERARCHICAL, LEGACY_MIXED, SINGLE_SCALE
from saltv2.growth.prat import prat
from saltv2.materials.brine import D_NACL

SHORT = dict(dumux1d.RUN_3D, t_end=1800.0, verbose=False)
SHORT.pop("d_sub", None)                     # these tests state the grain size themselves


@pytest.fixture(scope="module")
def default_run():
    return driver3d.run(**SHORT)


def test_the_default_grain_size_is_four_pore_radii(default_run):
    four_r = 4.0 * float(default_run["layer"].r[default_run["layer"].surface][0])
    stated = driver3d.run(**SHORT, d_sub=four_r)
    assert np.array_equal(stated["hist"], default_run["hist"])
    assert np.array_equal(stated["eng"].state, default_run["eng"].state)
    assert default_run["params"]["d_sub"] is None and stated["params"]["d_sub"] == four_r


def test_a_stated_grain_size_reaches_the_salt_closures(default_run):
    """The first hierarchy layer is one grain thick and its pores are 0.14 d_sub / 4 wide, so a
    smaller grain must give the structure finer pores from the first voxel on."""
    d_sub = 600e-6
    res = driver3d.run(**SHORT, d_sub=d_sub)
    born = np.isfinite(res["crust"].t_birth) & (res["crust"].generation == 1)
    assert born.any()
    # r_p is proportional to a_1 = 0.14 d_sub / 4, so its ratio to the default run is d_sub's
    born_default = np.isfinite(default_run["crust"].t_birth) & (default_run["crust"].generation == 1)
    ratio = res["crust"].r_birth[born].max() / default_run["crust"].r_birth[born_default].max()
    assert ratio == pytest.approx(d_sub / (4.0 * 200e-6), rel=1e-12)
    assert res["closure"]["d_sub"] == pytest.approx(d_sub, rel=1e-12)
    assert res["params"]["d_sub"] == d_sub


@pytest.mark.parametrize("bad", [0.0, -1e-4])
def test_a_grain_size_that_is_not_a_length_is_refused(bad):
    with pytest.raises(ValueError, match="d_sub"):
        driver3d.run(**SHORT, d_sub=bad)


def test_recording_less_often_changes_the_record_and_not_the_run(default_run):
    thin = driver3d.run(**SHORT, record_hours=0.1)
    assert thin["n_steps"] == default_run["n_steps"]
    assert np.array_equal(thin["eng"].state, default_run["eng"].state)
    assert thin["ledger"] == default_run["ledger"]
    assert len(thin["hist"]) < len(default_run["hist"])
    assert len(thin["pore_C"]) == len(thin["hist"]) == len(thin["M_liquid_in"])
    # every row that is written is a row of the full record, the last state included
    t_full = col(default_run["hist"], "t_h"); t_thin = col(thin["hist"], "t_h")
    rows = np.searchsorted(t_full, t_thin)
    assert np.array_equal(default_run["hist"][rows], thin["hist"])
    assert np.array_equal(default_run["pore_C"][rows], thin["pore_C"])
    assert t_thin[-1] == t_full[-1]
    assert np.all(np.diff(t_thin)[:-1] >= 0.1 - 1e-12)


@pytest.mark.parametrize("structure", [SINGLE_SCALE, HIERARCHICAL], ids=lambda s: s.name)
def test_a_stated_salt_structure_reaches_every_closure(structure):
    """Porosity and pore-size ratio travel together: the growth law's F, the birth porosity and
    the hierarchy of the voxels must all be the named structure's."""
    d_sub = 600e-6
    res = driver3d.run(**SHORT, d_sub=d_sub, salt_structure=structure)
    assert res["params"]["salt_structure"] == res["closure"]["salt_structure"] == structure.name
    assert res["closure"]["eps_s0"] == structure.eps_s
    assert res["closure"]["r_halving"] == structure.pore_ratio
    born = np.isfinite(res["crust"].t_birth)
    assert born.any()
    F = float(prat.F(prat.Da_star(eps_s=structure.eps_s, d_sub=d_sub, Ds=D_NACL)))
    assert np.allclose(res["crust"].eps_birth[born], 1.0 - (1.0 - structure.eps_s) / F, rtol=1e-12)
    # run() has already checked the salt balances (led.check raises); the record must show it
    for balance in ("salt_network", "salt_interface", "salt_precipitation"):
        assert abs(res["ledger"][f"relative_residual_{balance}"]) < 1e-10


def test_the_default_salt_structure_is_the_legacy_pairing(default_run):
    assert default_run["params"]["salt_structure"] == LEGACY_MIXED.name
    assert default_run["closure"]["eps_s0"] == 0.43 and default_run["closure"]["r_halving"] == 0.5
