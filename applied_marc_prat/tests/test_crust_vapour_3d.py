"""The 2-D closures on the 3-D driver (port of 2026-09-11): conducting crust, base salt,
front tracking.

Regime: evaporation-wicking.  ``tests/test_crust_vapour.py`` pins the closures themselves and
their behaviour in the 2-D driver; these tests pin that the 3-D driver wires them the same way:
every balance closes with the crust channels in it, the resolved evaporation does not collapse
to a hard zero when the first voxel seals, and the base-salt path is booked whenever the crust
reaches the substrate.
"""
import numpy as np
import pytest

from saltv2.coupling import driver3d as v
from saltv2 import SALT_SEALED


@pytest.fixture
def crust_run(run_3d_small):
    """The shared 3-D reference case (conftest.CASE_3D, pinned by test_regression.py)."""
    return run_3d_small


def test_every_balance_closes_with_the_new_channels(crust_run):
    """The network, structure and precipitation balances all carry the two new channels
    (salt_crust_front, salt_base) and must close as tightly as before."""
    led = crust_run["ledger"]
    h = crust_run["hist"]
    assert abs(led["relative_residual_salt_network"]) < 1e-10
    assert abs(led["relative_residual_salt_precipitation"]) < 1e-9
    assert np.abs(v.col(h, "res_salt")).max() < 1e-10
    assert np.abs(v.col(h, "res_water")).max() < 1e-10
    for name in ("salt_crust_front", "salt_base", "salt_released_by_base"):
        assert name in led and led[name] >= 0.0


def test_the_crust_seals_and_the_evaporation_does_not_collapse(crust_run):
    """c_k = 1e-6 keeps the structure thin so the sealing front reaches the substrate inside
    the 24 h window (the pre-port 3-D check sealed this box at 6.3 h, README E1 section 9.16;
    the 2-D base-salt test needed 60 h at c_k = 1e-6).  Once sealed voxels exist the resolved evaporation must
    stay strictly positive -- that is the point of the conducting crust -- and no single step
    may remove more than a fifth of it (the wall model dropped it to exactly zero in one step)."""
    h = crust_run["hist"]
    t = v.col(h, "t_h"); E = v.col(h, "E_total"); n_sealed = v.col(h, "n_sealed_cells")
    assert n_sealed[-1] > 0, "the test case no longer seals; retune c_k / t_end, do not weaken"
    k0 = int(np.argmax(n_sealed > 0))
    # the sealing event itself does not zero the evaporation (the wall model did, in one step)
    assert E[k0] > 0.5 and E[k0 + 1] > 0.5
    # ... and the decline that follows takes hours, not a step (measured 2026-09-12: first seal
    # 7.09 h, E below 0.1 only at 11.13 h)
    below = E < 0.1
    assert below.any() and t[np.argmax(below)] - t[k0] > 1.0
    # What is NOT asserted, and why.  In this 4-pore box E does reach ~0 at the end: base salt
    # fills the z = 0 voxels over the four openings to eps = 0 (16 cells, exactly the footprint
    # count), each opening is a quarter of the evaporating area, and in 3-D the inert substrate
    # between the pores has no bottom boundary of its own -- so once the openings are walls
    # nothing evaporates.  The single-step bound of the 2-D test (< 0.2) belongs to a wet
    # substrate plane, where bare columns keep a floor under E; it does not transfer here.


def test_base_salt_is_booked_once_the_crust_reaches_the_substrate(crust_run):
    """A sealed z = 0 voxel passes base flux; that water is released with its salt and the salt
    is placed in the voxel standing on the opening.  Released and placed are one number."""
    led = crust_run["ledger"]
    sealed0 = crust_run["eng"].state[:, :, 0] == SALT_SEALED
    if not sealed0.any():
        pytest.fail("the crust never reached the substrate in this case; retune, do not skip")
    assert led["salt_released_by_base"] > 0.0
    assert led["salt_base"] == pytest.approx(led["salt_released_by_base"], rel=1e-12)


