"""E2 2-D feedback model (saltv2/coupling/driver2d.py): consistency checks on the sealing and
leakage accounting that the Codex review (2026-09-04) asked for."""
import numpy as np

from saltv2.coupling import driver2d as fb
from saltv2.materials.salt import EPS_S
from tests.conftest import SMALL_2D as SMALL


def test_unlimited_supply_freezes_no_face():
    """Control for the supply-failure chain: with unlimited supply no face is ever FROZEN, so
    nothing feeds a reaction band and nothing seals."""
    res = fb.run(**SMALL, c_k=1.0)
    h = res["hist"]
    assert not (res["state"] == fb.SEALED).any()
    assert np.array_equal(h[:, 1], h[:, 2])              # E_tot is E_vap (no parameterised leak)
    assert h[:, 7].max() == 0                            # no frozen face


def test_unlimited_supply_does_not_seal_because_the_front_carries_the_share_away():
    """Task B.  The (1 - F) share is no longer dumped into the voxel behind an advancing face;
    it is carried by the porosity the new material is born with.  A front that keeps advancing
    therefore never closes its own reaction band, at any grain size -- the ratio of the band's
    filling time to the time the front needs to traverse it is
    (eps_s - eps_min)/(1 - eps_s) * sqrt(Da*) > 1 throughout the range of interest.
    What the share does instead is make the structure denser."""
    res = fb.run(**SMALL, c_k=1.0)
    assert not (res["state"] == fb.SEALED).any()
    assert res["m_skin"].max() == 0.0                     # no face ever stalled
    assert res["eps_birth"] < EPS_S                       # but the material is denser than 0.43
    assert res["hist"][-1, 1] > 0.5                       # and the surface is not shut down


def test_sealing_is_irreversible():
    """Sealed volume never decreases.  Until 2026-09-06 this test also asserted that a frozen
    face must come first; that is no longer the only route, because an advancing face leaves its
    (1 - F) share in the pores as well, so sealing can precede any freezing."""
    res = fb.run(**SMALL, c_k=1e-7)
    h = res["hist"]
    sealed_vol = h[:, 6]
    assert np.all(np.diff(sealed_vol) >= -1e-30)         # sealed volume never decreases
    assert (sealed_vol > 0).any()


def test_dt_factor_does_not_produce_zeno_steps_at_sealing_or_end():
    """Regression (2026-09-05): scaling the EVENT limiters (dt_seal, t_end - t) by r < 1 approaches
    each event geometrically, tens of ever smaller steps per sealing voxel.  Only the accuracy
    limiters may be scaled; every event must be reached in one step."""
    res = fb.run(**SMALL, c_k=1e-7, dt_factor=0.25)
    h = res["hist"]; dt = np.diff(h[:, 0]) * 3600
    n_sealed = int(round(h[-1, 6] / res["a"] ** 2))
    assert n_sealed > 0
    tiny = int((dt < 1e-3 * np.median(dt)).sum())
    assert tiny <= n_sealed + 1, (tiny, n_sealed)
    assert abs(h[-1, 0] - 12.0) < 1e-9
