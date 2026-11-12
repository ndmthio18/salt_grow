"""saltv2.coupling.interface.PoreOpenings: the footprint geometry and bookkeeping, no physics.

The closed loop through the E1 driver is pinned bit for bit by tests/test_regression.py; these
tests check the pieces on their own.
"""
import numpy as np
import pytest

from saltv2 import Grid, PoreOpenings


def _box(n_side=2, vpp=6, pitch=100e-6, r_over_pitch=0.3857):
    a = pitch / vpp
    g = Grid((n_side * vpp, n_side * vpp, 4), a)
    s = (np.arange(n_side) + 0.5) * pitch
    X, Y = np.meshgrid(s, s, indexing="ij")
    r = np.full(n_side * n_side, r_over_pitch * pitch)
    return g, PoreOpenings(g, X.ravel(), Y.ravel(), r, n_side * pitch)


def test_footprints_are_disjoint_squares_and_the_map_agrees():
    g, op = _box()
    assert op.foot.shape == (4,) + g.shape[:2]
    assert not np.any(op.foot.sum(axis=0) > 1)                 # r < pitch / 2: no overlap
    for i in range(op.n_pores):
        assert np.all(op.pore_id[op.foot[i]] == i)
    assert np.all((op.pore_id >= 0) == op.foot.any(axis=0))
    # vpp 6, r = 0.3857 pitch: every footprint is the same 4 x 4 block of cells
    assert np.all(op.foot.sum(axis=(1, 2)) == 16)
    assert np.allclose(op.voxel_area(), 16 * g.face_area)


def test_square_footprint_converges_to_the_network_opening_area():
    """The footprint is the top face of the QSPNM cube, A_open = (2 r)^2 (PoreLayer)."""
    side = 2 * 0.3857 * 100e-6
    for vpp in (6, 12, 24, 48):
        g, op = _box(n_side=1, vpp=vpp)
        ratio = op.voxel_area()[0] / side ** 2
        # cell-centre rasterisation: each side is off by less than one cell
        assert ((side - g.a) / side) ** 2 <= ratio <= ((side + g.a) / side) ** 2
    assert abs(ratio - 1.0) < 0.1                               # vpp 48


def test_footprint_wraps_across_the_periodic_edge():
    g, _ = _box(n_side=2)
    L = 2 * 100e-6
    op = PoreOpenings(g, [0.0], [0.0], [30e-6], L)             # centred on the corner
    f = op.foot[0]
    assert f[0, 0] and f[-1, -1] and f[0, -1] and f[-1, 0]
    assert f.sum() == 4 * f[: f.shape[0] // 2, : f.shape[1] // 2].sum()


def test_per_pore_sums_each_footprint():
    g, op = _box()
    field = np.arange(np.prod(g.shape[:2]), dtype=float).reshape(g.shape[:2])
    expect = np.array([field[op.foot[i]].sum() for i in range(op.n_pores)])
    assert np.allclose(op.per_pore(field), expect)


def test_nucleation_seed_is_one_reproducible_cell_inside_its_pore_footprint():
    _, op = _box()
    for i in range(op.n_pores):
        first = op.nucleation_seed(i)
        second = op.nucleation_seed(i)
        assert first.dtype == bool
        assert first.sum() == 1
        assert np.array_equal(first, second)
        assert np.all(first <= op.foot[i])


def test_bottom_condition_values_and_blocking():
    g, op = _box()
    X_open = np.array([0.011, 0.012, 0.013, 0.014])
    blocked = np.zeros(g.shape[:2], bool)
    blocked[op.foot[1]] = True                                  # pore 1 under wet salt
    bmask, bval = op.bottom_condition(X_open, blocked)
    for i in range(op.n_pores):
        assert np.all(bval[op.foot[i]] == X_open[i])
        assert np.all(bmask[op.foot[i]]) == (i != 1)
    assert not bmask[op.pore_id < 0].any()                      # inert substrate never evaporates
    with pytest.raises(ValueError):
        op.bottom_condition(X_open[:3], blocked)
