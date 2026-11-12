"""V0 acceptance tests (ideas/plan/11 section 4): geometry engine and
morphology observables, no physics."""
import numpy as np
import pytest

from saltv2 import GAS, SALT_WET, Grid, GeometryEngine, morphology as mo


A = 1e-4   # 100 um voxels


def step(eng, Vn, dt=None, f_max=0.5):
    """One occupancy step, the way the drivers compose it: rate, dt limit, advance."""
    rate = eng.volume_rate(Vn)
    if dt is None:
        dt = eng.dt_occupancy(rate, f_max)
    return eng.advance(rate, dt)


def flat_engine(nx=8, nz=40, h0=3):
    g = Grid((nx, nz), A)
    eng = GeometryEngine(g)
    eng.state[:, :h0] = SALT_WET
    return g, eng


# ---------------------------------------------------------------- flat layer
def test_flat_layer_height_equals_Vn_t():
    g, eng = flat_engine()
    Vn = 2.0e-8                      # m/s, uniform on every exposed face
    T = 5 * A / Vn                   # long enough for 5 cells
    t = 0.0
    while t < T - 1e-12:
        dt = min(0.3 * A / Vn, T - t)
        step(eng, Vn, dt=dt)
        t += dt
    h = g.column_max_height(eng.salt_mask()) + eng.f.max(axis=1) * A
    # every column identical, and equals initial height + Vn*T
    assert np.allclose(h, h[0])
    assert abs(h[0] - (3 * A + Vn * T)) < 1e-12 * A + 1e-9


def test_volume_conservation_with_overflow():
    g, eng = flat_engine(nx=6, nz=30, h0=2)
    V0 = eng.volume_placed
    rng = np.random.default_rng(0)
    for _ in range(40):
        Vn = 1e-8 * (0.5 + rng.random())     # random but positive
        step(eng, Vn, f_max=0.7)             # large steps -> overflow paths exercised
    assert abs(eng.conservation_residual(V0)) < 1e-12


# ---------------------------------------------------------------- symmetry
def test_single_seed_uniform_growth_is_mirror_symmetric():
    g = Grid((21, 21), A)
    eng = GeometryEngine(g)
    eng.seed(10, height_cells=1, width_cells=1)
    for _ in range(60):
        step(eng, 1e-8, f_max=0.5)
    s = eng.salt_mask()
    assert np.array_equal(s, s[::-1, :])            # mirror in x about the seed
    assert np.allclose(eng.f, eng.f[::-1, :])


def test_two_seeds_merge_and_cluster_count_drops():
    g = Grid((30, 20), A)
    eng = GeometryEngine(g)
    eng.seed(8, 1, 1)
    eng.seed(20, 1, 1)
    counts = []
    for _ in range(80):
        step(eng, 1e-8, f_max=0.5)
        counts.append(mo.n_clusters(eng.salt_mask(), g))
    assert counts[0] == 2
    assert counts[-1] == 1
    assert all(c2 <= c1 for c1, c2 in zip(counts, counts[1:]))


# ---------------------------------------------------------------- morphology
def test_periodic_cluster_counted_once():
    g = Grid((10, 5), A)
    s = np.zeros(g.shape, bool)
    s[0, 0] = s[9, 0] = True          # touches across the periodic edge
    assert mo.n_clusters(s, g) == 1
    assert not mo.spanning(s, g)


def test_spanning_requires_full_line_in_2d():
    g = Grid((10, 5), A)
    s = np.zeros(g.shape, bool)
    s[:, 0] = True
    assert mo.spanning(s, g)
    s[4, 0] = False                   # one bare column breaks the line
    assert not mo.spanning(s, g)


def test_height_stats():
    g = Grid((12, 6), A)
    s = np.zeros(g.shape, bool)
    s[2, :3] = True                   # height 3a
    s[8, :1] = True                   # height 1a
    hs = mo.height_stats(s, g)
    assert hs["h_max"] == pytest.approx(3 * A)
    assert hs["h_mean"] == pytest.approx(2 * A)


# ---------------------------------------------------------------- 3-D smoke
def test_3d_smoke_growth_and_labels():
    g = Grid((9, 9, 12), A)
    eng = GeometryEngine(g)
    eng.seed((4, 4), 1, 1)
    V0 = eng.volume_placed
    for _ in range(6):                # ~3 cells of radius in a 9-wide box
        step(eng, 1e-8, f_max=0.5)
    s = eng.salt_mask()
    assert np.array_equal(s, s[::-1, :, :]) and np.array_equal(s, s[:, ::-1, :])
    assert mo.n_clusters(s, g) == 1
    assert abs(eng.conservation_residual(V0)) < 1e-12
    assert not mo.spanning(s, g)
