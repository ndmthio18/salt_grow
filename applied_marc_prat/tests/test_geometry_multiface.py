"""P1a / P1b acceptance tests of plan/17 for the geometry engine.

P1a: the growth step limit is the occupancy gain of a candidate CELL, overflow
is redistributed locally (no transfer between disconnected clusters), and an
unsettled overflow raises instead of leaving cells above f = 1.
P1b: the regularisation length ell smooths the deposition rate over a physical
length while conserving the injected volume, and it damps a short-wavelength
perturbation of the front more strongly than a long-wavelength one.
"""
import numpy as np
import pytest

from saltv2 import GAS, SALT_WET, Grid, GeometryEngine
from saltv2.growth.regularize import smooth_on_mask

A = 1e-4


# ------------------------------------------------------------------ P1a item 1
def test_dt_from_cell_occupancy_bounds_a_multiface_cell():
    """A cell in a concave corner is fed by two faces: the face-based limit
    cfl*a/V_n lets it gain 2*cfl occupancy, the cell-based limit gives f_max."""
    g = Grid((9, 9), A)
    eng = GeometryEngine(g)
    eng.state[:4, 0] = SALT_WET          # horizontal arm
    eng.state[0, :4] = SALT_WET          # vertical arm -> cell (1,1) sees two salt faces
    Vn = 1e-8
    rate = eng.volume_rate(Vn)
    corner = (1, 1)
    assert rate[corner] == pytest.approx(2 * Vn * g.face_area)     # two faces
    f_max = 0.5
    dt_cell = eng.dt_occupancy(rate, f_max=f_max)
    dt_face = f_max * g.a / Vn                                     # the old, wrong limit
    assert dt_cell == pytest.approx(0.5 * dt_face)
    eng.advance(rate, dt_cell)
    assert eng.f.max() == pytest.approx(f_max)
    assert eng.state[corner] == GAS                                # did not flip
    # the old limit hands the corner cell 2*f_max = a whole cell in one step
    eng2 = GeometryEngine(g)
    eng2.state[:4, 0] = SALT_WET
    eng2.state[0, :4] = SALT_WET
    rate2 = eng2.volume_rate(Vn)
    assert rate2[corner] * dt_face / g.cell_volume == pytest.approx(2 * f_max)
    eng2.advance(rate2, dt_face)
    assert eng2.state[corner] == SALT_WET                          # overshoot: flipped


def test_every_candidate_gains_at_most_f_max_per_step_in_3d():
    g = Grid((13, 13, 13), A)
    eng = GeometryEngine(g)
    eng.seed((6, 6), height_cells=1, width_cells=1)
    rng = np.random.default_rng(3)
    f_max = 0.4
    for _ in range(30):
        Vn = 1e-8 * (0.2 + rng.random(g.shape))     # strongly non-uniform
        rate = eng.volume_rate(Vn)
        dt = eng.dt_occupancy(rate, f_max=f_max)
        before = eng.f.copy()
        rep = eng.advance(rate, dt)
        gained = (rate * dt / g.cell_volume)
        assert gained.max() <= f_max + 1e-12
        assert rep.dt == pytest.approx(dt)
        del before


# ------------------------------------------------------------------ P1a item 2
def test_overflow_stays_inside_its_own_cluster():
    """Two clusters separated by a wide gas gap.  The fast one overflows by 2.5
    cells; under the old global pooling that volume was shared over every
    receiver in the box, including the slow cluster's candidates.  Now the slow
    cluster may only gain what its own faces deposit."""
    g = Grid((40, 6), A)
    eng = GeometryEngine(g)
    eng.state[5, 0] = SALT_WET           # fast cluster
    eng.state[25, 0] = SALT_WET          # slow cluster, 20 cells away
    V0 = eng.volume_placed
    rate = np.zeros(g.shape)
    rate[5, 1] = 3.5 * g.cell_volume     # flips and overflows by 2.5 cells
    rate[25, 1] = 0.1 * g.cell_volume    # slow cluster's own deposition
    eng.advance(rate, 1.0)
    assert eng.state[5, 1] == SALT_WET
    assert eng.f[25, 1] == pytest.approx(0.1), "overflow leaked to the other cluster"
    assert np.all(eng.f[10:20] == 0.0) and np.all(eng.f[30:] == 0.0)
    assert abs(eng.conservation_residual(V0)) < 1e-12


def test_overflow_is_shared_only_with_own_gas_neighbours():
    """One flipped cell buried on three sides gives all its overflow to the one
    open neighbour; volume is conserved exactly."""
    g = Grid((7, 7), A)
    eng = GeometryEngine(g)
    eng.state[3, 0] = SALT_WET
    eng.state[2, 1] = SALT_WET
    eng.state[4, 1] = SALT_WET
    V0 = eng.volume_placed
    rate = np.zeros(g.shape)
    rate[3, 1] = 1.4 * g.cell_volume        # 1.4 cells worth in one second
    eng.advance(rate, 1.0)
    assert eng.state[3, 1] == SALT_WET
    assert eng.f[3, 2] == pytest.approx(0.4)          # the only gas neighbour
    assert eng.f.sum() == pytest.approx(0.4)
    assert abs(eng.conservation_residual(V0)) < 1e-12


def test_volume_conservation_under_local_redistribution_3d():
    g = Grid((11, 11, 14), A)
    eng = GeometryEngine(g)
    eng.seed((5, 5), 1, 3)
    V0 = eng.volume_placed
    rng = np.random.default_rng(7)
    for _ in range(40):
        Vn = 1e-8 * (0.5 + rng.random())
        rate = eng.volume_rate(Vn)
        eng.advance(rate, eng.dt_occupancy(rate, f_max=0.8))
    assert abs(eng.conservation_residual(V0)) < 1e-12


# ------------------------------------------------------------------ P1a item 3
def test_unsettled_overflow_raises():
    """A narrow gas channel fed with far more volume than one pass can place:
    each pass flips one cell and pushes the rest along, so 30 cells of volume
    cannot settle in 20 passes.  The engine must raise, not silently leave
    cells above f = 1 (plan/17 P1a item 3)."""
    g = Grid((3, 40), A)
    eng = GeometryEngine(g)
    eng.state[:, 0] = SALT_WET
    eng.state[0, :] = SALT_WET            # side walls confine the channel
    eng.state[2, :] = SALT_WET
    rate = np.zeros(g.shape)
    rate[1, 1] = 30.0 * g.cell_volume
    with pytest.raises(RuntimeError, match="did not settle"):
        eng.advance(rate, 1.0)


def test_enclosed_overflow_is_banked_not_lost():
    g = Grid((5, 5), A)
    eng = GeometryEngine(g)
    eng.state[:] = SALT_WET
    eng.state[2, 2] = GAS
    V0 = eng.volume_placed
    rate = np.zeros(g.shape)
    rate[2, 2] = 3.0 * g.cell_volume
    eng.advance(rate, 1.0)
    assert eng.bank == pytest.approx(2.0 * g.cell_volume)
    assert abs(eng.conservation_residual(V0)) < 1e-12


# ------------------------------------------------------------------ P1b
def test_smoothing_conserves_the_injected_volume():
    g = Grid((16, 16, 24), A)
    eng = GeometryEngine(g)
    eng.seed((8, 8), 1, 3)
    rng = np.random.default_rng(11)
    Vn = 1e-8 * rng.random(g.shape)
    raw = eng.volume_rate(Vn)
    for ell in (1.5 * A, 3 * A, 6 * A):
        sm = eng.volume_rate(Vn, smooth_length=ell)
        assert sm.sum() == pytest.approx(raw.sum(), rel=1e-12)
        assert np.all(sm[raw <= 0.0] == 0.0)          # stays on the active set
        assert sm.max() < raw.max()                   # peaks are damped


def test_regularisation_damps_short_wavelengths_more_than_long_ones():
    """The point of ell: on a flat front, a perturbation of wavelength L is
    damped by ~exp(-2 pi^2 ell^2 / L^2), so the short mode loses more."""
    g = Grid((64, 6), A)
    mask = np.zeros(g.shape, bool)
    mask[:, 3] = True                                  # a flat front line
    x = np.arange(64)
    ell = 4 * A
    ratios = []
    for n_periods in (16, 2):                          # L = 4a and 32a
        field = np.zeros(g.shape)
        field[:, 3] = 1.0 + 0.5 * np.cos(2 * np.pi * n_periods * x / 64)
        out = smooth_on_mask(field, mask, g, ell)
        amp_in = 0.5
        amp_out = 0.5 * (out[:, 3].max() - out[:, 3].min())
        ratios.append(amp_out / amp_in)
        assert out.sum() == pytest.approx(field.sum(), rel=1e-12)
    short, long = ratios
    assert short < 0.05 and long > 0.5
    # quantitative: the damping of a cosine by a Gaussian of width ell
    L_long = 64 / 2 * A
    assert long == pytest.approx(np.exp(-2 * np.pi ** 2 * ell ** 2 / L_long ** 2), rel=0.02)


def test_smoothing_does_not_grow_into_forbidden_cells():
    """Cells the caller excludes (e.g. above the supply-limited height) receive
    nothing, and the volume they would have taken stays on the allowed set."""
    g = Grid((10, 10), A)
    eng = GeometryEngine(g)
    eng.state[:, 0] = SALT_WET
    allowed = np.zeros(g.shape, bool)
    allowed[:, 1] = True                               # only the first gas row may grow
    Vn = 1e-8
    raw = eng.volume_rate(Vn)
    sm = eng.volume_rate(Vn, smooth_length=3 * A, smooth_mask=allowed)
    assert np.all(sm[:, 2:] == 0.0)
    assert sm.sum() == pytest.approx(raw.sum(), rel=1e-12)


def test_regularisation_damps_a_needle_tip():
    """The Laplacian instability acts at convex tips, whose top candidate cell is
    only diagonally adjacent to the flank candidates.  The kernel must therefore
    reach across that corner: a face-connected graph kernel on the candidate set
    leaves the tip untouched (found 2026-09-06, README E1 section 9.13)."""
    g = Grid((21, 12), A)
    eng = GeometryEngine(g)
    eng.state[10, :5] = SALT_WET                       # one-cell-wide column
    def velocity(ax, s, mask):
        v = np.full(g.shape, 1e-8)
        v[10, 5] = 1e-7                                # tip flux 10x the flanks
        return v
    raw = eng.volume_rate(velocity)
    for ell in (A, 2 * A, 3 * A):
        sm = eng.volume_rate(velocity, smooth_length=ell)
        assert sm.sum() == pytest.approx(raw.sum(), rel=1e-12)
        assert sm[10, 5] < 0.6 * raw[10, 5]            # tip damped
        assert sm[9, 4] > raw[9, 4]                    # flanks receive it
