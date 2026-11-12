"""V1c: iterative solver and vectorised connectivity reproduce the reference."""
import numpy as np
import pytest

from saltv2 import GAS, SALT_WET, Grid, VaporProperties, solve_vapor
from saltv2.growth.fv import top_connected

PROPS = VaporProperties(rho_D=2.86e-5, X_wet=0.01101, X_top=0.00293)


@pytest.mark.parametrize("shape", [(20, 24), (9, 9, 20)])
def test_cg_matches_direct(shape):
    g = Grid(shape, 1e-4)
    st = np.full(g.shape, GAS, np.int8)
    st[(slice(3, 6),) * (g.ndim - 1) + (slice(0, 4),)] = SALT_WET
    d = solve_vapor(g, st, PROPS)
    c = solve_vapor(g, st, PROPS, solver="cg", rtol=1e-12)
    assert np.allclose(c.X, d.X, equal_nan=True, rtol=1e-8, atol=1e-13)
    assert c.M_total == pytest.approx(d.M_total, rel=1e-8)
    w = solve_vapor(g, st, PROPS, solver="cg", x0=d.X, rtol=1e-12)
    assert w.M_total == pytest.approx(d.M_total, rel=1e-8)


def top_connected_bfs(grid, open_mask):
    """Reference implementation: a plain breadth-first search from the top layer over the
    grid's periodic lateral topology.  Lives in the test on purpose -- it exists only to check
    the vectorised ``fv.top_connected``."""
    from collections import deque
    open_mask = np.asarray(open_mask, dtype=bool)
    active = np.zeros(grid.shape, dtype=bool)
    queue = deque()
    for lateral in np.argwhere(open_mask[..., -1]):
        cell = tuple(int(i) for i in lateral) + (grid.n_z - 1,)
        active[cell] = True
        queue.append(cell)
    while queue:
        cell = queue.popleft()
        for axis, step in grid.face_directions():
            neighbour = list(cell)
            neighbour[axis] += step
            if axis in grid.lateral_axes:
                neighbour[axis] %= grid.shape[axis]
            elif neighbour[axis] < 0 or neighbour[axis] >= grid.shape[axis]:
                continue
            neighbour = tuple(neighbour)
            if open_mask[neighbour] and not active[neighbour]:
                active[neighbour] = True
                queue.append(neighbour)
    return active


def test_top_connected_matches_bfs_with_pockets_and_seams():
    g = Grid((14, 10, 12), 1e-4)
    rng = np.random.default_rng(1)
    st = np.full(g.shape, GAS, np.int8)
    st[rng.random(g.shape) < 0.35] = SALT_WET
    st[..., -1] = GAS
    a = top_connected(g, st == GAS); b = top_connected_bfs(g, st == GAS)
    assert np.array_equal(a, b)


def test_external_length_flat_surface_is_series_resistance():
    g = Grid((6, 6, 10), 1e-4)
    st = np.full(g.shape, GAS, np.int8)
    H = g.n_z * g.a; L_ext = 3.0 * H
    f = solve_vapor(g, st, PROPS, top_external_length=L_ext)
    expected = PROPS.rho_D * (PROPS.X_wet - PROPS.X_top) / (H + L_ext)
    assert np.all(f.base_flux == pytest.approx(expected, rel=1e-12))
    assert f.M_top == pytest.approx(f.M_base, rel=1e-11)
    same = solve_vapor(g, st, PROPS, top_external_length=0.0)
    ref = solve_vapor(g, st, PROPS)
    assert np.allclose(same.X, ref.X, equal_nan=True, rtol=0, atol=1e-15)


def test_external_length_reproduces_tall_dirichlet_box_above_a_structure():
    """A short box with L_ext = delta - z_box must match the tall box with
    Dirichlet at delta when the box top is >= one lateral period above the salt."""
    a = 1e-4; n_lat = 12; h_salt = 6; z_box = h_salt + n_lat; z_tall = 4 * n_lat
    tall = Grid((n_lat, n_lat, z_tall), a); short = Grid((n_lat, n_lat, z_box), a)
    st_tall = np.full(tall.shape, GAS, np.int8); st_tall[4:8, 4:8, :h_salt] = SALT_WET
    st_short = st_tall[:, :, :z_box].copy()
    ft = solve_vapor(tall, st_tall, PROPS)
    fs = solve_vapor(short, st_short, PROPS, top_external_length=(z_tall - z_box) * a)
    assert fs.M_total == pytest.approx(ft.M_total, rel=1e-6)
    for key in ft.face_fluxes:
        assert np.allclose(fs.face_fluxes[key], ft.face_fluxes[key][:, :, :z_box], rtol=1e-5, atol=1e-12)
    assert abs(fs.imbalance) < 1e-11


def test_porous_cap_is_series_resistance_and_default_unchanged():
    """A porous (sealed) salt cap of thickness L_d and factor f over a wet layer: flux =
    rho_D dX / (L_gas + L_d / f). porous=None reproduces the gas-only solution bit for bit."""
    from saltv2 import SALT_SEALED
    g = Grid((6, 6, 20), 1e-4)
    st = np.full(g.shape, GAS, np.int8)
    st[:, :, :3] = SALT_WET               # wet salt layer
    st[:, :, 3:8] = SALT_SEALED           # sealed cap, 5 cells
    f = 0.1
    fld = solve_vapor(g, st, PROPS, porous=(st == SALT_SEALED, f))
    L_gas = (g.n_z - 8) * g.a; L_d = 5 * g.a
    expected = PROPS.rho_D * (PROPS.X_wet - PROPS.X_top) / (L_gas + L_d / f)
    front = fld.face_flux(g.z_axis, -1)[:, :, 3]          # sealed cells just above the wet layer
    assert np.all(front == pytest.approx(expected, rel=1e-10))
    assert fld.M_top == pytest.approx(fld.M_salt, rel=1e-10)
    assert np.isfinite(fld.X[st == SALT_SEALED]).all()
    # without the porous argument the sealed cells are walls: no evaporation at all
    walls = solve_vapor(g, st, PROPS)
    assert walls.M_salt == 0.0
    # gas-only problem unchanged by porous=None versus an empty mask
    st2 = np.full(g.shape, GAS, np.int8); st2[2:4, 2:4, :3] = SALT_WET
    a = solve_vapor(g, st2, PROPS); b = solve_vapor(g, st2, PROPS, porous=(np.zeros(g.shape, bool), 1.0))
    assert np.array_equal(np.nan_to_num(a.X), np.nan_to_num(b.X)) and a.M_total == b.M_total


def test_adaptive_box_matches_the_full_box():
    """solve_vapor_adaptive cuts the box one lateral period above the salt and puts the removed
    height into the external segment: fluxes and total rate must match the full-box solve."""
    import numpy as np
    from saltv2 import Grid, GeometryEngine, VaporProperties
    from saltv2.growth.vapor import solve_vapor, solve_vapor_adaptive
    g = Grid((12, 12, 160), 1e-5)
    eng = GeometryEngine(g)
    eng.seed((3, 3), height_cells=9, width_cells=3); eng.seed((8, 7), height_cells=4, width_cells=2)
    eng.state[:, :, 0] = 1
    props = VaporProperties(rho_D=3e-5, X_wet=0.011, X_top=0.003)
    full = solve_vapor(g, eng.state, props, solver="cg", rtol=1e-10, top_external_length=5e-3)
    adap = solve_vapor_adaptive(g, eng.state, props, solver="cg", rtol=1e-10, top_external_length=5e-3)
    assert abs(adap.M_total - full.M_total) / full.M_total < 1e-6
    for key in full.face_fluxes:
        assert np.allclose(adap.face_fluxes[key], full.face_fluxes[key], rtol=1e-5, atol=1e-9 * full.face_fluxes[key].max())
    assert adap.active_gas[..., -1].all()
    assert np.isfinite(adap.X[..., -1]).all()

