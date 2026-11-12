"""V1b: bottom plane with bare pore openings (Dirichlet per column) and solid (no-flux)."""
import numpy as np
import pytest

from saltv2 import GAS, SALT_WET, Grid, VaporProperties, solve_vapor

A = 1.0e-4
PROPS = VaporProperties(rho_D=2.86e-5, X_wet=0.01101, X_top=0.00293)


def test_full_mask_reproduces_default():
    grid = Grid((10, 25), A)
    state = np.full(grid.shape, GAS, np.int8); state[4:6, :3] = SALT_WET
    ref = solve_vapor(grid, state, PROPS)
    mask = np.ones(grid.shape[:-1], bool); vals = np.full(grid.shape[:-1], PROPS.X_wet)
    new = solve_vapor(grid, state, PROPS, bottom=(mask, vals))
    assert np.allclose(new.X, ref.X, equal_nan=True, rtol=0, atol=1e-15)
    assert new.M_base == pytest.approx(ref.M_base, rel=1e-12)


def test_single_opening_in_solid_plane_is_conservative_and_smaller():
    grid = Grid((20, 20, 30), A)
    state = np.full(grid.shape, GAS, np.int8)
    mask = np.zeros(grid.shape[:-1], bool); mask[8:12, 8:12] = True
    vals = np.full(grid.shape[:-1], PROPS.X_wet)
    f = solve_vapor(grid, state, PROPS, bottom=(mask, vals))
    plane = solve_vapor(grid, state, PROPS)
    assert f.M_salt == 0.0
    assert 0.0 < f.M_base < plane.M_base
    assert abs(f.imbalance) < 1e-10
    assert np.all(f.base_flux[~mask] == 0.0)
    # an isolated opening evaporates faster per area than the full plane (edge effect)
    assert f.base_flux[mask].mean() > plane.base_flux.mean()


def test_per_column_values_lower_activity_gives_lower_flux():
    grid = Grid((12, 20), A)
    state = np.full(grid.shape, GAS, np.int8)
    mask = np.ones(grid.shape[:-1], bool)
    vals = np.full(grid.shape[:-1], PROPS.X_wet); vals[:6] = 0.8 * PROPS.X_wet
    f = solve_vapor(grid, state, PROPS, bottom=(mask, vals))
    assert f.base_flux[:6].mean() < f.base_flux[6:].mean()
    assert abs(f.imbalance) < 1e-10
