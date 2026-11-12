"""V1a acceptance tests: fixed-geometry quasi-steady vapour diffusion."""
import numpy as np
import pytest

from saltv2 import (
    GAS,
    SALT_WET,
    Grid,
    VaporProperties,
    solve_vapor,
)


A = 1.0e-4
PROPS = VaporProperties(rho_D=3.0e-5, X_wet=0.0198, X_top=0.00185)


@pytest.mark.parametrize("shape", [(12, 30), (7, 6, 30)])
def test_flat_wet_surface_matches_analytic_1d_flux(shape):
    grid = Grid(shape, A)
    state = np.full(shape, GAS, dtype=np.int8)
    field = solve_vapor(grid, state, PROPS)

    height = shape[-1] * A
    expected = PROPS.rho_D * (PROPS.X_wet - PROPS.X_top) / height
    assert np.all(field.base_flux == pytest.approx(expected, rel=1e-12))
    assert field.M_total / field.projected_area == pytest.approx(expected, rel=1e-12)
    assert field.M_salt == 0.0
    assert abs(field.imbalance) < 1e-12


def test_flat_salt_layer_uses_half_cell_dirichlet_distance():
    grid = Grid((10, 32), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    layer_cells = 5
    state[:, :layer_cells] = SALT_WET

    field = solve_vapor(grid, state, PROPS)
    gap = (grid.n_z - layer_cells) * A
    expected = PROPS.rho_D * (PROPS.X_wet - PROPS.X_top) / gap
    horizontal_top_flux = field.face_flux(grid.z_axis, -1)[:, layer_cells]

    assert np.all(horizontal_top_flux == pytest.approx(expected, rel=1e-12))
    assert field.M_base == 0.0
    assert abs(field.imbalance) < 1e-12


def test_zero_driving_force_gives_zero_flux():
    grid = Grid((16, 20), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    state[7:10, :3] = SALT_WET
    props = VaporProperties(rho_D=3.0e-5, X_wet=0.01, X_top=0.01)

    field = solve_vapor(grid, state, props)
    assert np.nanmax(np.abs(field.X - props.X_wet)) < 1e-14
    assert abs(field.M_total) < 1e-18
    assert abs(field.M_top) < 1e-18


def test_seed_flux_is_globally_conservative():
    grid = Grid((28, 30), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    state[12:16, :4] = SALT_WET

    field = solve_vapor(grid, state, PROPS)
    assert field.M_base > 0.0
    assert field.M_salt > 0.0
    assert field.M_top > 0.0
    assert abs(field.imbalance) < 1e-11


def test_3d_block_has_conservative_vertical_and_lateral_face_fluxes():
    grid = Grid((9, 8, 18), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    state[3:6, 3:5, :3] = SALT_WET

    field = solve_vapor(grid, state, PROPS)
    assert field.M_salt > 0.0
    assert field.face_flux(0, -1).sum() > 0.0
    assert field.face_flux(1, -1).sum() > 0.0
    assert field.face_flux(grid.z_axis, -1).sum() > 0.0
    assert abs(field.imbalance) < 1e-11


def test_periodic_translation_translates_solution_and_fluxes():
    grid = Grid((24, 26), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    state[4:7, :4] = SALT_WET
    shift = 9

    first = solve_vapor(grid, state, PROPS)
    second = solve_vapor(grid, np.roll(state, shift, axis=0), PROPS)

    assert np.allclose(second.X, np.roll(first.X, shift, axis=0), rtol=0, atol=1e-14)
    assert np.allclose(
        second.base_flux, np.roll(first.base_flux, shift, axis=0),
        rtol=0, atol=1e-14,
    )
    for direction, flux in first.face_fluxes.items():
        assert np.allclose(
            second.face_fluxes[direction], np.roll(flux, shift, axis=0),
            rtol=0, atol=1e-14,
        )


def test_sealed_gas_pocket_is_removed_from_linear_system():
    grid = Grid((12, 12), A)
    state = np.full(grid.shape, GAS, dtype=np.int8)
    trapped = (6, 5)
    for cell in ((5, 5), (7, 5), (6, 4), (6, 6)):
        state[cell] = SALT_WET

    field = solve_vapor(grid, state, PROPS)
    assert not field.active_gas[trapped]
    assert np.isnan(field.X[trapped])
    assert np.isfinite(field.X[field.active_gas]).all()
    assert abs(field.imbalance) < 1e-11


def test_two_nearby_patches_screen_each_other():
    grid = Grid((48, 32), A)
    one = np.full(grid.shape, GAS, dtype=np.int8)
    one[16:19, :4] = SALT_WET
    two = one.copy()
    two[26:29, :4] = SALT_WET

    M_one = solve_vapor(grid, one, PROPS).M_salt
    field_two = solve_vapor(grid, two, PROPS)

    assert field_two.M_salt < 2.0 * M_one
    assert abs(field_two.imbalance) < 1e-11
