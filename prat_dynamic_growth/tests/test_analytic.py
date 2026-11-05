from prat_growth.analytic import solve_stage_a_finite_volume
from prat_growth.parameters import PhysicalParameters


def test_stage_a_reproduces_general_growth_law():
    result = solve_stage_a_finite_volume(
        PhysicalParameters(), height=1e-3, n_cells=1280, grid_power=4.0
    )
    assert result.relative_velocity_error < 1e-5


def test_large_da_limit_matches_general_baseline():
    result = solve_stage_a_finite_volume(
        PhysicalParameters(), height=1e-3, n_cells=320, grid_power=4.0
    )
    rel = abs(result.velocity_exact - result.velocity_large_da) / result.velocity_exact
    assert rel < 1e-12
