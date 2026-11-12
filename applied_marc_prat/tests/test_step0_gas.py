"""Step 0 acceptance (verify/step0): gas-phase evaporation from one square surface pore,
saltv2.cases.step0_gas.  Gas only; seconds."""
import numpy as np
import pytest

from saltv2.cases import step0_gas as g


@pytest.fixture(scope="module")
def results():
    out = {}
    for case, h in g.H_BOX.items():
        out[("plane", case)] = g.solve_plane(h)[0]
        out[("pore", case)] = g.solve_pore(h)[0]
    return out


@pytest.mark.parametrize("case", list(g.H_BOX))
def test_V01_fully_wet_plane_reproduces_the_1d_solution(results, case):
    p = results[("plane", case)]
    assert abs(p["rel_err_1d"]) < 1e-9
    assert abs(p["imbalance"]) < 1e-9
    assert p["discrete_residual"] < 1e-12
    # the profile is the linear 1-D one at every cell centre
    z = (np.arange(len(p["X_profile"])) + 0.5) * p["a_voxel"]
    X_1d = p["X_wet"] - (p["X_wet"] - p["X_top"]) * z / g.DELTA_EFF
    assert np.allclose(p["X_profile"], X_1d, rtol=0.0, atol=1e-12 * p["X_wet"])


@pytest.mark.parametrize("case", list(g.H_BOX))
def test_V02_single_pore_vapour_mass_is_conserved(results, case):
    q = results[("pore", case)]
    assert abs(q["imbalance"]) < 1e-9
    assert q["discrete_residual"] < 1e-12


def test_V03_single_pore_rate_is_independent_of_the_numerical_domain_height(results):
    qa, qb = results[("pore", "A")], results[("pore", "B")]
    # The 1-D external segment is exact once X is laterally uniform at the box top; the slowest
    # lateral mode decays as exp(-2 pi z / L_box), 1.6e-8 at 2 mm.  Tolerance from that bound.
    bound = float(np.exp(-2.0 * np.pi * g.H_BOX["A"] / g.L_BOX))
    assert qa["top_nonuniformity"] < bound
    for key in ("E", "j_projected", "j_wet"):
        assert abs(qa[key] - qb[key]) / qb[key] < 10.0 * bound


def test_V04_rasterised_opening_area_and_physical_bounds(results):
    ras = g.rasterisation()
    assert ras["n_cells"] == 100 and ras["side_raster"] == pytest.approx(500e-6)
    assert ras["area_rel_err"] == pytest.approx(0.25 / 0.2916 - 1.0, abs=1e-12)
    for case in g.H_BOX:
        p, q = results[("plane", case)], results[("pore", case)]
        # an isolated wet patch evaporates less than the plane, more than its area share of it
        assert p["E"] * q["A_wet"] / q["A_box"] < q["E"] < p["E"]
        assert q["L_spread"] > 0.0
        assert q["base_flux_max_over_mean"] > 1.0        # edge enhancement of the isolated patch
