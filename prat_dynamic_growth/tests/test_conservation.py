import numpy as np
import pytest

from prat_growth.dynamic import DynamicALEModel, snapshot_diagnostics
from prat_growth.parameters import NumericalParameters, PhysicalParameters


@pytest.mark.parametrize("closure", ["B1", "B2"])
def test_short_run_global_conservation(closure):
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=80, t_end=30.0, output_interval=30.0, max_step=2.0)
    result = DynamicALEModel(p, n, closure).run()
    diagnostic = snapshot_diagnostics(result, -1)
    assert result.status == 0
    assert abs(diagnostic["salt_residual_normalized"]) < 1e-8
    assert abs(diagnostic["bulk_residual_normalized"]) < 1e-8


def test_reaction_law_allows_dissolution_without_clipping():
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=8, t_end=1.0, output_interval=1.0)
    model = DynamicALEModel(p, n, "B2")
    concentration = np.full(n.n_cells, p.C_eq - 1e-4)
    epsilon = np.full(n.n_cells, p.epsilon_ref)
    state = model.flux_state(epsilon, concentration, p.h0)
    assert np.all(state["reaction"] < 0.0)
