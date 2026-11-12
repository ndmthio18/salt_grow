import numpy as np
import pytest

from saltv2.cases import dumux5x5_x008 as case


def test_the_short_case_reconstructs_all_original_pores_and_throats():
    arrays, summary = case.run(dt=5.0, t_end=20.0, output_interval=10.0)

    assert arrays["time_s"].tolist() == [0.0, 10.0, 20.0]
    assert arrays["pore_x"].shape == (3, 30)
    assert arrays["pore_pressure_gauge_pa"].shape == (3, 30)
    assert arrays["throat_flow_m3_s"].shape == (3, 45)
    assert arrays["throat_radius_m"].shape == (3, 45)
    assert np.all(arrays["pore_x"][:, [0, 5, 10, 15, 20]] == pytest.approx(case.X_RES))
    assert np.all(arrays["pore_pressure_gauge_pa"][:, [0, 5, 10, 15, 20]] == 0.0)
    assert summary["actual_t_end_s"] == 20.0 and summary["stop_reason"] == "t_end"
    assert abs(summary["salt_balance_relative_residual"]) < 1e-10


def test_flux_replay_bypasses_the_python_boundary_condition():
    times = np.array([0.0, 20.0])
    imposed = np.arange(1.0, 6.0)[None, :] * 1e-13
    replay = (times, np.vstack([imposed, imposed]))
    arrays, summary = case.run(dt=5.0, t_end=10.0, output_interval=5.0,
                               precipitation=False, replay_flux=replay)

    assert np.all(arrays["top_evaporation_kg_s"][1:] == imposed)
    assert summary["boundary_mode"] == "dumux_flux_replay"
    assert summary["final_solid_mass_kg"] == 0.0

