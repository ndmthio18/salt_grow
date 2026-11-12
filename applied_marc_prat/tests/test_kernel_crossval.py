"""The two conservative network kernels (plan/19 section 6): same balances, different temporal
evaluation of the flux densities.  Neither bitwise nor 1e-12 agreement is required at a finite
dt; what must hold is (1) each closes the salt and water books with the common yardstick,
(2) their difference is first order in dt, (3) they extrapolate to the same limit.
Regime: evaporation-wicking, saturated rigid column, prescribed demand, no salt sink."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "runs" / "S_surface"))
from kernel_crossval import richardson, run_kernel                      # noqa: E402
from saltv2.cases import stage_t_wicking as st                          # noqa: E402

DTS = (60.0, 30.0, 15.0)


@pytest.fixture(scope="module")
def column_runs():
    layer = st.build_network("C", n_side=1, n_layers=10)
    q = st.demand(layer, "uniform")
    return {k: {dt: run_kernel(k, layer, q, dt) for dt in DTS} for k in ("network_step", "conservative_step")}


def test_both_kernels_conserve_salt_and_water_to_round_off(column_runs):
    for kern, runs in column_runs.items():
        for dt, r in runs.items():
            assert abs(r["rel_salt"]) < 1e-9, (kern, dt, r["rel_salt"])
            assert abs(r["rel_water"]) < 1e-8, (kern, dt, r["rel_water"])      # network_step: STORAGE_RTOL bound
            assert r["max_step"]["salt"] < 1e-18


def test_kernels_differ_at_first_order_in_dt_and_share_the_limit(column_runs):
    net, con = column_runs["network_step"], column_runs["conservative_step"]
    d = [abs(net[dt]["t_sat"] - con[dt]["t_sat"]) for dt in DTS]
    assert d[0] > 0.0
    assert 1.8 < d[0] / d[1] < 2.2 and 1.8 < d[1] / d[2] < 2.2          # halves with dt
    dw = [float(np.max(np.abs(net[dt]["w_snap"] - con[dt]["w_snap"]))) for dt in DTS]
    assert 1.8 < dw[0] / dw[1] < 2.2
    r_net = richardson(DTS, [net[dt]["t_sat"] for dt in DTS])
    r_con = richardson(DTS, [con[dt]["t_sat"] for dt in DTS])
    assert 0.9 < r_net["order"] < 1.1 and 0.9 < r_con["order"] < 1.1
    assert abs(r_net["limit"] - r_con["limit"]) < 0.1 * d[2]           # limits agree far below the finite-dt gap
