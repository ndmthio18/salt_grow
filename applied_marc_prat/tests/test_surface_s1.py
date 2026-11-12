"""S1 acceptance (plan/19 section 4.1): the single-pore case of saltv2.cases.surface_s1.
Regime: evaporation-wicking.  ~40 s."""
import numpy as np
import pytest

from saltv2.cases import surface_s1 as s1
from saltv2.materials.brine import C_SAT

T_SAT_D50_PRESCRIBED_H = 1.917636      # network_step, dt 30 s: runs/S_surface/results/kernel_crossval.md [V]


def test_vapour_solver_reproduces_the_1d_flux_on_a_fully_wet_plane():
    pc = s1.plane_check()
    assert abs(pc["rel_err"]) < 1e-9 and abs(pc["imbalance"]) < 1e-9
    assert pc["j_1d"] > 0.0


def test_prescribed_demand_reference_column_reproduces_the_kernel_anchor():
    summ, ser = s1.run(n_layers=50, feedback=False, growth=False, t_end_h=3.0, record_s=30.0)
    assert summ["stop"] == "first saturation"
    assert s1.t_sat_of(ser) / 3600.0 == pytest.approx(T_SAT_D50_PRESCRIBED_H, abs=1e-6)
    assert abs(summ["res_salt"]) < 1e-13 * summ["S_scale"]      # ~200 steps of round-off
    assert abs(summ["res_water"]) < 1e-10 * summ["L_scale"]


def test_feedback_lowers_the_flux_as_the_opening_enriches_and_books_close():
    summ, ser = s1.run(n_layers=10, feedback=True, growth=False, t_end_h=3.0, record_s=60.0)
    assert summ["stop"] == "first saturation"
    assert np.all(np.diff(ser["E"]) <= 1e-18)                  # a_w(C) falls monotonically
    assert ser["E"][-1] / ser["E"][0] == pytest.approx(0.97, abs=0.01)   # a_w 0.78 -> 0.753
    assert abs(summ["res_salt"]) < 1e-13 * summ["S_scale"]
    assert abs(summ["res_water"]) < 1e-10 * summ["L_scale"]
    assert abs(summ["res_mass"]) < 1e-10 * summ["L_scale"]


@pytest.fixture(scope="module")
def growth_saturated():
    return s1.run(n_layers=10, c0=C_SAT, feedback=True, growth=True, t_end_h=1.5, record_s=60.0)


def test_seed_is_sourced_and_growth_matches_prat_eq29_with_a_saturated_reservoir(growth_saturated):
    summ, ser = growth_saturated
    assert summ["t_nuc"] == summ["dt"]                          # at C_sat from the first step
    assert summ["m_seed"] is not None and summ["m_seed"] > 0.0
    # every crystal gram came out of the pore: solid = X delivered - the salt of the fill brine
    assert summ["M_front"] + summ["M_int"] == pytest.approx(summ["cum"]["X"] - C_SAT * ser["M_liq"][-1], rel=1e-12)
    # the record at the second step carries the seed (the salt beyond C_sat of the un-held first
    # step); from then on inflow is at C_sat and dh/dt is eq 29 to round-off
    ok = (ser["v_prat"] > 0.0) & (ser["t"] > 2 * summ["dt"])
    assert ok.sum() > 10 and np.allclose(ser["v"][ok], ser["v_prat"][ok], rtol=1e-9)
    assert summ["M_int"] / (summ["M_front"] + summ["M_int"]) == pytest.approx(1.0 - summ["F"], rel=1e-12)
    assert summ["h_end"] > 0.0


def test_four_books_close_globally_during_growth(growth_saturated):
    summ, ser = growth_saturated
    assert abs(summ["res_salt"]) < 1e-13 * summ["S_scale"]
    assert abs(summ["res_mass"]) < 1e-10 * summ["L_scale"]
    assert abs(summ["res_water"]) < 1e-10 * summ["L_scale"]
    assert np.max(np.abs(ser["res_salt"])) < 1e-13 * summ["S_scale"]
