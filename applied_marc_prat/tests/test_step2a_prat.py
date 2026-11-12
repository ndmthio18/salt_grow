"""Step 2A acceptance: Prat (2026) eq 29 against the independent V0 reaction-transport reference
(runs/S_surface/step2a_prat/layered_reference.py, loaded by path: review material, deliberately
not a saltv2 module).  One case, d_sub = 112.5 um, eps_new = 0.43, j = 1.6e-5 kg/m2/s.  ~5 s."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from saltv2.growth.prat import prat

_p = Path(__file__).resolve().parents[1] / "runs" / "S_surface" / "step2a_prat" / "layered_reference.py"
_spec = importlib.util.spec_from_file_location("layered_reference", _p)
lr = importlib.util.module_from_spec(_spec); sys.modules["layered_reference"] = lr; _spec.loader.exec_module(lr)

C_EQ, RHO, RHO_CR, EPS, KR, DS, BETA = 0.264, 1200.0, 2170.0, 0.43, 2.3e-3, 1.3e-9, 2.0
J_W, D_SUB = 1.6e-5, 112.5e-6
D_C = 0.14 * D_SUB; U_EVAP = J_W / RHO
V0 = dict(evolve_eps=True, liquid_consumption=True, local_props=False, transport_eps_fixed=True)


@pytest.fixture(scope="module")
def steady():
    return lr.homogeneous_steady(D_SUB, EPS, U_EVAP, lr.PratParams(), n_band=20, **V0)


def test_core_and_reference_analytic_agree_with_the_unrounded_specific_surface():
    av = 5.0 * (1 - EPS) / D_C
    a = lr.prat_analytic(D_C, EPS, U_EVAP)
    assert np.isclose(prat.Da_star(EPS, kr=KR, Ds=DS, beta=BETA, av=av), a["Da_star"], rtol=1e-12, atol=0.0)
    assert np.isclose(prat.F(a["Da_star"]), a["F"], rtol=1e-12, atol=0.0)
    assert np.isclose(prat.boundary_layer_thickness(EPS, kr=KR, Ds=DS, av=av), a["delta_C"], rtol=1e-12, atol=0.0)
    v29 = C_EQ * J_W / (RHO_CR * (1 - EPS) * (1 - C_EQ)) * prat.F(a["Da_star"])
    assert np.isclose(v29, a["v"], rtol=1e-12, atol=0.0)
    # the rounded coefficient 36 of core.specific_surface against 5 / 0.14: documented, below 1e-3 in v
    v36 = prat.growth_rate(J_W, EPS, D_SUB, KR, DS, BETA)
    assert 1e-4 < abs(v36 - v29) / v29 < 1e-3


def test_large_Da_approximation_holds_at_the_initial_height():
    Da = prat.Da(0.3e-3, EPS, D_SUB, KR, DS)
    assert Da > 1e4 and np.tanh(np.sqrt(Da)) == 1.0
    assert np.isclose(prat.growth_rate(J_W, EPS, D_SUB, KR, DS, BETA, h=0.3e-3), prat.growth_rate(J_W, EPS, D_SUB, KR, DS, BETA), rtol=1e-14, atol=0.0)


def test_reference_growth_velocity_and_partition_match_eq_29(steady):
    a = steady["analytic"]
    assert np.isclose(steady["v"], a["v"], rtol=5e-3, atol=0.0)              # verification target 0.5 %; measured ~1e-4
    assert np.isclose(steady["F_num"], a["F"], rtol=5e-3, atol=0.0)
    assert np.isclose(steady["F_prat"], steady["F_num"], rtol=1e-9, atol=0.0)   # v / v_max equals the precipitation share
    assert np.isclose(steady["C_top_star"], a["C_top_star"], rtol=5e-3, atol=0.0)
    prof = lr.prat_profile(steady["xi"], a["delta_C"], a["C_top_star"])
    assert np.max(np.abs((steady["C"] - C_EQ) - prof)) / a["C_top_star"] < 1e-3


def test_balances_and_porosity(steady):
    assert abs(steady["res_salt"]) < 1e-8 and abs(steady["res_solid"]) < 1e-6
    eps = steady["eps"]
    assert np.all(eps > 0.0) and np.all(eps <= EPS + 1e-15)
    closure = 1 - (1 - EPS) / steady["F_num"]
    assert np.isclose(steady["eps_final"], closure, rtol=1e-6, atol=0.0)
    # liquid: the brine volume consumed by eqs 3-4 is the rise of W through the band
    assert np.isclose(steady["W_L"] - steady["W_i"], steady["S_int"] / RHO, rtol=1e-6, atol=0.0)
    # water reaching the interface equals the evaporation flux (eqs 6 + 7); W_i carries the
    # previous Picard iterate of v, so the identity holds to the Picard tolerance (1e-10 on v)
    assert np.isclose(RHO * steady["W_i"] - steady["P_i"], J_W, rtol=1e-8, atol=0.0)
