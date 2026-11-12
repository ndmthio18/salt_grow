"""Acceptance tests of the independent 1-D reference (run from this directory, see README.md).

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
    cd reviews/prat_layered_validation && $PY -m pytest test_layered_reference.py -q      # ~1 min

Verified mathematical results (the reference reproduces Prat's own linearised problem) are tested
tightly; the physical variants are tested only for the sign and order of their departure, which is
what the report claims.
"""
import numpy as np
import pytest

from layered_reference import (DAY, Layering, PratParams, ReferenceSolver, eps_new_for_final,
                               homogeneous_steady, prat_analytic, prat_profile)

P = PratParams()
U = 1e-2 / DAY
HYBRID = dict(evolve_eps=True, local_props=False, transport_eps_fixed=True)
CONSISTENT = dict(evolve_eps=True, local_props=True)


@pytest.mark.parametrize("d_sub", [10e-6, 25e-6, 112.5e-6])
def test_hybrid_reproduces_prat_eq29_eq31_and_the_birth_porosity_closure(d_sub):
    """With eps held at its birth value inside the transport operator (Prat's linearisation) and
    the liquid balance of eqs 3-4 kept, the numerical interface speed, top supersaturation and
    final porosity match eqs 29, 31 and 1 - (1 - eps_s)/F to the discretisation error."""
    a = prat_analytic(P.dc_over_dsub * d_sub, 0.43, U, P)
    r = homogeneous_steady(d_sub, 0.43, U, P, n_band=40, **HYBRID)
    assert r["v"] == pytest.approx(a["v"], rel=3e-4)
    assert r["C_top_star"] == pytest.approx(a["C_top_star"], rel=3e-4)
    assert r["F_prat"] == pytest.approx(a["F"], rel=3e-4)
    assert r["eps_final"] == pytest.approx(1 - (1 - 0.43) / a["F"], abs=2e-4)
    assert abs(r["res_salt"]) < 1e-8 and abs(r["res_solid"]) < 1e-8


def test_profile_matches_eq19_in_the_band():
    a = prat_analytic(P.dc_over_dsub * 10e-6, 0.43, U, P)
    r = homogeneous_steady(10e-6, 0.43, U, P, n_band=80, **HYBRID)
    err = np.max(np.abs((r["C"] - P.C_eq) - prat_profile(r["xi"], a["delta_C"], a["C_top_star"]))) / a["C_top_star"]
    assert err < 2e-3


def test_steady_state_is_independent_of_dt_and_window():
    ref = homogeneous_steady(10e-6, 0.43, U, P, n_band=40, **CONSISTENT)
    for kw in (dict(dt_over_transit=0.1), dict(dt_over_transit=2.0), dict(L_over_deltaC=120.0)):
        r = homogeneous_steady(10e-6, 0.43, U, P, n_band=40, **CONSISTENT, **kw)
        assert r["v"] == pytest.approx(ref["v"], rel=1e-9)
        assert r["eps_final"] == pytest.approx(ref["eps_final"], rel=1e-9)


def test_first_order_grid_convergence():
    a = prat_analytic(P.dc_over_dsub * 10e-6, 0.43, U, P)
    errs = [abs(homogeneous_steady(10e-6, 0.43, U, P, n_band=nb, **HYBRID)["v"] / a["v"] - 1) for nb in (10, 20, 40)]
    assert errs[0] > errs[1] > errs[2]
    assert errs[0] / errs[1] > 1.8                                  # ~ first order while above the modelling plateau


def test_constant_porosity_model_does_not_reproduce_eq29():
    """Prat's eq 5 carries the (1 - C) factor that comes from the pore volume lost to precipitation
    (eq 3 with eq 2).  Holding eps constant in every equation therefore gives a SLOWER front than
    eq 29: the law presupposes the porosity reduction it says does not happen."""
    for d_sub, lo, hi in ((10e-6, 0.02, 0.05), (1080e-6, 0.002, 0.01)):
        a = prat_analytic(P.dc_over_dsub * d_sub, 0.43, U, P)
        r = homogeneous_steady(d_sub, 0.43, U, P, n_band=40, evolve_eps=False, local_props=False)
        assert lo < 1 - r["v"] / a["v"] < hi


def test_dropping_liquid_consumption_removes_the_one_minus_Ceq_factor():
    """Without eqs 3-4 the band sink acts on C with coefficient 1 instead of (1 - C_eq): delta_C
    shrinks by sqrt(1 - C_eq) and F falls accordingly (interface share / (share + a_v delta_C))."""
    a = prat_analytic(P.dc_over_dsub * 25e-6, 0.43, U, P)
    r = homogeneous_steady(25e-6, 0.43, U, P, n_band=40, liquid_consumption=False, **HYBRID)
    x = np.sqrt(a["Da_star"])                                    # beta(1-eps)/(a_v delta_C)
    F_expected = x / (x + 1 / np.sqrt(1 - P.C_eq))
    assert r["F_prat"] == pytest.approx(F_expected, rel=2e-3)


def test_consistent_model_is_faster_and_denser_than_eq29_and_converges():
    """Local a_v(eps) and D_eff(eps) in the densified band raise Da* locally: v above eq 29 and
    eps_final above the closure value; the departure shrinks with d_sub."""
    dep = []
    for d_sub in (10e-6, 112.5e-6, 1080e-6):
        a = prat_analytic(P.dc_over_dsub * d_sub, 0.43, U, P)
        r = homogeneous_steady(d_sub, 0.43, U, P, n_band=40, **CONSISTENT)
        dep.append(r["v"] / a["v"] - 1)
        assert r["eps_final"] > 1 - (1 - 0.43) / a["F"]
    assert dep[0] > dep[1] > dep[2] > 0
    assert 0.04 < dep[0] < 0.07
    v40 = homogeneous_steady(10e-6, 0.43, U, P, n_band=40, **CONSISTENT)["v"]
    v80 = homogeneous_steady(10e-6, 0.43, U, P, n_band=80, **CONSISTENT)["v"]
    assert abs(v80 / v40 - 1) < 1e-3


def test_eps_final_interpretation_is_always_feasible():
    """A birth porosity exists for every target final porosity: 1 - eps_new = F (1 - eps_final)."""
    e, hist, r = eps_new_for_final(10e-6, 0.43, U, P, n_band=40, **CONSISTENT)
    assert r["eps_final"] == pytest.approx(0.43, abs=1e-6)
    assert 0.43 < e < 1.0
    assert r["eps_final"] == pytest.approx(1 - (1 - e) / r["F_num"], abs=1e-6)


def test_layered_run_closes_the_pores_at_the_predicted_layer_and_conserves_salt():
    """27.5 um, r = 0.5, alpha = 2, c_h = 1, eps_new = 0.43.  The layer-wise Prat F crosses the crust
    limit 1 - eps_s at layer 4 (F_4 = 0.50): the hybrid reference (Prat's linearisation) closes its
    band porosity there; the consistent reference (local a_v, D_eff in the densified band) has a
    higher F on every layer and closes one layer later.  Run-integrated balances close either way."""
    lay = Layering(d_sub=27.5e-6, r=0.5, alpha=2.0, c_h=1.0)
    a1 = prat_analytic(lay.d_c(1, P), 0.43, U, P); a5 = prat_analytic(lay.d_c(5, P), 0.43, U, P)
    F = [prat_analytic(lay.d_c(n, P), 0.43, U, P)["F"] for n in range(1, 6)]
    assert F[2] > 1 - 0.43 > F[3]
    closed = {}
    for lab, kw in (("hybrid", HYBRID), ("consistent", CONSISTENT)):
        s = ReferenceSolver(lay, U, 0.43, P, L=60 * a1["delta_C"], dxi0=a5["delta_C"] / 30, **kw)
        s.run(t_end=1e9, dt=0.5 * a5["delta_C"] / a5["v"], h_end=lay.layer_bottom(7), record_every=10)
        assert s.closure_event is not None
        closed[lab] = s.closure_event["layer"]
        c = s.cumulative_balance()
        assert abs(c["res_salt"]) < 1e-8 and abs(c["res_solid"]) < 1e-6
        # closure happens within one band thickness of the layer's start
        assert 0 < s.closure_event["h"] - lay.layer_bottom(closed[lab]) < 3 * a1["delta_C"]
        # the quasi-steady F of every completed layer follows the layer-wise Prat value (hybrid) or sits above it
        Fp = np.array([q["F_prat"] for q in s.history]); n = np.array([q["layer_top"] for q in s.history])
        for k in range(1, closed[lab]):
            tail = np.where(n == k)[0][-20:]
            if lab == "hybrid":
                assert np.mean(Fp[tail]) == pytest.approx(F[k - 1], rel=5e-3)
            else:
                assert np.mean(Fp[tail]) > F[k - 1]
    assert closed["hybrid"] == 4 and closed["consistent"] == 5


def test_eps_final_interpretation_has_its_own_limit_in_the_consistent_model():
    """Hybrid (Prat's linearisation): a birth porosity reaching eps_final = 0.43 exists for every
    crystal size, as 1 - eps_new = F (1 - eps_final) promises.  Consistent model (a_v and D_eff at
    the local porosity): on the 0.24 um crystals of layer 5 (27.5 um beads) NO birth porosity in
    [0.43, 0.99] reaches 0.43 -- the a_v ~ (1 - eps) feedback closes the band -- so the crust
    limit survives the eps_final reading there."""
    d_eff = 27.5e-6 * 0.5 ** 4                                  # layer 5 crystal size through d_c = 0.14 d_sub
    e_h, hist_h, r_h = eps_new_for_final(d_eff, 0.43, U, P, n_band=40, **HYBRID)
    assert e_h is not None and r_h["eps_final"] == pytest.approx(0.43, abs=1e-5)
    a = prat_analytic(P.dc_over_dsub * d_eff, e_h, U, P)
    assert e_h == pytest.approx(1 - a["F"] * (1 - 0.43), abs=2e-3)   # the analytic identity at the found eps_new
    e_c, hist_c, r_c = eps_new_for_final(d_eff, 0.43, U, P, n_band=40, **CONSISTENT)
    assert e_c is None
    assert max(h[1] for h in hist_c) < 0.43                      # the scan's best final porosity stays below the target
