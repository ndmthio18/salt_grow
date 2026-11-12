"""Acceptance tests of the beta / hierarchy review (run from this directory, ~2 min).

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
    cd reviews/prat_beta_hierarchy && OMP_NUM_THREADS=1 $PY -m pytest test_beta_hierarchy.py -q
"""
import sys
from pathlib import Path

import numpy as np
import pytest

from common import (CONSISTENT, HYBRID, P0, U_EVAP, closure_layer_analytic, d_c_crit, homogeneous_steady, layer_bottom,
                    layer_plateaus, layered_run, prat_analytic, replace, transition_duration)
from beta_analysis import beta_of_cell

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "prat_layered_validation"))
import layered_reference as original                                   # noqa: E402  (the untouched solver)


# ------------------------------------------------------------------ the extended copy is the original
def test_extended_solver_reproduces_the_original_bit_for_bit():
    for kw in (HYBRID, CONSISTENT):
        a = original.homogeneous_steady(25e-6, 0.43, U_EVAP, original.PratParams(), n_band=40, **kw)
        b = homogeneous_steady(25e-6, 0.43, U_EVAP, P0, n_band=40, **kw)
        assert a["v"] == b["v"] and a["eps_final"] == b["eps_final"] and a["C_top_star"] == b["C_top_star"]
        assert np.array_equal(a["C"], b["C"]) and np.array_equal(a["eps"], b["eps"])


# ------------------------------------------------------------------ beta: derivation
def test_beta_two_follows_from_the_appendix_cell_and_is_dilation_invariant():
    A = lambda d: 24 * d ** 2; V = lambda d: 4 * d ** 3                   # eqs A1, A2
    assert beta_of_cell(A, V, 1.0) == pytest.approx(2.0, rel=1e-8)
    assert beta_of_cell(A, V, 123.4) == pytest.approx(2.0, rel=1e-8)   # dilation
    assert beta_of_cell(lambda d: 6 * d ** 2, lambda d: d ** 3) == pytest.approx(2.0, rel=1e-8)       # one free cube
    assert beta_of_cell(lambda d: 5 * d ** 2, lambda d: d ** 3) == pytest.approx(5 / 3, rel=1e-8)     # base attached
    # a needle whose length L sets h: beta = A/(dV/dL) = 4 L/d + 2 (42 at L = 10 d); uniform dilation of the same needle gives 1.4
    assert beta_of_cell(lambda d: 4 * d * 10 * d + 2 * d ** 2, lambda d: 10 * d ** 3) == pytest.approx(1.4, rel=1e-8)


def test_beta_enters_only_through_Da_star_and_the_interface_sink():
    """Da* ~ beta^2; dlnF/dlnbeta = 1 - F; d_c,crit ~ beta^-2."""
    for d_sub in (27.5e-6, 1080e-6):
        a2 = prat_analytic(0.14 * d_sub, 0.43, U_EVAP, P0); a4 = prat_analytic(0.14 * d_sub, 0.43, U_EVAP, replace(P0, beta=4.0))
        assert a4["Da_star"] == pytest.approx(4 * a2["Da_star"], rel=1e-12)
        b = 2.0; db = 1e-4
        Fp = prat_analytic(0.14 * d_sub, 0.43, U_EVAP, replace(P0, beta=b + db))["F"]; Fm = prat_analytic(0.14 * d_sub, 0.43, U_EVAP, replace(P0, beta=b - db))["F"]
        elasticity = (np.log(Fp) - np.log(Fm)) / (np.log(b + db) - np.log(b - db))
        assert elasticity == pytest.approx(1 - a2["F"], rel=1e-5)
    assert d_c_crit(0.43, replace(P0, beta=4.0)) == pytest.approx(d_c_crit(0.43, P0) / 4, rel=1e-12)


def test_closure_layer_shifts_by_one_per_sqrt2_in_beta_in_the_linearised_reference():
    """r = 0.5: n* = 2 + floor(log2(0.14 d_sub / d_c,crit)) and d_c,crit ~ beta^-2, so beta x sqrt(2) adds one layer.
    The hybrid reference closes at the analytic n* for each beta."""
    for b, n_expected in ((2 ** 0.5, 3), (2.0, 4), (2 ** 1.5, 5)):
        p = replace(P0, beta=float(b))
        assert closure_layer_analytic(27.5e-6, 0.43, 0.5, p) == n_expected
        run = layered_run(27.5e-6, 0.43, p=p, variant=HYBRID, n_layers=6, record_every=10)
        assert run["closure"] is not None and run["closure"]["layer"] == n_expected


def test_numerical_beta_sweep_matches_eq29_for_every_beta():
    for b in (1.0, 2.0, 4.0):
        p = replace(P0, beta=float(b)); a = prat_analytic(0.14 * 112.5e-6, 0.43, U_EVAP, p)
        r = homogeneous_steady(112.5e-6, 0.43, U_EVAP, p, n_band=40, **HYBRID)
        assert r["v"] == pytest.approx(a["v"], rel=3e-4) and r["C_top_star"] == pytest.approx(a["C_top_star"], rel=3e-4)


# ------------------------------------------------------------------ hierarchy
def test_alpha_and_c_h_do_not_change_the_closure_layer_only_its_height():
    ref = layered_run(27.5e-6, 0.43, n_layers=5, record_every=10)
    assert ref["closure"] is not None
    n_c = ref["closure"]["layer"]
    for kw in (dict(alpha=1.5), dict(alpha=3.0), dict(c_h=0.5), dict(c_h=2.0)):
        run = layered_run(27.5e-6, 0.43, n_layers=5, record_every=10, **kw)
        assert run["closure"]["layer"] == n_c
        hb = layer_bottom(n_c, 27.5e-6, kw.get("alpha", 2.0), kw.get("c_h", 1.0))
        assert 0 < run["closure"]["h"] - hb < 3 * ref["delta_C"][0]                 # within a band of the layer bottom
        for n in range(1, n_c):                                                       # same plateau F per layer
            assert layer_plateaus(run)[n]["F"] == pytest.approx(layer_plateaus(ref)[n]["F"], rel=2e-3)


def test_r_changes_the_closure_layer_as_the_analytic_formula_says():
    for r in (0.4, 0.6, 0.7):                                   # analytic n* = 3, 4, 6
        n_star = closure_layer_analytic(27.5e-6, 0.43, r)
        run = layered_run(27.5e-6, 0.43, r=r, n_layers=n_star + 1, variant=HYBRID, record_every=10)
        assert run["closure"] is not None and run["closure"]["layer"] == n_star
    # the consistent reference closes later, and by more than one layer once r is close to 1
    run = layered_run(27.5e-6, 0.43, r=0.6, n_layers=7, variant=CONSISTENT, record_every=10)
    assert run["closure"] is not None and run["closure"]["layer"] >= closure_layer_analytic(27.5e-6, 0.43, 0.6) + 2


def test_layer_transients_last_a_few_band_transits():
    run = layered_run(27.5e-6, 0.43, n_layers=4, record_every=1)
    for n in (2, 3, 4):
        tr = transition_duration(run, n)
        assert tr is not None and 0.2 < tr["tau_over_band"] < 30
        # ... which is far below the layer transit time
        idx = np.where(run["layer_top"] == n)[0]
        assert tr["tau_s"] < 0.05 * (run["t"][idx[-1]] - run["t"][idx[0]])


# ------------------------------------------------------------------ solver audit
def test_face_averaging_and_grid_do_not_move_the_plateaus_or_the_closure():
    ref = layered_run(27.5e-6, 0.43, n_layers=5, record_every=10)
    for kw in (dict(face_average="harmonic"), dict(n_band=80), dict(dt_over_band=0.25), dict(stretch=1.015)):
        run = layered_run(27.5e-6, 0.43, n_layers=5, record_every=10, **kw)
        for n in (2, 3, 4):
            assert layer_plateaus(run)[n]["F"] == pytest.approx(layer_plateaus(ref)[n]["F"], rel=2e-3)
        assert run["closure"]["layer"] == ref["closure"]["layer"]
        assert abs(run["closure"]["h"] - ref["closure"]["h"]) < 3 * ref["delta_C"][0]


def test_closure_detection_on_the_first_picard_iterate_is_dt_sensitive_and_the_converged_check_is_not():
    """Audit finding: the original solver abandons the step as soon as ANY Picard iterate has a band porosity
    <= eps_stop.  The first iterate uses the lagged sink, so a large step overshoots: at dt = 2 delta_C/v the
    27.5 um run closes on layer 4 instead of 5.  Testing the converged step instead (closure_check =
    "converged", extended copy only) gives layer 5 for dt = 0.5, 2 and 4 delta_C/v."""
    layers = {}
    for dtf in (0.5, 2.0):
        for chk in ("iterate", "converged"):
            run = layered_run(27.5e-6, 0.43, n_layers=5, record_every=10, dt_over_band=dtf, closure_check=chk)
            layers[(dtf, chk)] = run["closure"]["layer"]
    assert layers[(0.5, "iterate")] == layers[(0.5, "converged")] == 5
    assert layers[(2.0, "converged")] == 5 and layers[(2.0, "iterate")] == 4


def test_infeasibility_claim_holds_up_to_eps_new_0_99_on_two_grids_and_is_qualified_beyond():
    """On the 0.24 um crystals (layer 5 of 27.5 um beads) no birth porosity <= 0.99 reaches eps_final = 0.43 in the
    consistent model on either grid; the linearised model reaches it at 0.93.  Beyond 0.99 the dense scan of
    solver_audit.py finds a grid-dependent sliver (eps_new >= 0.994-0.998) where a_v = 5(1-eps)/d_c -> 0 makes the
    band inert; a packing formula has no meaning there, so the claim is stated for eps_new <= 0.99 only."""
    d_sub = 0.241e-6 / P0.dc_over_dsub
    best = {}
    for nb in (40, 80):
        ef = []
        for e in (0.5, 0.7, 0.9, 0.97, 0.99):
            try:
                ef.append(homogeneous_steady(d_sub, e, U_EVAP, P0, n_band=nb, **CONSISTENT)["eps_final"])
            except RuntimeError:
                ef.append(0.0)
        best[nb] = max(ef)
    assert best[40] < 0.43 and best[80] < 0.43
    r = homogeneous_steady(d_sub, 0.93, U_EVAP, P0, n_band=40, **HYBRID)
    assert r["eps_final"] > 0.43
