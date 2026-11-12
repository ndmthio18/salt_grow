"""Interface conservation of ``pnm.network_step`` with a held (crystal-touching) surface pore
(plan/19 S1, priority 1).  Regime: evaporation-wicking, saturated rigid column, no vapour field.

Four books, each closed independently, per step and cumulatively, with the common yardstick
``pnm.step_budget`` at the densities the kernel used for its fluxes (the old ones):

    dissolved salt   S = sum m                 dS = salt_in - X
    solid salt       X = what left the brine of the held pores (>= 0 by the complementarity)
    total liquid     L = sum V rho(w)          dL = liquid_in - q_water - X     (S1 convention)
    water            W = L - S                 dW = (liquid_in - salt_in) - q_water

S1 convention (``crystal_liquid_sink=True``): ``q_mass`` is the WATER a pore loses, X every salt it
loses.  The production default keeps X out of the liquid book: whatever salt share ``q_mass``
does not declare is created as water, X per step when ``q_mass`` is water only (here), and
X - C_sat q_mass dt when the brine the faces draw is in ``q_mass`` (the production driver).
The magnitude is pinned here, not hidden.  Tolerances are relative to the inventories (liquid
L ~ 2e-6 kg, dissolved salt S ~ 5e-7 kg of the 10-pore column): round-off is 1e-16 of them.

A column at C_sat fed by a reservoir at C_sat has no transient at all, so the exact
evaporated-water-to-salt conversion of plan/19 section 3.2, X = C_sat / (1 - C_sat) q dt, must
hold from the first step to round-off; a reservoir below saturation delivers less than that.
"""
import numpy as np
import pytest

from saltv2.cases import stage_t_wicking as st
from saltv2.materials import brine
from saltv2.materials.brine import C_SAT
from saltv2.pnm import network_step, step_budget

DT = 30.0            # [s]
N_STEPS = 60


def march(C_res, w0, sink, n_layers=10, dt=DT, n_steps=N_STEPS):
    layer = st.build_network("C", n_side=1, n_layers=n_layers)
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    q = st.demand(layer, "uniform")                     # pure water lost by the surface pore
    held = np.zeros(layer.N, bool); held[layer.surface] = True
    V = np.asarray(layer.V, float)
    w = np.full(layer.N, w0)
    m = brine.salt_mass(w, V)
    rho_res = float(brine.density(C_res))
    books = dict(dissolved=[], liquid=[], water=[], X=[], X_over_q=[])
    L_ref = float((V * brine.density(C_SAT)).sum()); S_ref = C_SAT * L_ref
    for _ in range(n_steps):
        w_old = w; rho_old = brine.density(w_old)
        S0 = float(m.sum()); L0 = float((V * rho_old).sum())
        step = network_step(layer, m, C_res, q, g_lat, g_res, d_lat, d_res, dt,
                            crystal=held, w_eq=C_SAT, crystal_liquid_sink=sink)
        m = step.m; w = brine.mass_fraction_from_salt_mass(m, V)
        X = float(step.crystallised.sum())
        bud = step_budget(layer, w, rho_old, step.Q_lat, step.Q_res, d_lat, d_res, C_res, rho_res)
        liq_in = float(bud["liquid_res"].sum()) * dt; salt_in = float(bud["salt_res"].sum()) * dt
        S1 = float(m.sum()); L1 = float((V * brine.density(w)).sum()); out = float(q.sum()) * dt
        books["dissolved"].append((S1 - S0) - (salt_in - X))
        books["liquid"].append((L1 - L0) - (liq_in - out - X))
        books["water"].append(((L1 - S1) - (L0 - S0)) - ((liq_in - salt_in) - out))
        books["X"].append(X); books["X_over_q"].append(X / out)
    b = {k: np.array(v) for k, v in books.items()}
    b["tol_L"] = 1e-15 * L_ref; b["tol_S"] = 1e-15 * S_ref
    return b, float(q.sum()) * dt


def test_four_books_close_with_the_liquid_sink_saturated_column():
    b, out = march(C_SAT, C_SAT, sink=True)
    assert np.all(b["X"] > 0.0)
    assert np.max(np.abs(b["dissolved"])) < b["tol_S"]     # salt equation exact
    assert np.max(np.abs(b["liquid"])) < b["tol_L"]        # within STORAGE_RTOL, per step
    assert np.max(np.abs(b["water"])) < b["tol_L"]
    assert abs(b["dissolved"].sum()) < 10 * b["tol_S"] and abs(b["water"].sum()) < 10 * b["tol_L"]
    # exact conversion: no transient, inflow at C_sat, so X = C_sat/(1 - C_sat) q dt every step
    assert np.allclose(b["X_over_q"], C_SAT / (1.0 - C_SAT), rtol=1e-12, atol=0.0)


def test_without_the_sink_water_is_created_at_the_rate_X():
    b, out = march(C_SAT, C_SAT, sink=False)
    assert np.max(np.abs(b["dissolved"])) < b["tol_S"]     # the salt equation is still exact ...
    assert np.allclose(b["water"], b["X"], rtol=1e-9, atol=b["tol_L"])    # ... but water is created
    assert np.allclose(b["liquid"], b["X"], rtol=1e-9, atol=b["tol_L"])
    assert np.allclose(b["X_over_q"], C_SAT, rtol=1e-12)   # only the water-borne salt arrives
    assert b["water"].sum() > 0.2 * N_STEPS * out          # 26 % of the evaporated water here


def test_books_close_in_a_transient_with_a_subsaturated_reservoir():
    """Reservoir at C0 = 0.25 below C_sat, column at C_sat: diffusion drains salt downwards, the
    held pore gets less than the conversion, the complementarity keeps X >= 0, books still close."""
    b, out = march(0.25, C_SAT, sink=True, n_steps=120)
    assert np.all(b["X"] >= 0.0)
    assert np.max(np.abs(b["dissolved"])) < b["tol_S"]
    assert np.max(np.abs(b["water"])) < b["tol_L"] and np.max(np.abs(b["liquid"])) < b["tol_L"]
    assert np.all(b["X_over_q"] < C_SAT / (1.0 - C_SAT))
    assert b["X_over_q"][-1] > 0.5 * C_SAT / (1.0 - C_SAT)   # most of the advected salt still arrives


def test_seed_mass_is_the_first_crystallised_salt_not_an_injection():
    """A structure's first mass is the X of the step in which the pore is first held: nothing is
    added that the pore did not deliver.  The cumulative solid book equals the salt the dissolved
    book lost beyond its inflow."""
    b, out = march(C_SAT, C_SAT, sink=True, n_steps=5)
    assert b["X"][0] > 0.0 and b["X"][0] == pytest.approx(C_SAT / (1.0 - C_SAT) * out, rel=1e-12)
