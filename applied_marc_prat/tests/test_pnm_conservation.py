"""``pnm.network_step`` conserves the liquid AND the salt of the network in a transient.

Regime: evaporation-wicking -- a saturated column fed by a Dirichlet reservoir, its top pore
losing liquid.  Both balances are checked from the OUTSIDE: what crosses the reservoir throat is
recomputed here from the flows and states the step returns, and the inventories are read from
the states through ``materials.brine``; nothing is taken from the step's own bookkeeping except
where that bookkeeping is the thing under test (``storage``).

The case is deliberately violent (a 10 wt% column driven to saturation in minutes, with a time
step that jumps by up to 15x): the terms at stake are V w d rho and V d rho, which vanish in a
steady state -- the reason no steady check ever saw that they were missing.
"""
import numpy as np
import pytest

from saltv2.materials import brine
from saltv2.pnm import DGF_DUMUX_1D, PoreLayer, network_step, salt_step, solve_pressure

C_RES = 0.10
SINK = 4e-10                 # liquid mass the top pore loses [kg/s]


def column():
    layer = PoreLayer.from_dgf(DGF_DUMUX_1D)
    return layer, layer.conductances(), layer.diffusive(), int(layer.surface[0])


def reservoir_salt(layer, rho, rho_res, Q_res, d_res, w_new, dt):
    """Salt that crossed the reservoir throats in one step [kg], by salt_step's own rule:
    upwind advection plus diffusion at the mean density, implicit in the pore's mass fraction."""
    rho_face = 0.5 * (rho + rho_res)
    gain = (np.maximum(Q_res, 0.0) * rho_res + d_res * rho_face) * C_RES
    loss = (np.maximum(-Q_res, 0.0) * rho + d_res * rho_face) * w_new
    return dt * float((gain - loss).sum())


def march(scheme, n_steps=300, crystal_at=None):
    layer, (g_lat, g_res), (d_lat, d_res), top = column()
    V = layer.V
    rho_res = float(brine.density(C_RES))
    m = brine.salt_mass(C_RES, V)
    book = dict(salt_in=0.0, crystallised=0.0, liquid_in=0.0, liquid_out=0.0, stored=0.0,
                m0=float(m.sum()), L0=float((V * brine.density(C_RES)).sum()))
    rng = np.random.default_rng(1)
    crystal = None
    for _ in range(n_steps):
        dt = float(rng.uniform(0.2, 3.0))
        q = np.zeros(layer.N); q[top] = SINK
        rho = brine.density(brine.mass_fraction_from_salt_mass(m, V))
        if crystal_at is not None:
            crystal = np.zeros(layer.N, bool); crystal[top] = True
        if scheme == "network_step":
            step = network_step(layer, m, C_RES, q, g_lat, g_res, d_lat, d_res, dt,
                                crystal=crystal, w_eq=crystal_at)
            m_new, X, _, Q_lat, Q_res, storage = step[:6]
            w_new = m_new / (V * (rho + storage * dt / V))          # the w' the step solved for
            book["stored"] += dt * float(storage.sum())
        else:                                                       # the scheme until 2026-09-21
            _, Q_lat, Q_res = solve_pressure(layer, q, g_lat, g_res, rho, rho_res)
            m_new = salt_step(layer, m, C_RES, Q_lat, Q_res, d_lat, d_res, dt)
            X = np.zeros(layer.N)
            w_new = brine.mass_fraction_from_salt_mass(m_new, V)
        book["salt_in"] += reservoir_salt(layer, rho, rho_res, Q_res, d_res, w_new, dt)
        if scheme == "network_step":                 # the step's own metering is the same number
            assert step.salt_in == pytest.approx(
                reservoir_salt(layer, rho, rho_res, Q_res, d_res, w_new, dt), rel=1e-9, abs=1e-30)
        book["crystallised"] += float(X.sum())
        book["liquid_in"] += dt * float((np.where(Q_res > 0.0, rho_res, rho) * Q_res).sum())
        book["liquid_out"] += dt * SINK
        m = m_new
    w_end = brine.mass_fraction_from_salt_mass(m, V)
    book.update(d_salt=float(m.sum()) - book["m0"],
                d_liquid=float((V * brine.density(w_end)).sum()) - book["L0"], w_end=w_end)
    return book


def test_salt_and_liquid_are_conserved_through_a_strong_transient():
    b = march("network_step")
    assert b["w_end"].max() > 2.5 * C_RES                    # the transient is a real one
    assert b["d_salt"] == pytest.approx(b["salt_in"], rel=1e-12)
    # liquid: reservoir in - sink = what the pores hold more, read from the states ...
    assert b["d_liquid"] == pytest.approx(b["liquid_in"] - b["liquid_out"], rel=1e-9)
    # ... and that is the storage the step booked
    assert b["stored"] == pytest.approx(b["d_liquid"], rel=1e-9)
    assert b["d_liquid"] > 1e-3 * b["liquid_in"]              # storage is not a rounding matter here


def test_the_scheme_it_replaces_did_not_conserve_either():
    """The teeth of the test above: the same march with ``solve_pressure`` + ``salt_step`` gains
    salt and liquid that nothing delivered."""
    b = march("old")
    assert abs(b["d_salt"] - b["salt_in"]) > 1e-2 * b["salt_in"]
    assert b["liquid_in"] - b["liquid_out"] == pytest.approx(0.0, abs=1e-6 * b["liquid_in"])
    assert b["d_liquid"] > 1e-2 * b["liquid_in"]


def test_a_held_pore_crystallises_exactly_what_the_balance_leaves_over():
    """With the top pore held at w_eq the salt that crystallises is the row residual of the
    CONSERVATIVE operator: in = d(pores) + crystallised, to round-off."""
    w_eq = 0.12
    b = march("network_step", crystal_at=w_eq)
    assert b["crystallised"] > 0.1 * b["salt_in"]
    assert b["w_end"].max() == pytest.approx(w_eq, rel=1e-12)
    assert b["d_salt"] + b["crystallised"] == pytest.approx(b["salt_in"], rel=1e-12)
    assert b["stored"] == pytest.approx(b["d_liquid"], rel=1e-9)


def test_a_steady_state_stores_nothing_and_is_the_old_steady_state():
    layer, (g_lat, g_res), (d_lat, d_res), top = column()
    V = layer.V
    q = np.zeros(layer.N); q[top] = 2e-11
    m = brine.salt_mass(C_RES, V)
    for _ in range(4000):
        m, _, _, Q_lat, Q_res, storage = network_step(layer, m, C_RES, q, g_lat, g_res, d_lat,
                                                      d_res, 50.0)[:6]
    assert np.max(np.abs(storage)) < 1e-12 * q[top]
    # one step of the old scheme from here does not move it: the two share their steady states
    rho = brine.density(brine.mass_fraction_from_salt_mass(m, V))
    _, Q_lat, Q_res = solve_pressure(layer, q, g_lat, g_res, rho, float(brine.density(C_RES)))
    m_old = salt_step(layer, m, C_RES, Q_lat, Q_res, d_lat, d_res, 50.0)
    assert np.allclose(m_old, m, rtol=1e-11, atol=0.0)
