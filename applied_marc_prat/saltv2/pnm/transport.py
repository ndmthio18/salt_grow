"""Dissolved-salt transport between the pores: backward-Euler upwind advection + diffusion.

``salt_step`` is the transport solve at a given flow; ``network_step`` is what a time loop calls:
flow and salt of one step, consistent with each other, conserving both masses."""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from ..materials import brine
from .flow import solve_pressure


def salt_step(layer, m, C_res, Q_lat, Q_res, d_lat, d_res, dt, crystal=None, w_eq=None):
    """Backward-Euler upwind advection + diffusion of the dissolved salt between pores;
    reservoir at fixed mass fraction ``C_res``.  Unconditionally stable.

    NOT conservative on its own: the storage term carries the density of the OLD state (see
    Formulation), so the salt a pore gains through its own density change, V w d rho, has no
    flux behind it.  A time loop calls ``network_step``, which closes that; this function is the
    linear solve inside it and the stand-alone form of the scheme every record before
    2026-09-21 was made with.

    State.  ``m`` is the dissolved NaCl MASS per pore [kg]; the pore's liquid volume is
    ``layer.V`` and the brine density follows the salt it holds (``materials.brine``, Batzle & Wang
    as in DuMuX), so the mass fraction is ``brine.mass_fraction_from_salt_mass(m, V)``.

    Formulation.  The MASS-averaged reference system, which is what the DuMuX reference case
    turned out to use as well (measured from its own steady profile: the mass balance returns its
    stated D = 1.54e-9 m^2/s per throat to five digits, the molar one would need 1.37e-9):

        V rho(w) d w / dt = sum_throats [ Q rho^up w^up + d rho^face (w_j - w_i) ] - S

    with w the NaCl mass fraction, Q the volumetric liquid flow [m^3/s] and d = D A / L the
    throat's diffusive conductance [m^3/s] (``layer.diffusive``).  The densities are lagged by one
    step (evaluated at w^n), which keeps the system linear; the step is limited elsewhere so that
    a pore changes little per step.  Before 2026-09-17 the density was the constant RHO_L = 1200
    everywhere, which is 8 % too high for a 15 wt% brine.

    ``crystal`` (bool per pore) marks pores whose brine touches crystal and ``w_eq`` is the
    equilibrium mass fraction there.  Such a pore cannot exceed ``w_eq``: it is held AT ``w_eq``
    and what the throats carry into it beyond that crystallises,

        X_i = dt * (net salt inflow at w_i = w_eq) - V rho (w_eq - w_i^n)   >= 0      [kg]

    solved implicitly with the rest of the network, so the neighbours see ``w_eq`` and not an
    overshoot that is clipped afterwards.  Where X_i would be negative -- salt diffuses away
    faster than it arrives, the crystal would have to dissolve, and nothing here can -- the pore
    is left free instead and falls below ``w_eq`` (a complementarity: X_i >= 0, w_i <= w_eq, one
    of them an equality).  With ``crystal`` given the return value is ``(m_new, X)``.
    """
    w_new, X = _salt_solve(layer, m, C_res, Q_lat, Q_res, d_lat, d_res, dt, crystal, w_eq, None,
                           brine.T_REF)
    if crystal is not None:
        return brine.salt_mass(w_new, layer.V), X
    return brine.salt_mass(w_new, layer.V)


def _salt_solve(layer, m, C_res, Q_lat, Q_res, d_lat, d_res, dt, crystal, w_eq, rho_new, T):
    """(w_new, X): the linear system of ``salt_step``, with the storage term

        V (rho_new w_new - rho w) / dt        instead of        V rho (w_new - w) / dt

    when ``rho_new`` [kg/m^3 per pore], a guess of the brine density at the END of the step, is
    given.  The rows stay normalised by the old liquid mass V rho, so ``rho_new`` enters as the
    factor rho_new / rho on the diagonal and ``rho_new = None`` is the old scheme bit for bit.
    The throat fluxes use the old densities either way, as ``solve_pressure`` does."""
    N = layer.N
    V = np.asarray(layer.V, float)
    w = brine.mass_fraction_from_salt_mass(m, V, T=T)        # mass fraction now
    rho = brine.density(w, T=T)                              # [kg/m^3], lagged
    rho_res = float(brine.density(C_res, T=T))
    store = V * rho                                          # liquid mass per pore [kg]
    diag = np.zeros(N)
    rhs = w.copy()
    # reservoir throat: upwind advection plus diffusion against the fixed reservoir
    inflow = np.maximum(Q_res, 0.0); outflow = np.maximum(-Q_res, 0.0)
    rho_face_res = 0.5 * (rho + rho_res)
    rhs += dt * (inflow * rho_res + d_res * rho_face_res) * C_res / store
    diag += dt * (outflow * rho + d_res * rho_face_res) / store
    # lateral throats, Q > 0: a -> b.  Loss from a: (qa rho_a + d rho_face) w_a; gain to a from b:
    # (qb rho_b + d rho_face) w_b.  Vectorised 2026-09-21 (a column of 1273 pores spent most of
    # its time in the throat-by-throat loop); ``np.add.at`` is sequential and the two ends are
    # interleaved a0, b0, a1, b1, ..., so the diagonal accumulates in the order that loop used
    # and the matrix is the same to the last bit (as in ``flow.solve_pressure``).
    a_idx = np.array([t[0] for t in layer.lat], int)
    b_idx = np.array([t[1] for t in layer.lat], int)
    Q_lat = np.asarray(Q_lat, float); d_lat = np.asarray(d_lat, float)
    qa = np.maximum(Q_lat, 0.0); qb = np.maximum(-Q_lat, 0.0)
    rho_face = 0.5 * (rho[a_idx] + rho[b_idx])
    out_a = qa * rho[a_idx] + d_lat * rho_face           # what a loses to b, per unit w_a
    out_b = qb * rho[b_idx] + d_lat * rho_face
    ab = np.empty(2 * a_idx.size, int); ab[0::2] = a_idx; ab[1::2] = b_idx
    ba = np.empty(2 * a_idx.size, int); ba[0::2] = b_idx; ba[1::2] = a_idx
    loss = np.empty(2 * a_idx.size); loss[0::2] = dt * out_a / store[a_idx]
    loss[1::2] = dt * out_b / store[b_idx]
    gain = np.empty(2 * a_idx.size); gain[0::2] = -dt * out_b / store[a_idx]
    gain[1::2] = -dt * out_a / store[b_idx]
    np.add.at(diag, ab, loss)
    unit = 1.0 if rho_new is None else np.asarray(rho_new, float) / rho
    self_idx = np.arange(N)
    A = coo_matrix((np.concatenate([gain, unit + diag]),
                    (np.concatenate([ab, self_idx]), np.concatenate([ba, self_idx]))),
                   shape=(N, N)).tocsr()
    X = np.zeros(N)
    if crystal is None:
        w_new = spsolve(A, rhs)
    else:
        held = np.asarray(crystal, bool).copy()
        row_of = np.repeat(np.arange(N), np.diff(A.indptr))      # row of every stored entry
        while True:
            # a held row becomes  w_i = w_eq:  its entries go, a unit diagonal stays
            A_held = A.copy()
            in_held = held[row_of]
            A_held.data[in_held] = np.where(A_held.indices[in_held] == row_of[in_held], 1.0, 0.0)
            A_held.eliminate_zeros()
            rhs_held = np.where(held, w_eq, rhs)
            w_new = spsolve(A_held, rhs_held)
            X = np.where(held, store * (rhs - A @ w_new), 0.0)       # what the held rows leave over
            if not np.any(X < 0.0):
                break
            held &= ~(X < 0.0)
    if np.any(w_new < 0.0):
        # the structure asked for more salt than the pore holds: report it, do
        # not silently clip a negative mass into the next step
        raise RuntimeError(
            f"salt sink emptied {int((w_new < 0).sum())} pore(s) in one step "
            f"(most negative mass fraction {w_new.min():.3e}); reduce dt")
    return w_new, X


class NetworkStep(NamedTuple):
    """What ``network_step`` returns."""
    m: np.ndarray            # dissolved NaCl per pore after the step [kg]
    crystallised: np.ndarray  # X: what left the brine of the held pores [kg per pore]
    p: np.ndarray            # pressures of the step [Pa]
    Q_lat: np.ndarray        # volumetric throat flows of the step [m^3/s]
    Q_res: np.ndarray
    storage: np.ndarray      # liquid mass the pores gained through their density [kg/s per pore]
    salt_in: float           # salt that crossed the reservoir throats, net [kg]


# ``network_step``: how well the density the storage terms were built with has to agree with the
# density of the state they produce [N].  The salt balance is exact whatever this is (see there);
# it bounds the liquid-mass residual, per step and relative to the pore's liquid.
STORAGE_RTOL = 1e-12
STORAGE_MAX_ITER = 25
# With ``crystal_liquid_sink`` the salt sink X of a held pore is iterated as well; each pass
# contracts its error by about the inflow mass fraction (~0.26), so 1e-12 needs ~25 passes [N].
SINK_MAX_ITER = 100


def network_step(layer, m, C_res, q_mass, g_lat, g_res, d_lat, d_res, dt, crystal=None,
                 w_eq=None, flow=None, T=brine.T_REF, crystal_liquid_sink=False):
    """One backward-Euler step of the liquid AND the salt of the network, each conserved:

        V (rho' - rho) / dt      = sum_throats rho^up Q            - q_mass        (liquid)
        V (rho' w' - rho w) / dt = sum_throats [rho^up Q w^up + d rho^face dw] - X / dt   (salt)

    with rho' = rho(w') the density of the NEW state.  Until 2026-09-21 the first left-hand side
    was zero (``solve_pressure`` has no storage term) and the second was V rho (w' - w) / dt
    (``salt_step``): a pore whose brine gets denser gained V w d rho of salt and V d rho of
    liquid that no throat had carried -- bounded, -0.21 % of the reservoir salt on the DuMuX
    column, and absent in any steady state, which is why no steady check saw it.

    Method.  rho' is iterated: with a guess rho'_k the storage V (rho'_k - rho) / dt joins
    ``q_mass`` in ``solve_pressure`` (whose signature and assembly are untouched), the salt rows
    take rho'_k on their storage term, and rho'_(k+1) = rho(w'); the throat densities stay those
    of the old state in both solves, so what one pore loses the next one gains.  The state
    returned is  m' = V rho'_k w'  of the LAST pair, which satisfies the discrete salt balance
    to round-off for any tolerance; ``STORAGE_RTOL`` only bounds how far rho'_k is from
    rho(w'), i.e. the liquid residual.  A steady state leaves the loop on its first pass with
    zero storage.  It raises rather than returning an unconverged pair.

    ``q_mass`` [kg/s per pore, > 0 = leaving] as for ``solve_pressure``; ``crystal`` / ``w_eq``
    as for ``salt_step``; ``flow`` = (p, Q_lat, Q_res) already solved for this ``q_mass`` and this
    ``m`` saves the first pressure solve (storage is zero on the first pass); ``T`` [K] is the
    temperature every brine density of the step is taken at.

    ``crystal_liquid_sink`` (OPT-IN, 2026-10-26, plan/19 S1): the salt X that leaves the brine of
    a held pore is a mass sink of its LIQUID balance as well (the pore stays liquid-filled, so the
    network has to deliver X of brine behind it).  With it, ``q_mass`` is the WATER a pore loses
    and X every salt it loses, and the four books (dissolved salt, solid salt, total liquid,
    water) close to round-off (tests/test_interface_conservation.py); X is then iterated with
    rho' (it raises the inflow, which carries salt).  The default keeps the production
    convention, in which the brine the faces draw -- water AND salt -- is already in ``q_mass``
    and the interface buffer of the driver absorbs X - q_salt dt: there the water book is off by
    X - C_sat q_mass dt per step in a transient (zero in a steady state).

    Returns a ``NetworkStep``: the new salt masses, what crystallised in the held pores (zeros
    without ``crystal``), the flow that belongs to the step, the liquid mass the pores gained
    through their density, and the net salt that crossed the reservoir throats -- advection
    upwind plus diffusion at the mean density, with the w' the step solved for, i.e. exactly the
    reservoir terms of the rows, so that  sum(m_new - m) = salt_in - sum(X)  to round-off.
    """
    V = np.asarray(layer.V, float)
    rho = brine.density(brine.mass_fraction_from_salt_mass(m, V, T=T), T=T)
    rho_res = float(brine.density(C_res, T=T))
    q_mass = np.asarray(q_mass, float)
    rho_new = rho
    X = np.zeros(layer.N)
    n_iter = SINK_MAX_ITER if crystal_liquid_sink else STORAGE_MAX_ITER
    for k in range(n_iter):
        storage = V * (rho_new - rho) / dt
        demand = q_mass + storage
        if crystal_liquid_sink:
            demand = demand + X / dt                 # the previous iterate's X, a liquid sink
        if k == 0 and flow is not None:
            p, Q_lat, Q_res = flow
        else:
            p, Q_lat, Q_res = solve_pressure(layer, demand, g_lat, g_res, rho, rho_res)
        w_new, X_new = _salt_solve(layer, m, C_res, Q_lat, Q_res, d_lat, d_res, dt, crystal, w_eq,
                                   rho_new, T)
        rho_next = brine.density(w_new, T=T)
        settled = np.max(np.abs(rho_next - rho_new) / rho_next) <= STORAGE_RTOL
        if crystal_liquid_sink:
            settled = settled and (np.max(np.abs(X_new - X))
                                   <= STORAGE_RTOL * max(float(np.max(np.abs(X_new))), 1e-300))
        X = X_new
        if settled:
            rho_face = 0.5 * (rho + rho_res)
            salt_in = dt * float(np.sum(
                (np.maximum(Q_res, 0.0) * rho_res + d_res * rho_face) * C_res
                - (np.maximum(-Q_res, 0.0) * rho + d_res * rho_face) * w_new))
            return NetworkStep(rho_new * V * w_new, X, p, Q_lat, Q_res, storage, salt_in)
        rho_new = rho_next
    raise RuntimeError(
        f"network_step: the new-state density did not settle in {n_iter} passes "
        f"(max relative change {float(np.max(np.abs(rho_next - rho_new) / rho_next)):.3e} > "
        f"{STORAGE_RTOL:.0e}); reduce dt")
