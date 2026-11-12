"""Liquid pressure of the network: Kirchhoff on the liquid MASS balance, reservoir at p = 0."""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

# ``solve_pressure``: how still the pressure field has to be before a run of upwind sweeps that
# never reproduces its own density pattern is accepted as converged anyway [N].  Set at double
# precision, not at an accuracy the physics would care about -- the case it exists for sits at
# 1.7e-14 (see the else branch there); anything a flow reversal actually moves is far larger.
UPWIND_FIXED_POINT_RTOL = 1e-10
TINY = 1e-300                # guard for a divide by the norm of an all-zero field


def solve_pressure(layer, q_mass, g_lat, g_res, rho, rho_res, max_sweeps=20):
    """Kirchhoff on the liquid MASS balance, reservoir at p = 0:

        sum_j rho^up_ij g_ij (p_j - p_i) + rho^up_res,i g_res,i (0 - p_i) = q_mass,i

    ``q_mass`` [kg/s per pore, > 0 = liquid leaving] is what the pore loses: evaporated water,
    brine drawn into the growing structure, liquid displaced by in-pore precipitate.  ``rho``
    [kg/m^3 per pore] and ``rho_res`` are the brine densities; the throat fluxes returned are
    VOLUMETRIC (Q = g dp, what ``salt_step`` and the dt limiters want), only the balance is
    mass-weighted.

    Why not the volumetric balance.  Until 2026-09-17 this solved ``sum_j g_ij (p_j - p_i) =
    q_vol,i``, i.e. it conserved liquid VOLUME.  The brine is not incompressible in composition:
    rho runs from 998 to 1233 kg/m^3 over the concentrations this project uses, so a chain of
    throats carrying one steady mass flux must carry a volume flux that CHANGES along it.  DuMuX
    gets this right -- in its 1-D reference case rho^up Q is constant along the column to six
    digits while the volume flux varies by 3 % -- and the approximation was worth 0.61 % of the
    closed-loop deviation against it, the single largest term.

    The upwind density makes the system non-linear through the flow directions, so the sign
    pattern is iterated to a fixed point (usually one extra sweep: the flow is reservoir -> pore
    -> structure and does not reverse).  It raises rather than returning a half-converged field:
    a pattern that never repeats is accepted only when the pressure field itself has stopped
    moving, to ``UPWIND_FIXED_POINT_RTOL`` (see the ``else`` branch of the sweep loop).
    A uniform ``rho`` reduces it exactly to the old volumetric solve with q_vol = q_mass / rho.
    """
    N = layer.N
    rho = np.broadcast_to(np.asarray(rho, float), (N,))
    rho_res = float(rho_res)
    a_idx = np.array([t[0] for t in layer.lat], int)
    b_idx = np.array([t[1] for t in layer.lat], int)
    g_lat = np.asarray(g_lat, float)
    def upwind(p):
        """Upwind density of every throat at the pressures ``p`` (at p = 0 the lateral throats
        take the mean and the reservoir supplies: the volumetric guess the first sweep starts
        from)."""
        r_lat = np.where(p[a_idx] > p[b_idx], rho[a_idx],
                         np.where(p[b_idx] > p[a_idx], rho[b_idx],
                                  0.5 * (rho[a_idx] + rho[b_idx])))
        return r_lat, np.where(p < 0.0, rho_res, rho)

    # Assembly indices.  ``ab`` interleaves the two ends throat by throat -- a0, b0, a1, b1, ... --
    # so that ``np.add.at``, which is sequential, accumulates the diagonal in the order the
    # throat-by-throat loop this replaced used.  Summation order is not associative in floating
    # point, and the run records are compared bit for bit.
    ab = np.empty(2 * a_idx.size, int); ab[0::2] = a_idx; ab[1::2] = b_idx
    ba = np.empty(2 * a_idx.size, int); ba[0::2] = b_idx; ba[1::2] = a_idx
    self_idx = np.arange(N)
    p = np.zeros(N)
    r_lat, r_res = upwind(p)
    p_prev, moved = None, np.inf
    for _ in range(max_sweeps):
        gr = g_lat * r_lat
        off = np.repeat(gr, 2)
        diag = g_res * r_res
        np.add.at(diag, ab, off)
        A = coo_matrix((np.concatenate([-off, diag]),
                        (np.concatenate([ab, self_idx]), np.concatenate([ba, self_idx]))),
                       shape=(N, N)).tocsr()
        p = spsolve(A, -np.asarray(q_mass, float))     # suction-negative pressures
        # Converged when the densities the NEXT sweep would use are the ones this sweep used:
        # the matrix would be unchanged, so p is already the fixed point.  The test is on the
        # DENSITIES, not on the flow directions -- between two pores holding the same brine the
        # direction can flip on rounding for ever while nothing about the system changes.
        nxt = upwind(p)
        if np.array_equal(nxt[0], r_lat) and np.array_equal(nxt[1], r_res):
            break
        moved = (np.inf if p_prev is None
                 else float(np.max(np.abs(p - p_prev)) / max(float(np.max(np.abs(p))), TINY)))
        p_prev = p
        # The same acceptance as the else branch below, taken as soon as it is true rather than
        # after all ``max_sweeps``: on a network whose density pattern never repeats, every solve
        # would otherwise run the full 20 sweeps for nothing (the resolved column of 2026-09-18 did
        # not finish the 24 h test case in 90 minutes; the single-layer version of it takes 19 s).
        # It cannot pre-empt the test above where that test works: measured over the 2209 solves of
        # that case, every one exits on sweep 1 or 2, and at the exit sweep |dp| / |p| is 2.6e-2
        # median and 0.195 at worst -- nine orders above the tolerance.
        if moved <= UPWIND_FIXED_POINT_RTOL:
            break
        r_lat, r_res = nxt
    else:
        # The density test above can fail for ever on a network large enough to hold two pores at
        # the same concentration to rounding: their Delta rho is then ~1e-14 relative, the flow
        # direction between them flips on the last bit, and the pick oscillates although the
        # SOLUTION has stopped moving.  Measured on the 1006-layer column of 2026-09-18: from the
        # third sweep on, |dp| / |p| stays at 1.7e-14 while 110-150 of the 12068 picks keep
        # flipping.  So before refusing, ask the question that matters -- has the iterate stopped?
        # A genuinely unconverged fixed point (a real flow reversal hunting between two states)
        # moves p by a physical amount and is still refused.  Runs that reach the test above are
        # untouched: this branch is only ever entered where the old code raised.
        if not moved <= UPWIND_FIXED_POINT_RTOL:
            raise RuntimeError(
                f"network upwind densities did not settle in {max_sweeps} sweeps and the pressure "
                f"field is still moving: max|dp| / max|p| = {moved:.3e} > {UPWIND_FIXED_POINT_RTOL:.0e}")
    Q_lat = g_lat * (p[a_idx] - p[b_idx])              # a -> b, volumetric [m^3/s]
    Q_res = g_res * (0.0 - p)                          # reservoir -> pore
    return p, Q_lat, Q_res
