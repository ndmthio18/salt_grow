"""The quasi-1-D hierarchical supply check that gates growth (ideas/plan/13 v2).

Wijnhorst topology: the pore size shrinks by the structure's ``pore_ratio`` and the layer
thickness doubles per layer, permeability reduced by the factor c_k.  A face can keep growing
while the capillary threshold of the top layer exceeds gravity + the viscous drop over the stack
+ the suction at its footprint.  The stack is the one ``crust.BirthClosure`` gives the voxels,
from the same constants and the same ``crust.SaltStructure``.
"""
from __future__ import annotations

import numpy as np

from ..materials.brine import GAMMA, MU, RHO_L
from ..materials.constants import G
from .crust import A_MIN, ALPHA, C_H, SaltStructure


def _layer(n, a1, d1, c_k, structure: SaltStructure):
    """(a_n, k_n) of layer n: the same expressions as ``BirthClosure.pore_radius`` and
    ``.permeability`` (tests/test_P3_state_kernel.py holds the two together)."""
    r = structure.pore_ratio
    a_n = max(a1 * r ** (n - 1), A_MIN)
    return a_n, c_k * (2.0 * r ** (1 - n)) * a_n ** 3 / (12 * d1)


def supply_need(h, U, d_sub, s_foot, c_k, structure: SaltStructure):
    """(need, Pc_n): cumulative viscous drop over the stack of ``structure`` up to height h
    plus gravity plus footprint suction, and the capillary threshold of the top layer.  The
    structure is an argument without a default: the gate must be asked about the salt that is
    actually being grown."""
    d1 = 0.14 * d_sub; a1 = d1 / 4; h1 = C_H * d_sub
    dp = 0.0; h_cum = 0.0; n = 0; a_n = a1
    while h_cum < h and n < 40:
        n += 1
        a_n, k_n = _layer(n, a1, d1, c_k, structure)
        h_n = min(h1 * ALPHA ** (n - 1), h - h_cum)
        dp += MU * U * h_n / k_n
        h_cum += h_n
    Pc = 2 * GAMMA / a_n
    return s_foot + RHO_L * G * h + dp, Pc


def first_failure_height(U, d_sub, s_foot, c_k, h_top, structure: SaltStructure):
    """Lowest height h in (0, h_top] at which need(h) > Pc(h) for the stack of supply_need, or inf.

    need(h) is continuous and linear within each layer (slope rho g + mu U / k_n) while Pc jumps up
    by 1 / pore_ratio at every layer boundary, so need - Pc is only piecewise monotonic: it can
    fail, recover at the next boundary and fail again.  A single bisection over [0, h_top] may return a later
    crossing; the growth process is stopped by the first one.  This scans the layers from the
    substrate upwards and returns the first crossing exactly (linear solve inside the layer)."""
    d1 = 0.14 * d_sub; a1 = d1 / 4; h1 = C_H * d_sub
    if U <= 0.0:                                              # no viscous demand (condensing faces): only
        U = 0.0                                               # gravity and the footprint suction remain
    need0 = s_foot; h_cum = 0.0
    for n in range(1, 41):
        a_n, k_n = _layer(n, a1, d1, c_k, structure)
        Pc = 2 * GAMMA / a_n
        h_layer = h1 * ALPHA ** (n - 1)
        h_end = min(h_cum + h_layer, h_top)
        slope = RHO_L * G + MU * U / k_n                      # d need / d h inside layer n
        if need0 + RHO_L * G * h_cum > Pc:                    # already failing at the layer bottom
            return h_cum
        h_cross = h_cum + (Pc - need0 - RHO_L * G * h_cum) / slope
        if h_cross <= h_end:
            return h_cross
        need0 += MU * U * (h_end - h_cum) / k_n               # viscous drop accumulated over this layer
        h_cum = h_end
        if h_cum >= h_top:
            return np.inf
    return np.inf
