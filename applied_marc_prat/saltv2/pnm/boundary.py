"""Zero-dimensional vapour closures for a pore opening (reference cases only)."""
from __future__ import annotations

import numpy as np

from ..materials.air import D_V, RHO_G, W_SAT
from ..materials.brine import a_w


def boundary_layer_evaporation(C, area, delta, rh):
    """Water leaving a pore opening through a stagnant boundary layer [kg/s].

    This is the DuMuX reference case's top boundary condition (``problem.hh::diffWaterFlux``,
    mass form), which is a zero-dimensional vapour closure:

        m_dot = rho_g D_v A (X_pore - X_ref) / delta,      X = a_w(C) W_sat and rh W_sat

    Every factor already exists in ``saltv2.materials``, so nothing new is introduced -- the point
    of having it as a function is that a run can be driven by the pore's OWN concentration instead
    of an imposed flux, which closes the loop evaporation -> pressure -> transport -> a_w.

    Two differences from DuMuX remain:
    they average rho_g and D_v between the pore's gas state and the free-flow reference, while
    these are constants here (about 1-2 %); and the constants are 20 C values, so the reference
    case must be run with ``SpatialParams.Temperature = 293.15`` (its DEFAULT is 283.15 K, where
    p_sat is a factor 1.9 smaller -- the trap this function was written for).

    The sign is not clamped: a pore fresher than the far field would condense, exactly as their
    Neumann condition does.  ``area`` [m^2] is their ``poreCrossSectionalArea`` = (2 r)^2 for a
    cubic pore, i.e. what ``saltv2.coupling.interface.PoreOpenings`` calls A_open.

    NOTE: this is NOT the production 3-D closure.  There the vapour field is solved in the voxel
    domain, so the openings screen each other and the crust conducts; a boundary layer per
    opening cannot represent either.  It is the reference-case closure, for reference cases,
    which is why ``delta`` and ``rh`` have no defaults: they are the case's, not the kernel's
    (the DuMuX 1-D values are in ``saltv2.cases.dumux1d``).
    """
    return RHO_G * D_V * np.asarray(area, float) * W_SAT * (a_w(C) - rh) / delta
