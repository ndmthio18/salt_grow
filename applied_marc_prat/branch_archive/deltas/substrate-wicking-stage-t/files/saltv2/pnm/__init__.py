"""The substrate pore network, aligned with the DuMuX pnm2pncmin2 reference case.

    dgf            parse_dgf_2d, write_column_dgf
    network        PoreLayer (generated lattice, or PoreLayer.from_dgf), throat conductances
    flow           solve_pressure: mass-conserving Kirchhoff solve with upwind densities
    transport      salt_step: backward-Euler upwind advection + diffusion of the pore salt
    conservative   conservative_step: OPT-IN storage-resolving liquid + salt step (Stage-T),
                   throat_kind / drop_throats (paired networks), step_residuals, hydrostatic_head
    precipitation  PorePrecipitation: the OPT-IN in-pore precipitation closure (off by default)
    boundary       boundary_layer_evaporation: the reference case's top boundary condition

Imports ``saltv2.materials`` only; it knows nothing about the voxel domain.
"""
from .boundary import boundary_layer_evaporation
from .conservative import (SUBSURFACE_LATERAL, SURFACE_LATERAL, VERTICAL, conservative_step,
                           drop_throats, hydrostatic_head, net_inflow, step_budget, step_residuals,
                           throat_kind)
from .dgf import parse_dgf_2d, write_column_dgf
from .flow import solve_pressure
from .network import DGF_DUMUX_1D, PoreLayer
from .precipitation import PoreClogged, PorePrecipitation
from .transport import salt_step

__all__ = ["PoreLayer", "DGF_DUMUX_1D", "parse_dgf_2d", "write_column_dgf", "solve_pressure",
           "salt_step", "boundary_layer_evaporation", "PorePrecipitation", "PoreClogged",
           "conservative_step", "step_budget", "step_residuals", "net_inflow", "throat_kind",
           "drop_throats", "hydrostatic_head", "SURFACE_LATERAL", "SUBSURFACE_LATERAL",
           "VERTICAL"]
