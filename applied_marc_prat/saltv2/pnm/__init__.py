"""The substrate pore network, aligned with the DuMuX pnm2pncmin2 reference case.

    dgf            parse_dgf_2d, write_column_dgf
    network        PoreLayer (generated lattice, or PoreLayer.from_dgf), throat conductances
    flow           solve_pressure: mass-conserving Kirchhoff solve with upwind densities
    transport      network_step: one step of liquid and salt, both conserved (pore storage
                   included); salt_step: the transport solve inside it
    conservative   the Stage-T OPT-IN step (conservative_step, kept only for cross-validation
                   against network_step) and its utilities: throat_kind / drop_throats (paired
                   networks), step_residuals / step_budget / net_inflow, hydrostatic_head
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
from .network import DGF_DUMUX_1D, DGF_DUMUX_5X5_VARTHROAT, PoreLayer
from .precipitation import PoreClogged, PorePrecipitation
from .transport import network_step, salt_step

__all__ = ["PoreLayer", "DGF_DUMUX_1D", "DGF_DUMUX_5X5_VARTHROAT", "parse_dgf_2d",
           "write_column_dgf", "solve_pressure",
           "salt_step", "network_step", "boundary_layer_evaporation", "PorePrecipitation", "PoreClogged",
           "conservative_step", "step_budget", "step_residuals", "net_inflow", "throat_kind",
           "drop_throats", "hydrostatic_head", "SURFACE_LATERAL", "SUBSURFACE_LATERAL", "VERTICAL"]
