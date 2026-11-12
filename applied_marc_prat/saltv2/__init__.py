"""saltv2 -- evaporation-wicking salt growth: a voxel gas / salt domain on a substrate pore network.

    materials/   property values: brine, moist air, the salt, constants.  Imports nothing else.
    pnm/         the substrate pore network, aligned with DuMuX: geometry, pressure, salt
                 transport, the opt-in in-pore precipitation.  Imports materials only.
    growth/      the voxel domain: lattice and states, occupancy growth, vapour diffusion, the
                 per-voxel salt properties and sealing, the hierarchical supply check, the
                 Prat (2026) growth law.  Imports materials only.
    coupling/    where the two meet: pore openings on the voxel grid, the mass ledger, and the
                 closed loops (driver3d, driver2d).  Imports the three above.
    cases/       every value that would change with the experiment.  Imported by launchers and
                 tests, never by the four packages above.

``pnm`` and ``growth`` do not import each other.  The names below are the voxel-side API that
tests and scripts have always imported from the package root.
"""
from .growth.grid import GAS, SALT_SEALED, SALT_WET, Grid
from .growth.geometry import GeometryEngine
from .growth.vapor import VaporProperties, solve_vapor
from .growth.crust import BirthClosure, CrustField, effective_vapour_diffusivity
from .growth.regularize import smooth_on_mask
from .growth import morphology
from .coupling.ledger import MassLedger
from .coupling.interface import PoreOpenings

__all__ = ["GAS", "SALT_WET", "SALT_SEALED", "Grid", "GeometryEngine",
           "VaporProperties", "solve_vapor", "BirthClosure", "CrustField",
           "effective_vapour_diffusivity",
           "MassLedger", "smooth_on_mask", "PoreOpenings", "morphology"]
