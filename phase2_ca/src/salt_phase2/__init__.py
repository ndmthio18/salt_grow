"""Independent CPU Phase 2A salt-growth cellular automaton."""

from .ca import GrowthStepResult, MorphologyDiagnostics, SaltGrowthCA
from .state import Coordinate, Material, PhysicalState, initialize_physical_state

__all__ = [
    "Coordinate",
    "GrowthStepResult",
    "Material",
    "MorphologyDiagnostics",
    "PhysicalState",
    "SaltGrowthCA",
    "initialize_physical_state",
]
