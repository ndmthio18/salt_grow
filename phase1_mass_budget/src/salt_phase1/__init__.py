"""Independent Phase 1 evaporation and salt-mass accounting tools."""

from .evaporation import (
    EvaporationParameters,
    MassBudgetResult,
    compute_mass_budget,
    saturation_vapor_concentration_kg_m3,
    saturation_vapor_pressure_pa,
)
from .voxel import AllocationResult, MassToVoxelAllocator, allocate_salt_mass

__all__ = [
    "AllocationResult",
    "EvaporationParameters",
    "MassBudgetResult",
    "MassToVoxelAllocator",
    "allocate_salt_mass",
    "compute_mass_budget",
    "saturation_vapor_concentration_kg_m3",
    "saturation_vapor_pressure_pa",
]
