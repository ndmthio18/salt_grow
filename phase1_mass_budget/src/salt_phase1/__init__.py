"""Independent Phase 1 evaporation and salt-mass accounting tools."""

from .capacity import (
    CapacityDiagnostic,
    diagnose_domain_capacity,
    required_domain_height_m,
)
from .evaporation import (
    EvaporationParameters,
    MassBudgetResult,
    compute_mass_budget,
    saturation_vapor_concentration_kg_m3,
    saturation_vapor_pressure_pa,
)
from .voxel import (
    AllocationProposal,
    AllocationResult,
    MassToVoxelAllocator,
    TransactionResult,
    TransactionalMassToVoxelAllocator,
    allocate_salt_mass,
)

__all__ = [
    "AllocationProposal",
    "AllocationResult",
    "CapacityDiagnostic",
    "EvaporationParameters",
    "MassBudgetResult",
    "MassToVoxelAllocator",
    "TransactionResult",
    "TransactionalMassToVoxelAllocator",
    "allocate_salt_mass",
    "compute_mass_budget",
    "diagnose_domain_capacity",
    "required_domain_height_m",
    "saturation_vapor_concentration_kg_m3",
    "saturation_vapor_pressure_pa",
]
