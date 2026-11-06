"""Conservative conversion from precipitated NaCl mass to CA solid voxels."""

from __future__ import annotations

from dataclasses import dataclass
import math


def _positive_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


def _nonnegative_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")


@dataclass(frozen=True)
class AllocationResult:
    """Mass allocation for one call, excluding all pre-existing seed voxels."""

    incoming_salt_mass_kg: float
    previous_residual_mass_kg: float
    available_mass_kg: float
    voxel_mass_kg: float
    new_solid_voxels: int
    allocated_mass_kg: float
    residual_mass_kg: float


def allocate_salt_mass(
    precipitated_salt_mass_kg: float,
    residual_mass_kg: float,
    solid_density_kg_m3: float,
    voxel_size_m: float,
) -> AllocationResult:
    """Allocate whole coarse-grained NaCl voxels and retain the mass remainder.

    Initial crystal seeds are deliberately absent from the inputs and outputs;
    this function accounts only for newly precipitated salt budget.
    """

    _nonnegative_finite("precipitated_salt_mass_kg", precipitated_salt_mass_kg)
    _nonnegative_finite("residual_mass_kg", residual_mass_kg)
    _positive_finite("solid_density_kg_m3", solid_density_kg_m3)
    _positive_finite("voxel_size_m", voxel_size_m)

    voxel_mass = solid_density_kg_m3 * voxel_size_m**3
    available = precipitated_salt_mass_kg + residual_mass_kg
    quotient, remainder = divmod(available, voxel_mass)
    new_voxels = int(quotient)
    allocated = new_voxels * voxel_mass

    # Preserve the algebraic identity using the same operands as allocation.
    # divmod already guarantees 0 <= remainder < voxel_mass for valid inputs.
    remainder = available - allocated
    if remainder < 0.0 and math.isclose(remainder, 0.0, abs_tol=math.ulp(available)):
        remainder = 0.0
    if new_voxels < 0 or allocated < 0.0 or remainder < 0.0:
        raise ArithmeticError("mass allocation produced an invalid negative result")

    return AllocationResult(
        incoming_salt_mass_kg=precipitated_salt_mass_kg,
        previous_residual_mass_kg=residual_mass_kg,
        available_mass_kg=available,
        voxel_mass_kg=voxel_mass,
        new_solid_voxels=new_voxels,
        allocated_mass_kg=allocated,
        residual_mass_kg=remainder,
    )


class MassToVoxelAllocator:
    """Stateful conservative ledger for repeated physical time increments.

    ``initial_seed_voxel_count`` is retained only as separate morphology
    metadata.  It is never added to cumulative newly precipitated salt.
    """

    def __init__(
        self,
        *,
        solid_density_kg_m3: float,
        voxel_size_m: float,
        initial_seed_voxel_count: int = 0,
    ) -> None:
        _positive_finite("solid_density_kg_m3", solid_density_kg_m3)
        _positive_finite("voxel_size_m", voxel_size_m)
        if isinstance(initial_seed_voxel_count, bool) or initial_seed_voxel_count < 0:
            raise ValueError("initial_seed_voxel_count must be a non-negative integer")
        if int(initial_seed_voxel_count) != initial_seed_voxel_count:
            raise ValueError("initial_seed_voxel_count must be an integer")
        self.solid_density_kg_m3 = solid_density_kg_m3
        self.voxel_size_m = voxel_size_m
        self.initial_seed_voxel_count = int(initial_seed_voxel_count)
        self.residual_mass_kg = 0.0
        self.cumulative_new_voxels = 0
        self.cumulative_budget_mass_kg = 0.0

    @property
    def voxel_mass_kg(self) -> float:
        return self.solid_density_kg_m3 * self.voxel_size_m**3

    @property
    def cumulative_allocated_mass_kg(self) -> float:
        return self.cumulative_new_voxels * self.voxel_mass_kg

    @property
    def total_morphology_voxel_count(self) -> int:
        return self.initial_seed_voxel_count + self.cumulative_new_voxels

    def allocate(self, precipitated_salt_mass_kg: float) -> AllocationResult:
        result = allocate_salt_mass(
            precipitated_salt_mass_kg,
            self.residual_mass_kg,
            self.solid_density_kg_m3,
            self.voxel_size_m,
        )
        self.residual_mass_kg = result.residual_mass_kg
        self.cumulative_new_voxels += result.new_solid_voxels
        self.cumulative_budget_mass_kg += precipitated_salt_mass_kg
        return result

    def closure_error_kg(self) -> float:
        """Return allocated plus residual minus cumulative input budget."""

        return (
            self.cumulative_allocated_mass_kg
            + self.residual_mass_kg
            - self.cumulative_budget_mass_kg
        )
