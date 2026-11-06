"""Conservative conversion from precipitated NaCl mass to CA solid voxels."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal


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


@dataclass(frozen=True)
class AllocationProposal:
    """Non-mutating mass-based request awaiting a geometry decision."""

    transaction_id: int
    incoming_available_precipitation_mass_kg: float
    previous_deferred_mass_kg: float
    released_solution_mass_kg: float
    available_mass_kg: float
    voxel_mass_kg: float
    requested_solid_voxels: int
    requested_solid_mass_kg: float
    subvoxel_mass_kg: float


@dataclass(frozen=True)
class TransactionResult:
    """Committed result after the geometry layer confirms created voxels."""

    transaction_id: int
    available_mass_kg: float
    requested_solid_voxels: int
    realized_solid_voxels: int
    committed_solid_mass_kg: float
    unrealized_requested_mass_kg: float
    subvoxel_mass_kg: float
    unallocated_mass_kg: float
    deferred_precipitation_mass_kg: float
    mass_remaining_in_solution_kg: float
    unallocated_state: str
    status: str


class TransactionalMassToVoxelAllocator:
    """Two-stage ledger for future coupling to a spatial growth algorithm.

    :meth:`prepare` converts mass into a *request* without changing the ledger.
    The Phase 2 geometry layer may create at most that many eligible voxels and
    must then call :meth:`commit` with the actual count.  Only realized voxels
    enter cumulative solid mass.  All other mass remains explicitly either in
    a deferred-precipitation pool or in solution.

    Mass placed in the deferred pool is automatically reconsidered by the next
    transaction.  Mass placed in solution is not reconsidered until a caller
    explicitly releases some or all of it through ``released_solution_mass_kg``.
    This distinction prevents a physical eligibility failure from being
    silently relabeled as precipitated solid.
    """

    def __init__(self, *, solid_density_kg_m3: float, voxel_size_m: float) -> None:
        _positive_finite("solid_density_kg_m3", solid_density_kg_m3)
        _positive_finite("voxel_size_m", voxel_size_m)
        self.solid_density_kg_m3 = solid_density_kg_m3
        self.voxel_size_m = voxel_size_m
        self.deferred_precipitation_mass_kg = 0.0
        self.mass_remaining_in_solution_kg = 0.0
        self.cumulative_input_mass_kg = 0.0
        self.cumulative_committed_solid_mass_kg = 0.0
        self.cumulative_realized_solid_voxels = 0
        self._next_transaction_id = 1
        self._pending: AllocationProposal | None = None

    @property
    def voxel_mass_kg(self) -> float:
        return self.solid_density_kg_m3 * self.voxel_size_m**3

    def prepare(
        self,
        available_precipitation_mass_kg: float,
        *,
        released_solution_mass_kg: float = 0.0,
    ) -> AllocationProposal:
        """Create a voxel request without committing any solid mass."""

        if self._pending is not None:
            raise RuntimeError("a geometry decision is already pending")
        _nonnegative_finite(
            "available_precipitation_mass_kg", available_precipitation_mass_kg
        )
        _nonnegative_finite(
            "released_solution_mass_kg", released_solution_mass_kg
        )
        if released_solution_mass_kg > self.mass_remaining_in_solution_kg:
            raise ValueError("released_solution_mass_kg exceeds mass held in solution")

        available = (
            available_precipitation_mass_kg
            + self.deferred_precipitation_mass_kg
            + released_solution_mass_kg
        )
        requested = int(available // self.voxel_mass_kg)
        requested_mass = requested * self.voxel_mass_kg
        subvoxel_mass = available - requested_mass
        if subvoxel_mass < 0.0 and math.isclose(
            subvoxel_mass, 0.0, abs_tol=math.ulp(max(available, self.voxel_mass_kg))
        ):
            subvoxel_mass = 0.0

        proposal = AllocationProposal(
            transaction_id=self._next_transaction_id,
            incoming_available_precipitation_mass_kg=available_precipitation_mass_kg,
            previous_deferred_mass_kg=self.deferred_precipitation_mass_kg,
            released_solution_mass_kg=released_solution_mass_kg,
            available_mass_kg=available,
            voxel_mass_kg=self.voxel_mass_kg,
            requested_solid_voxels=requested,
            requested_solid_mass_kg=requested_mass,
            subvoxel_mass_kg=subvoxel_mass,
        )
        self._pending = proposal
        return proposal

    def commit(
        self,
        proposal: AllocationProposal,
        *,
        realized_solid_voxels: int,
        unallocated_state: Literal["deferred", "solution"] = "deferred",
    ) -> TransactionResult:
        """Commit only geometry-confirmed voxels and retain all other mass."""

        if self._pending is None or proposal is not self._pending:
            raise RuntimeError("proposal is not the allocator's pending transaction")
        if isinstance(realized_solid_voxels, bool):
            raise ValueError("realized_solid_voxels must be a non-negative integer")
        if int(realized_solid_voxels) != realized_solid_voxels or realized_solid_voxels < 0:
            raise ValueError("realized_solid_voxels must be a non-negative integer")
        realized = int(realized_solid_voxels)
        if realized > proposal.requested_solid_voxels:
            raise ValueError("geometry cannot realize more voxels than were requested")
        if unallocated_state not in {"deferred", "solution"}:
            raise ValueError("unallocated_state must be 'deferred' or 'solution'")

        committed_mass = realized * proposal.voxel_mass_kg
        unrealized_requested_mass = (
            proposal.requested_solid_voxels - realized
        ) * proposal.voxel_mass_kg
        unallocated_mass = proposal.available_mass_kg - committed_mass
        solution_after_release = (
            self.mass_remaining_in_solution_kg - proposal.released_solution_mass_kg
        )

        if unallocated_state == "deferred":
            deferred = unallocated_mass
            solution = solution_after_release
        else:
            deferred = 0.0
            solution = solution_after_release + unallocated_mass

        if proposal.requested_solid_voxels == 0:
            status = "zero_growth_insufficient_mass"
        elif realized == 0:
            status = "no_spatial_capacity"
        elif realized < proposal.requested_solid_voxels:
            status = "partial_spatial_capacity"
        else:
            status = "committed"

        self.deferred_precipitation_mass_kg = deferred
        self.mass_remaining_in_solution_kg = solution
        self.cumulative_input_mass_kg += (
            proposal.incoming_available_precipitation_mass_kg
        )
        self.cumulative_committed_solid_mass_kg += committed_mass
        self.cumulative_realized_solid_voxels += realized
        self._pending = None
        self._next_transaction_id += 1

        return TransactionResult(
            transaction_id=proposal.transaction_id,
            available_mass_kg=proposal.available_mass_kg,
            requested_solid_voxels=proposal.requested_solid_voxels,
            realized_solid_voxels=realized,
            committed_solid_mass_kg=committed_mass,
            unrealized_requested_mass_kg=unrealized_requested_mass,
            subvoxel_mass_kg=proposal.subvoxel_mass_kg,
            unallocated_mass_kg=unallocated_mass,
            deferred_precipitation_mass_kg=deferred,
            mass_remaining_in_solution_kg=solution,
            unallocated_state=unallocated_state,
            status=status,
        )

    def closure_error_kg(self) -> float:
        """Return solid plus retained pools minus all supplied mass."""

        return (
            self.cumulative_committed_solid_mass_kg
            + self.deferred_precipitation_mass_kg
            + self.mass_remaining_in_solution_kg
            - self.cumulative_input_mass_kg
        )
