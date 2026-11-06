"""Deterministic CPU cellular automaton for the bounded Phase 2A prototype."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from salt_phase1.voxel import TransactionResult, TransactionalMassToVoxelAllocator

from .state import Coordinate, Material, PhysicalState


_FACE_OFFSETS: tuple[Coordinate, ...] = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)


def _face_neighbor_count(mask: NDArray[np.bool_]) -> NDArray[np.uint8]:
    """Count non-periodic six-face neighbors for every cell."""

    count = np.zeros(mask.shape, dtype=np.uint8)
    count[1:, :, :] += mask[:-1, :, :]
    count[:-1, :, :] += mask[1:, :, :]
    count[:, 1:, :] += mask[:, :-1, :]
    count[:, :-1, :] += mask[:, 1:, :]
    count[:, :, 1:] += mask[:, :, :-1]
    count[:, :, :-1] += mask[:, :, 1:]
    return count


@dataclass(frozen=True)
class MorphologyDiagnostics:
    seed_salt_voxels: int
    new_salt_voxels: int
    total_salt_voxels: int
    supply_connected_salt_voxels: int
    eligible_growth_voxels: int
    footprint_voxels: int
    occupied_salt_layers: int
    maximum_height_voxels_above_substrate: int
    seed_salt_mass_kg: float
    new_salt_mass_kg: float
    total_salt_mass_kg: float
    retained_unallocated_mass_kg: float
    mass_closure_error_kg: float


@dataclass(frozen=True)
class GrowthStepResult:
    transaction: TransactionResult
    candidate_voxel_count: int
    selected_coordinates: tuple[Coordinate, ...]
    morphology: MorphologyDiagnostics


class SaltGrowthCA:
    """Minimal mass-conservative CA with a conceptual surface-film supply rule.

    With ``enforce_wetting=True``, salt is supply-connected only when it is in
    the six-neighbor salt component of a liquid-accessible seed.  With the
    switch disabled, all salt may act as a parent; this is an intentionally
    nonphysical control, not a second wetting model.
    """

    def __init__(
        self,
        state: PhysicalState,
        *,
        solid_density_kg_m3: float,
        voxel_size_m: float,
        enforce_wetting: bool = True,
        tie_break_seed: int = 0,
    ) -> None:
        state.validate()
        if not isinstance(enforce_wetting, bool):
            raise ValueError("enforce_wetting must be boolean")
        if isinstance(tie_break_seed, bool) or int(tie_break_seed) != tie_break_seed:
            raise ValueError("tie_break_seed must be an integer")
        self.state = state
        self.enforce_wetting = enforce_wetting
        self.tie_break_seed = int(tie_break_seed)
        self.allocator = TransactionalMassToVoxelAllocator(
            solid_density_kg_m3=solid_density_kg_m3,
            voxel_size_m=voxel_size_m,
        )
        self._initial_seed_voxel_count = state.initial_seed_voxel_count
        initial_salt_count = int(np.count_nonzero(state.material == Material.SALT))
        if initial_salt_count != self._initial_seed_voxel_count:
            raise ValueError("initial salt voxels must all be declared seeds")

    @property
    def voxel_mass_kg(self) -> float:
        return self.allocator.voxel_mass_kg

    @property
    def initial_seed_mass_kg(self) -> float:
        return self._initial_seed_voxel_count * self.voxel_mass_kg

    def supply_connected_salt_mask(self) -> NDArray[np.bool_]:
        """Return salt connected through six-face contacts to a wetted seed."""

        salt = self.state.material == Material.SALT
        anchors = salt & self.state.seed_mask & self.state.liquid_accessible
        connected = np.zeros(self.state.shape, dtype=bool)
        queue: deque[Coordinate] = deque(
            tuple(int(value) for value in coordinate) for coordinate in np.argwhere(anchors)
        )
        for coordinate in queue:
            connected[coordinate] = True

        shape = self.state.shape
        while queue:
            x, y, z = queue.popleft()
            for dx, dy, dz in _FACE_OFFSETS:
                neighbor = x + dx, y + dy, z + dz
                if not (
                    0 <= neighbor[0] < shape[0]
                    and 0 <= neighbor[1] < shape[1]
                    and 0 <= neighbor[2] < shape[2]
                ):
                    continue
                if salt[neighbor] and not connected[neighbor]:
                    connected[neighbor] = True
                    queue.append(neighbor)
        return connected

    def eligible_growth_mask(self) -> NDArray[np.bool_]:
        """Return empty sites adjacent to a permitted salt parent.

        The target's ``liquid_accessible`` value is intentionally not used: an
        air voxel can grow when an adjacent salt surface is film supplied.
        """

        if self.enforce_wetting:
            parent_salt = self.supply_connected_salt_mask()
        else:
            parent_salt = self.state.material == Material.SALT
        return (self.state.material == Material.EMPTY) & (
            _face_neighbor_count(parent_salt) > 0
        )

    def _tie_break_priority(self, coordinate: Coordinate, transaction_id: int) -> int:
        """Return a reproducible mixed coordinate priority without axis ordering."""

        x, y, z = coordinate
        return (
            x * 73_856_093
            ^ y * 19_349_663
            ^ z * 83_492_791
            ^ transaction_id * 2_654_435_761
            ^ self.tie_break_seed
        ) & 0xFFFFFFFF

    def _ordered_candidates(
        self, eligible: NDArray[np.bool_], *, transaction_id: int
    ) -> list[Coordinate]:
        salt_coordination = _face_neighbor_count(self.state.material == Material.SALT)
        coordinates = [
            tuple(int(value) for value in coordinate) for coordinate in np.argwhere(eligible)
        ]
        return sorted(
            coordinates,
            key=lambda coordinate: (
                -int(salt_coordination[coordinate]),
                self._tie_break_priority(coordinate, transaction_id),
                coordinate[0],
                coordinate[1],
                coordinate[2],
            ),
        )

    def grow(
        self,
        available_precipitation_mass_kg: float,
        *,
        released_solution_mass_kg: float = 0.0,
        unallocated_state: Literal["deferred", "solution"] = "deferred",
    ) -> GrowthStepResult:
        """Execute one transactional mass-limited geometry update."""

        self.assert_invariants()
        proposal = self.allocator.prepare(
            available_precipitation_mass_kg,
            released_solution_mass_kg=released_solution_mass_kg,
        )
        eligible = self.eligible_growth_mask()
        candidates = self._ordered_candidates(
            eligible, transaction_id=proposal.transaction_id
        )
        selected = tuple(candidates[: proposal.requested_solid_voxels])
        for coordinate in selected:
            self.state.material[coordinate] = Material.SALT

        try:
            transaction = self.allocator.commit(
                proposal,
                realized_solid_voxels=len(selected),
                unallocated_state=unallocated_state,
            )
        except Exception:
            for coordinate in selected:
                self.state.material[coordinate] = Material.EMPTY
            raise

        self.assert_invariants()
        return GrowthStepResult(
            transaction=transaction,
            candidate_voxel_count=len(candidates),
            selected_coordinates=selected,
            morphology=self.morphology_diagnostics(),
        )

    def morphology_diagnostics(self) -> MorphologyDiagnostics:
        salt = self.state.material == Material.SALT
        salt_coordinates = np.argwhere(salt)
        total_salt = int(salt_coordinates.shape[0])
        new_salt = self.allocator.cumulative_realized_solid_voxels
        connected = int(np.count_nonzero(self.supply_connected_salt_mask()))
        eligible = int(np.count_nonzero(self.eligible_growth_mask()))
        if total_salt:
            footprint = len({(int(x), int(y)) for x, y, _ in salt_coordinates})
            occupied_layers = len({int(z) for _, _, z in salt_coordinates})
            maximum_height = (
                int(np.max(salt_coordinates[:, 2])) - self.state.substrate_layers + 1
            )
        else:
            footprint = 0
            occupied_layers = 0
            maximum_height = 0
        retained = (
            self.allocator.deferred_precipitation_mass_kg
            + self.allocator.mass_remaining_in_solution_kg
        )
        return MorphologyDiagnostics(
            seed_salt_voxels=self._initial_seed_voxel_count,
            new_salt_voxels=new_salt,
            total_salt_voxels=total_salt,
            supply_connected_salt_voxels=connected,
            eligible_growth_voxels=eligible,
            footprint_voxels=footprint,
            occupied_salt_layers=occupied_layers,
            maximum_height_voxels_above_substrate=maximum_height,
            seed_salt_mass_kg=self.initial_seed_mass_kg,
            new_salt_mass_kg=self.allocator.cumulative_committed_solid_mass_kg,
            total_salt_mass_kg=(
                self.initial_seed_mass_kg
                + self.allocator.cumulative_committed_solid_mass_kg
            ),
            retained_unallocated_mass_kg=retained,
            mass_closure_error_kg=self.allocator.closure_error_kg(),
        )

    def assert_invariants(self) -> None:
        self.state.validate()
        salt_voxels = int(np.count_nonzero(self.state.material == Material.SALT))
        expected = (
            self._initial_seed_voxel_count
            + self.allocator.cumulative_realized_solid_voxels
        )
        if salt_voxels != expected:
            raise RuntimeError(
                "salt geometry is inconsistent with seed plus committed mass ledger"
            )
        tolerance = 32.0 * math.ulp(
            max(self.allocator.cumulative_input_mass_kg, self.voxel_mass_kg)
        )
        if abs(self.allocator.closure_error_kg()) > tolerance:
            raise RuntimeError("transactional salt-mass ledger does not close")
