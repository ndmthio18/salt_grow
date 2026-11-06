"""CPU cellular automaton for the bounded Phase 2A/2A.1 prototypes."""

from __future__ import annotations

from collections import deque
import copy
from dataclasses import dataclass
import math
from typing import Literal, Protocol

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
    growth_mode: str
    update_mode: str


@dataclass(frozen=True)
class AttachmentWeights:
    """Uncalibrated relative weights for coordination-dependent sampling."""

    w1: float = 1.0
    w2: float = 1.0
    w3: float = 1.0
    w4: float = 1.0

    def __post_init__(self) -> None:
        for name, value in (
            ("w1", self.w1),
            ("w2", self.w2),
            ("w3", self.w3),
            ("w4", self.w4),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")

    def values_for(self, coordination: NDArray[np.uint8]) -> NDArray[np.float64]:
        """Map coordination 1, 2, 3, and >=4 to relative weights."""

        return np.select(
            [coordination == 1, coordination == 2, coordination == 3],
            [self.w1, self.w2, self.w3],
            default=self.w4,
        ).astype(float)


class CandidateWeightModel(Protocol):
    """Transactional optional weighting model for eligible growth sites."""

    def weights_for_candidates(
        self,
        salt_mask: NDArray[np.bool_],
        coordinates: list[Coordinate],
    ) -> NDArray[np.float64]: ...

    def record_selection(self, coordinate: Coordinate) -> None: ...

    def notify_insertion(self, coordinate: Coordinate) -> None: ...

    def snapshot_state(self) -> object: ...

    def restore_state(self, snapshot: object) -> None: ...


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
        growth_mode: Literal["deterministic", "stochastic"] = "deterministic",
        update_mode: Literal["synchronous", "sequential"] = "synchronous",
        attachment_weights: AttachmentWeights | None = None,
        rng_seed: int = 0,
        candidate_weight_model: CandidateWeightModel | None = None,
    ) -> None:
        state.validate()
        if not isinstance(enforce_wetting, bool):
            raise ValueError("enforce_wetting must be boolean")
        if isinstance(tie_break_seed, bool) or int(tie_break_seed) != tie_break_seed:
            raise ValueError("tie_break_seed must be an integer")
        if growth_mode not in {"deterministic", "stochastic"}:
            raise ValueError("growth_mode must be 'deterministic' or 'stochastic'")
        if update_mode not in {"synchronous", "sequential"}:
            raise ValueError("update_mode must be 'synchronous' or 'sequential'")
        if isinstance(rng_seed, bool) or int(rng_seed) != rng_seed:
            raise ValueError("rng_seed must be an integer")
        self.state = state
        self.enforce_wetting = enforce_wetting
        self.tie_break_seed = int(tie_break_seed)
        self.growth_mode = growth_mode
        self.update_mode = update_mode
        self.attachment_weights = attachment_weights or AttachmentWeights()
        self.rng_seed = int(rng_seed)
        self.rng = np.random.default_rng(self.rng_seed)
        self.candidate_weight_model = candidate_weight_model
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
        self,
        eligible: NDArray[np.bool_],
        *,
        transaction_id: int,
        event_index: int = 0,
    ) -> list[Coordinate]:
        salt_coordination = _face_neighbor_count(self.state.material == Material.SALT)
        coordinates = [
            tuple(int(value) for value in coordinate) for coordinate in np.argwhere(eligible)
        ]
        return sorted(
            coordinates,
            key=lambda coordinate: (
                -int(salt_coordination[coordinate]),
                self._tie_break_priority(
                    coordinate, transaction_id * 1_000_003 + event_index
                ),
                coordinate[0],
                coordinate[1],
                coordinate[2],
            ),
        )

    def _stochastic_candidate(
        self,
        coordinates: list[Coordinate],
        salt_coordination: NDArray[np.uint8],
    ) -> Coordinate:
        coordination = np.array(
            [salt_coordination[coordinate] for coordinate in coordinates],
            dtype=np.uint8,
        )
        weights = self.attachment_weights.values_for(coordination)
        if self.candidate_weight_model is not None:
            weights *= self.candidate_weight_model.weights_for_candidates(
                self.state.material == Material.SALT,
                coordinates,
            )
        probabilities = weights / np.sum(weights)
        index = int(self.rng.choice(len(coordinates), p=probabilities))
        selected = coordinates[index]
        if self.candidate_weight_model is not None:
            self.candidate_weight_model.record_selection(selected)
        return selected

    def _synchronous_selection(
        self,
        eligible: NDArray[np.bool_],
        *,
        requested: int,
        transaction_id: int,
    ) -> tuple[Coordinate, ...]:
        coordinates = [
            tuple(int(value) for value in coordinate)
            for coordinate in np.argwhere(eligible)
        ]
        realized = min(requested, len(coordinates))
        if realized == 0:
            return ()
        if self.growth_mode == "deterministic":
            ordered = self._ordered_candidates(
                eligible, transaction_id=transaction_id
            )
            return tuple(ordered[:realized])

        salt_coordination = _face_neighbor_count(
            self.state.material == Material.SALT
        )
        coordination = np.array(
            [salt_coordination[coordinate] for coordinate in coordinates],
            dtype=np.uint8,
        )
        weights = self.attachment_weights.values_for(coordination)
        if self.candidate_weight_model is not None:
            weights *= self.candidate_weight_model.weights_for_candidates(
                self.state.material == Material.SALT,
                coordinates,
            )
        probabilities = weights / np.sum(weights)
        indices = self.rng.choice(
            len(coordinates), size=realized, replace=False, p=probabilities
        )
        selected = tuple(coordinates[int(index)] for index in indices)
        if self.candidate_weight_model is not None:
            for coordinate in selected:
                self.candidate_weight_model.record_selection(coordinate)
        return selected

    def _sequential_selection(
        self, *, requested: int, transaction_id: int
    ) -> tuple[int, tuple[Coordinate, ...]]:
        initial_candidate_count = int(np.count_nonzero(self.eligible_growth_mask()))
        selected: list[Coordinate] = []
        for event_index in range(requested):
            eligible = self.eligible_growth_mask()
            coordinates = [
                tuple(int(value) for value in coordinate)
                for coordinate in np.argwhere(eligible)
            ]
            if not coordinates:
                break
            if self.growth_mode == "deterministic":
                coordinate = self._ordered_candidates(
                    eligible,
                    transaction_id=transaction_id,
                    event_index=event_index,
                )[0]
            else:
                salt_coordination = _face_neighbor_count(
                    self.state.material == Material.SALT
                )
                coordinate = self._stochastic_candidate(
                    coordinates, salt_coordination
                )
            self._insert_tentative_voxel(coordinate)
            selected.append(coordinate)
        return initial_candidate_count, tuple(selected)

    def _insert_tentative_voxel(self, coordinate: Coordinate) -> None:
        """Insert one uncommitted voxel; the enclosing step owns rollback."""

        if self.state.material[coordinate] != Material.EMPTY:
            raise RuntimeError("tentative growth target is not empty")
        self.state.material[coordinate] = Material.SALT
        if self.candidate_weight_model is not None:
            self.candidate_weight_model.notify_insertion(coordinate)

    def grow(
        self,
        available_precipitation_mass_kg: float,
        *,
        released_solution_mass_kg: float = 0.0,
        unallocated_state: Literal["deferred", "solution"] = "deferred",
    ) -> GrowthStepResult:
        """Execute one transactional mass-limited geometry update."""

        self.assert_invariants()
        material_before = self.state.material.copy()
        rng_state_before = copy.deepcopy(self.rng.bit_generator.state)
        allocator_before = self.allocator.snapshot_state()
        weight_model_before = (
            self.candidate_weight_model.snapshot_state()
            if self.candidate_weight_model is not None
            else None
        )
        selected: tuple[Coordinate, ...] = ()
        try:
            proposal = self.allocator.prepare(
                available_precipitation_mass_kg,
                released_solution_mass_kg=released_solution_mass_kg,
            )
            if self.update_mode == "synchronous":
                eligible = self.eligible_growth_mask()
                candidate_count = int(np.count_nonzero(eligible))
                selected = self._synchronous_selection(
                    eligible,
                    requested=proposal.requested_solid_voxels,
                    transaction_id=proposal.transaction_id,
                )
                for coordinate in selected:
                    self.state.material[coordinate] = Material.SALT
                    if self.candidate_weight_model is not None:
                        self.candidate_weight_model.notify_insertion(coordinate)
            else:
                candidate_count, selected = self._sequential_selection(
                    requested=proposal.requested_solid_voxels,
                    transaction_id=proposal.transaction_id,
                )
            transaction = self.allocator.commit(
                proposal,
                realized_solid_voxels=len(selected),
                unallocated_state=unallocated_state,
            )
        except Exception:
            # ``selected`` can still be empty when sequential selection fails
            # after internal insertions, so restore complete snapshots instead.
            self.state.material[...] = material_before
            self.allocator.restore_state(allocator_before)
            self.rng.bit_generator.state = rng_state_before
            if self.candidate_weight_model is not None:
                self.candidate_weight_model.restore_state(weight_model_before)
            self.assert_invariants()
            raise

        self.assert_invariants()
        return GrowthStepResult(
            transaction=transaction,
            candidate_voxel_count=candidate_count,
            selected_coordinates=selected,
            morphology=self.morphology_diagnostics(),
            growth_mode=self.growth_mode,
            update_mode=self.update_mode,
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
