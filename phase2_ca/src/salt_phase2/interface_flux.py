"""Conservative interface-area allocation for Phase 2A.3b-3.

The continuum law in this module is deliberately independent of CA candidate
counts.  Integer voxels are only a transactional geometric realization of the
already allocated physical patch volumes.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import math
from time import perf_counter
from typing import Hashable, Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree

from .anisotropic_growth import OrientationCandidateWeightModel
from .surface_orientation import estimate_plane_fit_normals, signed_distance_from_mask


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class PatchAllocation:
    """One exact continuum allocation step in SI units."""

    total_volume_m3: float
    patch_volumes_m3: FloatArray
    normal_advances_m: FloatArray
    weights_m2: FloatArray
    closure_error_m3: float


def allocate_interface_volume(
    total_volume_m3: float,
    patch_areas_m2: ArrayLike,
    kinetic_factors: ArrayLike,
) -> PatchAllocation:
    """Allocate volume in proportion to ``area * kappa``.

    The last addition of the floating residual makes the returned vector sum
    to the supplied floating-point volume exactly.  This is accounting, not a
    physical redistribution or calibration.
    """

    areas = np.asarray(patch_areas_m2, dtype=float)
    kappa = np.asarray(kinetic_factors, dtype=float)
    if areas.ndim != 1 or areas.size == 0 or kappa.shape != areas.shape:
        raise ValueError("areas and kinetic factors must be matching non-empty vectors")
    if not math.isfinite(total_volume_m3) or total_volume_m3 < 0.0:
        raise ValueError("total_volume_m3 must be finite and nonnegative")
    if np.any(~np.isfinite(areas)) or np.any(areas <= 0.0):
        raise ValueError("patch areas must be finite and positive")
    if np.any(~np.isfinite(kappa)) or np.any(kappa <= 0.0):
        raise ValueError("kinetic factors must be finite and positive")
    weights = areas * kappa
    volumes = total_volume_m3 * weights / math.fsum(float(v) for v in weights)
    residual = total_volume_m3 - math.fsum(float(v) for v in volumes)
    volumes[int(np.argmax(weights))] += residual
    closure = math.fsum(float(v) for v in volumes) - total_volume_m3
    return PatchAllocation(
        total_volume_m3=float(total_volume_m3),
        patch_volumes_m3=volumes,
        normal_advances_m=volumes / areas,
        weights_m2=weights,
        closure_error_m3=float(closure),
    )


@dataclass(frozen=True)
class SurfaceElements:
    """Non-overlapping surface elements linked to exterior candidate voxels."""

    areas_m2: FloatArray
    centroids_m: FloatArray
    normals: FloatArray
    candidate_coordinates: IntArray
    normal_valid: NDArray[np.bool_]
    method: str
    runtime_s: float

    @property
    def total_area_m2(self) -> float:
        return math.fsum(float(value) for value in self.areas_m2)


_CUBE_CORNERS = np.asarray(
    [[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=int
)
_TETRAHEDRA = (
    (0, 1, 3, 7),
    (0, 3, 2, 7),
    (0, 2, 6, 7),
    (0, 6, 4, 7),
    (0, 4, 5, 7),
    (0, 5, 1, 7),
)
_TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


def _triangulate_polygon(points: FloatArray, outward: FloatArray) -> list[FloatArray]:
    centre = np.mean(points, axis=0)
    normal = outward / np.linalg.norm(outward)
    helper = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(helper, normal))) > 0.8:
        helper = np.array([0.0, 1.0, 0.0])
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    second = np.cross(normal, first)
    offsets = points - centre
    angles = np.arctan2(offsets @ second, offsets @ first)
    ordered = points[np.argsort(angles)]
    return [np.asarray([ordered[0], ordered[i], ordered[i + 1]]) for i in range(1, len(ordered) - 1)]


def _candidate_coordinates(mask: NDArray[np.bool_]) -> IntArray:
    candidates: set[tuple[int, int, int]] = set()
    for axis in range(3):
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        lhs = mask[tuple(left)]
        rhs = mask[tuple(right)]
        for index in np.argwhere(lhs != rhs):
            coordinate = index.copy()
            if lhs[tuple(index)]:
                coordinate[axis] += 1
            candidates.add(tuple(int(v) for v in coordinate))
    if not candidates:
        raise ValueError("mask has no internal interface candidates")
    return np.asarray(sorted(candidates), dtype=np.int64)


def triangulated_interface(mask: ArrayLike, voxel_size_m: float) -> SurfaceElements:
    """Reconstruct the internal signed-distance zero surface with tetrahedra.

    Cell centres are interpolation vertices.  Domain-boundary faces are not
    invented, so a half-space benchmark returns only its internal plane.
    Every triangle belongs to exactly one tetrahedron and is assigned to its
    nearest exterior face-connected candidate voxel.
    """

    started = perf_counter()
    solid = np.asarray(mask, dtype=bool)
    if solid.ndim != 3 or min(solid.shape) < 3 or not np.any(solid) or np.all(solid):
        raise ValueError("mask must be a mixed three-dimensional array")
    if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
        raise ValueError("voxel_size_m must be finite and positive")
    origin = (-0.5 * np.asarray(solid.shape, dtype=float) + 0.5) * voxel_size_m
    level_set = signed_distance_from_mask(solid, voxel_size_m)
    triangles: list[FloatArray] = []
    normals: list[FloatArray] = []
    areas: list[float] = []
    active = np.argwhere(
        np.logical_and(
            np.maximum.reduce([
                solid[:-1, :-1, :-1], solid[1:, :-1, :-1],
                solid[:-1, 1:, :-1], solid[:-1, :-1, 1:],
                solid[1:, 1:, :-1], solid[1:, :-1, 1:],
                solid[:-1, 1:, 1:], solid[1:, 1:, 1:],
            ]),
            ~np.minimum.reduce([
                solid[:-1, :-1, :-1], solid[1:, :-1, :-1],
                solid[:-1, 1:, :-1], solid[:-1, :-1, 1:],
                solid[1:, 1:, :-1], solid[1:, :-1, 1:],
                solid[:-1, 1:, 1:], solid[1:, 1:, 1:],
            ]),
        )
    )
    for base in active:
        cube_indices = base + _CUBE_CORNERS
        cube_values = level_set[tuple(cube_indices.T)]
        cube_solid = solid[tuple(cube_indices.T)]
        cube_points = origin + cube_indices * voxel_size_m
        for tetrahedron in _TETRAHEDRA:
            ids = np.asarray(tetrahedron, dtype=int)
            values = cube_values[ids]
            inside = values <= 0.0
            if np.all(inside) or not np.any(inside):
                continue
            points = cube_points[ids]
            crossings = []
            for first, second in _TET_EDGES:
                if inside[first] != inside[second]:
                    fraction = values[first] / (values[first] - values[second])
                    crossings.append(points[first] + fraction * (points[second] - points[first]))
            crossings = np.unique(np.asarray(crossings), axis=0)
            tetra_solid = cube_solid[ids]
            outward = np.mean(points[~tetra_solid], axis=0) - np.mean(points[tetra_solid], axis=0)
            for triangle in _triangulate_polygon(crossings, outward):
                vector = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
                magnitude = float(np.linalg.norm(vector))
                if magnitude <= 1.0e-30:
                    continue
                normal = vector / magnitude
                if float(np.dot(normal, outward)) < 0.0:
                    normal = -normal
                triangles.append(triangle)
                normals.append(normal)
                areas.append(0.5 * magnitude)
    if not triangles:
        raise RuntimeError("triangulation produced no surface elements")
    centroids = np.asarray([np.mean(triangle, axis=0) for triangle in triangles])
    candidates = _candidate_coordinates(solid)
    candidate_centres = origin + candidates * voxel_size_m
    nearest = cKDTree(candidate_centres).query(centroids, k=1)[1]
    return SurfaceElements(
        areas_m2=np.asarray(areas),
        centroids_m=centroids,
        normals=np.asarray(normals),
        candidate_coordinates=candidates[np.asarray(nearest, dtype=int)],
        normal_valid=np.ones(len(areas), dtype=bool),
        method="marching_tetrahedra_signed_distance_zero",
        runtime_s=perf_counter() - started,
    )


def corrected_digital_interface(
    mask: ArrayLike,
    voxel_size_m: float,
    fitting_radius_m: float,
) -> SurfaceElements:
    """Partition unique internal salt/empty faces and apply L1 area correction.

    Each digital face occurs once and links to exactly one empty voxel.  A
    multi-face candidate therefore owns several distinct face elements; it is
    not silently reclassified as a physical crystallographic corner.
    """

    started = perf_counter()
    solid = np.asarray(mask, dtype=bool)
    if solid.ndim != 3 or min(solid.shape) < 3 or not np.any(solid) or np.all(solid):
        raise ValueError("mask must be a mixed three-dimensional array")
    origin = -0.5 * np.asarray(solid.shape, dtype=float) * voxel_size_m
    centres: list[FloatArray] = []
    candidates: list[IntArray] = []
    for axis in range(3):
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        lhs = solid[tuple(left)]
        rhs = solid[tuple(right)]
        for index in np.argwhere(lhs != rhs):
            centre = origin + (index.astype(float) + 0.5) * voxel_size_m
            centre[axis] += 0.5 * voxel_size_m
            candidate = index.copy()
            if lhs[tuple(index)]:
                candidate[axis] += 1
            centres.append(centre)
            candidates.append(candidate.astype(np.int64))
    points = np.asarray(centres)
    estimate = estimate_plane_fit_normals(
        solid, voxel_size_m, fitting_radius_m, points_m=points
    )
    l1 = np.sum(np.abs(estimate.normals), axis=1)
    areas = np.full(points.shape[0], voxel_size_m**2, dtype=float)
    areas[estimate.valid] /= l1[estimate.valid]
    return SurfaceElements(
        areas_m2=areas,
        centroids_m=points,
        normals=estimate.normals,
        candidate_coordinates=np.asarray(candidates, dtype=np.int64),
        normal_valid=estimate.valid,
        method="unique_faces_l1_corrected_plane_fit",
        runtime_s=perf_counter() - started,
    )


@dataclass(frozen=True)
class CandidatePatch:
    coordinate: tuple[int, int, int]
    area_m2: float
    normal: tuple[float, float, float]
    element_count: int
    normal_valid: bool


def bundle_elements_by_candidate(elements: SurfaceElements) -> tuple[CandidatePatch, ...]:
    """Create disjoint candidate-owned patch bundles without duplicating area."""

    groups: dict[tuple[int, int, int], list[int]] = {}
    for index, coordinate in enumerate(elements.candidate_coordinates):
        groups.setdefault(tuple(int(v) for v in coordinate), []).append(index)
    result = []
    for coordinate in sorted(groups):
        indices = np.asarray(groups[coordinate], dtype=int)
        area = math.fsum(float(elements.areas_m2[i]) for i in indices)
        valid = elements.normal_valid[indices]
        vector = np.sum(
            elements.normals[indices][valid] * elements.areas_m2[indices][valid, None],
            axis=0,
        ) if np.any(valid) else np.full(3, np.nan)
        magnitude = float(np.linalg.norm(vector))
        normal = vector / magnitude if magnitude > 1.0e-30 else np.full(3, np.nan)
        result.append(CandidatePatch(
            coordinate=coordinate,
            area_m2=area,
            normal=tuple(float(v) for v in normal),
            element_count=len(indices),
            normal_valid=bool(np.all(valid) and np.all(np.isfinite(normal))),
        ))
    if not math.isclose(
        math.fsum(patch.area_m2 for patch in result),
        elements.total_area_m2,
        rel_tol=2.0e-15,
        abs_tol=0.0,
    ):
        raise RuntimeError("candidate bundling changed total interface area")
    return tuple(result)


@dataclass(frozen=True)
class LedgerUpdate:
    transaction_id: int
    incoming_mass_kg: float
    allocated_mass_kg: tuple[float, ...]
    realized_voxels: int
    realized_mass_kg: float
    reservoir_mass_kg: float
    retained_mass_kg: float
    closure_error_kg: float


class ConservativePatchLedger:
    """Transactional patch reservoirs with unique candidate consumption."""

    def __init__(self, density_kg_m3: float, voxel_size_m: float) -> None:
        if not math.isfinite(density_kg_m3) or density_kg_m3 <= 0.0:
            raise ValueError("density must be finite and positive")
        if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
            raise ValueError("voxel size must be finite and positive")
        self.density_kg_m3 = float(density_kg_m3)
        self.voxel_size_m = float(voxel_size_m)
        self.voxel_mass_kg = self.density_kg_m3 * self.voxel_size_m**3
        self.cumulative_input_mass_kg = 0.0
        self.cumulative_allocated_mass_kg = 0.0
        self.cumulative_realized_mass_kg = 0.0
        self.cumulative_retained_mass_kg = 0.0
        self.patch_reservoirs_kg: dict[Hashable, float] = {}
        self.consumed_candidates: set[Hashable] = set()
        self.transaction_id = 0

    def snapshot_state(self) -> object:
        return copy.deepcopy(self.__dict__)

    def restore_state(self, snapshot: object) -> None:
        if not isinstance(snapshot, dict):
            raise TypeError("invalid ledger snapshot")
        self.__dict__ = copy.deepcopy(snapshot)

    def closure_error_kg(self) -> float:
        return self.cumulative_input_mass_kg - (
            self.cumulative_realized_mass_kg
            + self.cumulative_retained_mass_kg
            + math.fsum(self.patch_reservoirs_kg.values())
        )

    def assert_closure(self) -> None:
        scale = max(self.cumulative_input_mass_kg, self.voxel_mass_kg)
        if abs(self.closure_error_kg()) > 64.0 * math.ulp(scale):
            raise RuntimeError("patch mass ledger does not close")
        if not math.isclose(
            self.cumulative_allocated_mass_kg,
            self.cumulative_input_mass_kg,
            rel_tol=0.0,
            abs_tol=64.0 * math.ulp(scale),
        ):
            raise RuntimeError("allocated mass differs from incoming mass")

    def update(
        self,
        incoming_mass_kg: float,
        patch_ids: Sequence[Hashable],
        patch_areas_m2: ArrayLike,
        kinetic_factors: ArrayLike,
        candidate_queues: Sequence[Sequence[Hashable]],
        *,
        retain_if_exhausted: bool = True,
        inject_failure_after_voxels: int | None = None,
    ) -> LedgerUpdate:
        if len(patch_ids) != len(candidate_queues) or len(set(patch_ids)) != len(patch_ids):
            raise ValueError("patch ids must be unique and match candidate queues")
        if not math.isfinite(incoming_mass_kg) or incoming_mass_kg < 0.0:
            raise ValueError("incoming mass must be finite and nonnegative")
        before = self.snapshot_state()
        realized = 0
        try:
            allocation = allocate_interface_volume(
                incoming_mass_kg / self.density_kg_m3,
                patch_areas_m2,
                kinetic_factors,
            )
            allocated = allocation.patch_volumes_m3 * self.density_kg_m3
            self.transaction_id += 1
            self.cumulative_input_mass_kg += incoming_mass_kg
            self.cumulative_allocated_mass_kg += math.fsum(float(v) for v in allocated)
            for patch_id, mass in zip(patch_ids, allocated, strict=True):
                self.patch_reservoirs_kg[patch_id] = (
                    self.patch_reservoirs_kg.get(patch_id, 0.0) + float(mass)
                )
            for patch_id, queue in zip(patch_ids, candidate_queues, strict=True):
                available = [candidate for candidate in queue if candidate not in self.consumed_candidates]
                while (
                    self.patch_reservoirs_kg[patch_id]
                    >= self.voxel_mass_kg * (1.0 - 8.0 * np.finfo(float).eps)
                    and available
                ):
                    candidate = available.pop(0)
                    self.consumed_candidates.add(candidate)
                    self.patch_reservoirs_kg[patch_id] -= self.voxel_mass_kg
                    if self.patch_reservoirs_kg[patch_id] < 0.0:
                        self.patch_reservoirs_kg[patch_id] = 0.0
                    self.cumulative_realized_mass_kg += self.voxel_mass_kg
                    realized += 1
                    if inject_failure_after_voxels == realized:
                        raise RuntimeError("injected patch realization failure")
                if retain_if_exhausted and not available and self.patch_reservoirs_kg[patch_id] >= self.voxel_mass_kg:
                    retained = self.patch_reservoirs_kg[patch_id]
                    self.cumulative_retained_mass_kg += retained
                    self.patch_reservoirs_kg[patch_id] = 0.0
            self.assert_closure()
        except Exception:
            self.restore_state(before)
            raise
        return LedgerUpdate(
            transaction_id=self.transaction_id,
            incoming_mass_kg=float(incoming_mass_kg),
            allocated_mass_kg=tuple(float(v) for v in allocated),
            realized_voxels=realized,
            realized_mass_kg=realized * self.voxel_mass_kg,
            reservoir_mass_kg=math.fsum(self.patch_reservoirs_kg.values()),
            retained_mass_kg=self.cumulative_retained_mass_kg,
            closure_error_kg=self.closure_error_kg(),
        )


def frozen_front_realization(
    total_volume_m3: float,
    patch_areas_m2: ArrayLike,
    kinetic_factors: ArrayLike,
    voxel_size_m: float,
    candidate_capacities: ArrayLike,
) -> dict[str, FloatArray | float | int]:
    """Deterministic one-update integer realization of frozen physical patches."""

    allocation = allocate_interface_volume(total_volume_m3, patch_areas_m2, kinetic_factors)
    capacities = np.asarray(candidate_capacities, dtype=np.int64)
    if capacities.shape != allocation.patch_volumes_m3.shape or np.any(capacities < 0):
        raise ValueError("candidate capacities must match patches and be nonnegative")
    voxel_volume = float(voxel_size_m) ** 3
    requested = np.floor(
        allocation.patch_volumes_m3 / voxel_volume + 8.0 * np.finfo(float).eps
    ).astype(np.int64)
    realized = np.minimum(requested, capacities)
    realized_volumes = realized.astype(float) * voxel_volume
    reservoirs = allocation.patch_volumes_m3 - requested.astype(float) * voxel_volume
    retained = (requested - realized).astype(float) * voxel_volume
    closure = total_volume_m3 - math.fsum(
        float(v) for v in np.concatenate((realized_volumes, reservoirs, retained))
    )
    return {
        "realized_voxels": realized,
        "requested_voxels": requested,
        "realized_volumes_m3": realized_volumes,
        "reservoir_volumes_m3": reservoirs,
        "retained_volumes_m3": retained,
        "voxel_only_advances_m": realized_volumes / np.asarray(patch_areas_m2, dtype=float),
        "allocated_advances_m": allocation.normal_advances_m,
        "closure_error_m3": float(closure),
        "capacity_limited_patches": int(np.count_nonzero(requested > capacities)),
    }


class AreaCandidateWeightModel(OrientationCandidateWeightModel):
    """CA adapter whose event probability is corrected patch area times kappa.

    This adapter is only the bounded single-grain realization requested in
    Phase 2A.3b-3.  The independent continuum allocation and patch ledger above
    remain the governing closure; this class does not turn ``kappa`` into a
    calibrated material rate constant.
    """

    def weights_for_candidates(
        self,
        salt_mask: NDArray[np.bool_],
        coordinates: list[tuple[int, int, int]],
    ) -> FloatArray:
        self.weight_calls += 1
        if self.update_strategy == "full":
            records = self._compute(salt_mask, coordinates)
            self._cache = {record.coordinate: record for record in records}
        else:
            missing = [coordinate for coordinate in coordinates if coordinate not in self._cache]
            for record in self._compute(salt_mask, missing):
                self._cache[record.coordinate] = record
            records = tuple(self._cache[coordinate] for coordinate in coordinates)
        self._last_records = {record.coordinate: record for record in records}
        self.candidate_evaluations += len(records)
        weights = []
        for record in records:
            self.class_counts[record.site_class] += 1
            self.anisotropy_eligible_evaluations += int(record.anisotropy_eligible)
            normal = np.asarray(record.combined_normal, dtype=float)
            if np.all(np.isfinite(normal)):
                projection = float(np.sum(np.abs(normal)))
                corrected_area = (
                    record.adjacent_salt_faces * self.voxel_size_m**2 / projection
                )
            else:
                corrected_area = record.adjacent_salt_faces * self.voxel_size_m**2
            weights.append(corrected_area * record.anisotropy_weight)
        result = np.asarray(weights, dtype=float)
        if np.any(~np.isfinite(result)) or np.any(result <= 0.0):
            raise RuntimeError("candidate area-flux weights must be finite and positive")
        return result
