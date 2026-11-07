"""Dynamic conservative interface-flux realization for Phase 2A.3b-4.

Incoming volume is integrated over distinct corrected digital face elements.
Fractional deposited solid belongs persistently to a physical candidate
location.  Geometry changes only when that owner reaches one voxel volume.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import math
from time import perf_counter
from typing import Hashable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree

from .interface_flux import (
    CandidatePatch,
    SurfaceElements,
    bundle_elements_by_candidate,
    corrected_digital_interface,
)
from .surface_orientation import cubic_attachment_weight


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]
Coordinate = tuple[int, int, int]
Owner = tuple[str, Coordinate]
FaceKey = tuple[int, int, int, int]


class LocalCorrectedInterfaceCache:
    """Unique-face plane-fit patches with conservative local invalidation."""

    def __init__(self, mask: BoolArray, voxel_size_m: float, fitting_radius_m: float, version: int) -> None:
        self.shape = tuple(int(v) for v in mask.shape)
        self.voxel_size_m = float(voxel_size_m)
        self.fitting_radius_m = float(fitting_radius_m)
        self.version = int(version)
        self.geometry: dict[FaceKey, tuple[FloatArray, Coordinate, FloatArray]] = {}
        self.normals: dict[FaceKey, FloatArray] = {}
        self.areas_m2: dict[FaceKey, float] = {}
        self.valid: dict[FaceKey, bool] = {}
        self.recomputed_elements = 0
        self._build_geometry(mask)
        self._recompute(set(self.geometry))

    def _origin(self) -> FloatArray:
        return -0.5 * np.asarray(self.shape, dtype=float) * self.voxel_size_m

    def _face_geometry(
        self, mask: BoolArray, key: FaceKey
    ) -> tuple[FloatArray, Coordinate, FloatArray] | None:
        axis, i, j, k = key
        lower = np.asarray([i, j, k], dtype=int)
        upper = lower.copy()
        upper[axis] += 1
        if np.any(lower < 0) or np.any(upper >= np.asarray(self.shape)):
            return None
        lhs = bool(mask[tuple(lower)])
        rhs = bool(mask[tuple(upper)])
        if lhs == rhs:
            return None
        centre = self._origin() + (lower.astype(float) + 0.5) * self.voxel_size_m
        centre[axis] += 0.5 * self.voxel_size_m
        digital_normal = np.zeros(3, dtype=float)
        if lhs:
            candidate = tuple(int(v) for v in upper)
            digital_normal[axis] = 1.0
        else:
            candidate = tuple(int(v) for v in lower)
            digital_normal[axis] = -1.0
        return centre, candidate, digital_normal

    def _build_geometry(self, mask: BoolArray) -> None:
        for axis in range(3):
            left = [slice(None)] * 3
            right = [slice(None)] * 3
            left[axis] = slice(0, -1)
            right[axis] = slice(1, None)
            transitions = mask[tuple(left)] != mask[tuple(right)]
            for index in np.argwhere(transitions):
                key = (axis, *(int(v) for v in index))
                geometry = self._face_geometry(mask, key)
                if geometry is not None:
                    self.geometry[key] = geometry

    def _recompute(self, keys: set[FaceKey]) -> None:
        if not self.geometry:
            return
        ordered = sorted(self.geometry)
        points = np.asarray([self.geometry[key][0] for key in ordered])
        tree = cKDTree(points)
        index = {key: position for position, key in enumerate(ordered)}
        sigma = 0.5 * self.fitting_radius_m
        for key in sorted(keys & set(self.geometry)):
            point = self.geometry[key][0]
            neighbours = tree.query_ball_point(point, self.fitting_radius_m)
            if len(neighbours) < 6:
                self.normals[key] = np.full(3, np.nan)
                self.areas_m2[key] = self.voxel_size_m**2
                self.valid[key] = False
                continue
            local = points[np.asarray(neighbours, dtype=int)]
            offsets = local - point
            distances = np.linalg.norm(offsets, axis=1)
            weights = np.exp(-0.5 * (distances / sigma) ** 2)
            weight_sum = float(np.sum(weights))
            centroid = np.sum(local * weights[:, None], axis=0) / weight_sum
            centred = local - centroid
            covariance = (centred * weights[:, None]).T @ centred / weight_sum
            eigenvalues, eigenvectors = np.linalg.eigh(covariance)
            if not np.all(np.isfinite(eigenvalues)) or eigenvalues[1] <= 1.0e-24:
                self.normals[key] = np.full(3, np.nan)
                self.areas_m2[key] = self.voxel_size_m**2
                self.valid[key] = False
                continue
            normal = eigenvectors[:, 0]
            digital_normal = self.geometry[key][2]
            if float(np.dot(normal, digital_normal)) < 0.0:
                normal = -normal
            normal /= np.linalg.norm(normal)
            self.normals[key] = normal
            self.areas_m2[key] = self.voxel_size_m**2 / float(np.sum(np.abs(normal)))
            self.valid[key] = True
        self.recomputed_elements += len(keys & set(self.geometry))

    def update_after_insertions(
        self, mask: BoolArray, coordinates: Sequence[Coordinate], version: int
    ) -> int:
        changed: set[FaceKey] = set()
        for coordinate in coordinates:
            base = np.asarray(coordinate, dtype=int)
            for axis in range(3):
                lower = base.copy()
                lower[axis] -= 1
                changed.add((axis, *(int(v) for v in lower)))
                changed.add((axis, *(int(v) for v in base)))
        for key in changed:
            geometry = self._face_geometry(mask, key)
            if geometry is None:
                self.geometry.pop(key, None)
                self.normals.pop(key, None)
                self.areas_m2.pop(key, None)
                self.valid.pop(key, None)
            else:
                self.geometry[key] = geometry
        if not self.geometry:
            self.version = int(version)
            return 0
        ordered = sorted(self.geometry)
        points = np.asarray([self.geometry[key][0] for key in ordered])
        tree = cKDTree(points)
        origin = self._origin()
        insertion_centres = np.asarray([
            origin + (np.asarray(coordinate, dtype=float) + 0.5) * self.voxel_size_m
            for coordinate in coordinates
        ])
        radius = self.fitting_radius_m + math.sqrt(3.0) * self.voxel_size_m
        affected_indices: set[int] = set()
        for centre in insertion_centres:
            affected_indices.update(tree.query_ball_point(centre, radius))
        affected = {ordered[index] for index in affected_indices} | (changed & set(self.geometry))
        before = self.recomputed_elements
        self._recompute(affected)
        self.version = int(version)
        return self.recomputed_elements - before

    def to_surface_elements(self) -> SurfaceElements:
        ordered = sorted(self.geometry)
        return SurfaceElements(
            areas_m2=np.asarray([self.areas_m2[key] for key in ordered]),
            centroids_m=np.asarray([self.geometry[key][0] for key in ordered]),
            normals=np.asarray([self.normals[key] for key in ordered]),
            candidate_coordinates=np.asarray([self.geometry[key][1] for key in ordered], dtype=np.int64),
            normal_valid=np.asarray([self.valid[key] for key in ordered], dtype=bool),
            method="unique_faces_l1_corrected_local_plane_fit",
            runtime_s=0.0,
        )


@dataclass
class FrontRegion:
    """One evolving solid mask and its predeclared active front window."""

    name: str
    solid_mask: BoolArray
    reference_normal: FloatArray
    centre_m: FloatArray
    radius_m: float | None = None
    minimum_alignment: float = 0.5
    version: int = 0

    def __post_init__(self) -> None:
        self.solid_mask = np.asarray(self.solid_mask, dtype=bool).copy()
        if self.solid_mask.ndim != 3 or min(self.solid_mask.shape) < 3:
            raise ValueError("solid mask must be a three-dimensional array")
        if not np.any(self.solid_mask) or np.all(self.solid_mask):
            raise ValueError("solid mask must contain solid and empty cells")
        normal = np.asarray(self.reference_normal, dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            raise ValueError("reference normal must be a finite vector")
        magnitude = float(np.linalg.norm(normal))
        if magnitude <= 0.0:
            raise ValueError("reference normal must be nonzero")
        self.reference_normal = normal / magnitude
        self.centre_m = np.asarray(self.centre_m, dtype=float)
        if self.centre_m.shape != (3,) or not np.all(np.isfinite(self.centre_m)):
            raise ValueError("centre must be a finite three-vector")
        if self.radius_m is not None and (
            not math.isfinite(self.radius_m) or self.radius_m <= 0.0
        ):
            raise ValueError("radius must be finite and positive")
        if not 0.0 <= self.minimum_alignment <= 1.0:
            raise ValueError("minimum alignment must lie in [0, 1]")


@dataclass(frozen=True)
class FrontSnapshot:
    name: str
    reconstructed_area_m2: float
    area_weighted_position_m: float
    active_candidate_count: int
    active_fractional_volume_m3: float
    inactive_fractional_volume_m3: float
    realized_voxels: int
    cumulative_allocated_volume_m3: float


@dataclass(frozen=True)
class DynamicUpdate:
    transaction_id: int
    incoming_mass_kg: float
    incoming_volume_m3: float
    interface_allocated_volume_m3: float
    deferred_volume_m3: float
    realized_voxels: int
    event_batches: int
    reconstruction_count: int
    local_patch_update_count: int
    local_recomputed_elements: int
    active_fractional_volume_m3: float
    inactive_fractional_volume_m3: float
    closure_error_m3: float
    runtime_s: float


def _validate_rotation(rotation: ArrayLike) -> FloatArray:
    matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
        raise ValueError("grain_to_lab must be a finite 3 by 3 matrix")
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1.0e-12):
        raise ValueError("grain_to_lab must be orthogonal")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1.0e-12):
        raise ValueError("grain_to_lab must be a proper rotation")
    return matrix


def _subset_elements(elements: SurfaceElements, keep: NDArray[np.bool_]) -> SurfaceElements:
    return SurfaceElements(
        areas_m2=elements.areas_m2[keep],
        centroids_m=elements.centroids_m[keep],
        normals=elements.normals[keep],
        candidate_coordinates=elements.candidate_coordinates[keep],
        normal_valid=elements.normal_valid[keep],
        method=elements.method + "_front_window",
        runtime_s=elements.runtime_s,
    )


def candidate_face_flux_audit(
    elements: SurfaceElements,
    grain_to_lab: ArrayLike,
    alpha: float,
) -> list[dict[str, object]]:
    """Compare face-integrated and candidate-bundle ``A*kappa`` weights."""

    rotation = _validate_rotation(grain_to_lab)
    if not np.all(elements.normal_valid):
        raise ValueError("face-flux audit requires valid reconstructed normals")
    face_kappa = cubic_attachment_weight(elements.normals, rotation, alpha)
    face_weights = elements.areas_m2 * face_kappa
    groups: dict[Coordinate, list[int]] = {}
    for index, coordinate in enumerate(elements.candidate_coordinates):
        groups.setdefault(tuple(int(v) for v in coordinate), []).append(index)
    bundles = {patch.coordinate: patch for patch in bundle_elements_by_candidate(elements)}
    rows = []
    for coordinate in sorted(groups):
        indices = np.asarray(groups[coordinate], dtype=int)
        face_integrated = math.fsum(float(face_weights[index]) for index in indices)
        bundle = bundles[coordinate]
        if not bundle.normal_valid:
            raise RuntimeError("candidate bundle has invalid normal")
        bundle_kappa = float(
            cubic_attachment_weight(np.asarray(bundle.normal)[None, :], rotation, alpha)[0]
        )
        adapter = bundle.area_m2 * bundle_kappa
        rows.append({
            "coordinate": coordinate,
            "element_count": len(indices),
            "face_integrated_weight_m2": face_integrated,
            "adapter_weight_m2": adapter,
            "signed_relative_discrepancy": adapter / face_integrated - 1.0,
            "absolute_relative_discrepancy": abs(adapter / face_integrated - 1.0),
        })
    return rows


def perpendicular_face_counterexample(alpha: float) -> dict[str, float]:
    """Exact equal-area +x/+y counterexample to normal-before-kappa bundling."""

    face_normal = np.asarray([[1.0, 0.0, 0.0]])
    bundle_normal = np.asarray([[1.0, 1.0, 0.0]]) / math.sqrt(2.0)
    face_kappa = float(cubic_attachment_weight(face_normal, np.eye(3), alpha)[0])
    bundle_kappa = float(cubic_attachment_weight(bundle_normal, np.eye(3), alpha)[0])
    face_integrated = 2.0 * face_kappa
    adapter = 2.0 * bundle_kappa
    return {
        "alpha": float(alpha),
        "face_integrated_weight_per_unit_area": face_integrated,
        "adapter_weight_per_unit_area": adapter,
        "adapter_to_face_ratio": adapter / face_integrated,
        "exact_ratio_exp_minus_alpha_over_2": math.exp(-0.5 * alpha),
    }


class DynamicInterfaceFluxSystem:
    """Event-driven evolving fronts with persistent candidate-owned volume."""

    def __init__(
        self,
        fronts: Sequence[FrontRegion],
        *,
        voxel_size_m: float,
        fitting_radius_m: float,
        solid_density_kg_m3: float,
        grain_to_lab: ArrayLike,
        alpha: float,
        rng_seed: int = 0,
    ) -> None:
        if not fronts or len({front.name for front in fronts}) != len(fronts):
            raise ValueError("front names must be non-empty and unique")
        if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
            raise ValueError("voxel size must be finite and positive")
        if not math.isfinite(fitting_radius_m) or fitting_radius_m <= voxel_size_m:
            raise ValueError("fitting radius must exceed one voxel")
        if not math.isfinite(solid_density_kg_m3) or solid_density_kg_m3 <= 0.0:
            raise ValueError("solid density must be finite and positive")
        if not math.isfinite(alpha):
            raise ValueError("alpha must be finite")
        shapes = {front.solid_mask.shape for front in fronts}
        if len(shapes) != 1:
            raise ValueError("all fronts must use the same grid shape")
        self.fronts = {front.name: copy.deepcopy(front) for front in fronts}
        self.voxel_size_m = float(voxel_size_m)
        self.fitting_radius_m = float(fitting_radius_m)
        self.solid_density_kg_m3 = float(solid_density_kg_m3)
        self.voxel_volume_m3 = self.voxel_size_m**3
        self.voxel_mass_kg = self.solid_density_kg_m3 * self.voxel_volume_m3
        self.grain_to_lab = _validate_rotation(grain_to_lab).copy()
        self.alpha = float(alpha)
        self.rng = np.random.default_rng(rng_seed)
        self.local_reservoirs_m3: dict[Owner, float] = {}
        self.cumulative_local_allocated_m3: dict[Owner, float] = {}
        self.cumulative_input_volume_m3 = 0.0
        self.cumulative_interface_allocated_volume_m3 = 0.0
        self.cumulative_deferred_volume_m3 = 0.0
        self.cumulative_realized_voxels = 0
        self.realized_by_front = {name: 0 for name in self.fronts}
        self.allocated_by_front_m3 = {name: 0.0 for name in self.fronts}
        self.transaction_id = 0
        self._surface_cache: dict[str, LocalCorrectedInterfaceCache] = {}
        self.reconstruction_count = 0
        self.local_patch_update_count = 0
        self.local_recomputed_elements = 0
        self.history: list[DynamicUpdate] = []
        self._initial_positions_m = {
            name: self._front_geometry(name)[1] for name in self.fronts
        }
        self.assert_invariants()

    def snapshot_state(self) -> object:
        return copy.deepcopy({
            "fronts": self.fronts,
            "local_reservoirs_m3": self.local_reservoirs_m3,
            "cumulative_local_allocated_m3": self.cumulative_local_allocated_m3,
            "cumulative_input_volume_m3": self.cumulative_input_volume_m3,
            "cumulative_interface_allocated_volume_m3": self.cumulative_interface_allocated_volume_m3,
            "cumulative_deferred_volume_m3": self.cumulative_deferred_volume_m3,
            "cumulative_realized_voxels": self.cumulative_realized_voxels,
            "realized_by_front": self.realized_by_front,
            "allocated_by_front_m3": self.allocated_by_front_m3,
            "transaction_id": self.transaction_id,
            "surface_cache": self._surface_cache,
            "reconstruction_count": self.reconstruction_count,
            "local_patch_update_count": self.local_patch_update_count,
            "local_recomputed_elements": self.local_recomputed_elements,
            "history": self.history,
            "rng_state": self.rng.bit_generator.state,
            "initial_positions_m": self._initial_positions_m,
        })

    def restore_state(self, snapshot: object) -> None:
        if not isinstance(snapshot, dict):
            raise TypeError("invalid dynamic interface snapshot")
        restored = copy.deepcopy(snapshot)
        self.fronts = restored["fronts"]
        self.local_reservoirs_m3 = restored["local_reservoirs_m3"]
        self.cumulative_local_allocated_m3 = restored["cumulative_local_allocated_m3"]
        self.cumulative_input_volume_m3 = restored["cumulative_input_volume_m3"]
        self.cumulative_interface_allocated_volume_m3 = restored[
            "cumulative_interface_allocated_volume_m3"
        ]
        self.cumulative_deferred_volume_m3 = restored["cumulative_deferred_volume_m3"]
        self.cumulative_realized_voxels = restored["cumulative_realized_voxels"]
        self.realized_by_front = restored["realized_by_front"]
        self.allocated_by_front_m3 = restored["allocated_by_front_m3"]
        self.transaction_id = restored["transaction_id"]
        self._surface_cache = restored["surface_cache"]
        self.reconstruction_count = restored["reconstruction_count"]
        self.local_patch_update_count = restored["local_patch_update_count"]
        self.local_recomputed_elements = restored["local_recomputed_elements"]
        self.history = restored["history"]
        self.rng.bit_generator.state = restored["rng_state"]
        self._initial_positions_m = restored["initial_positions_m"]

    def _full_surface(self, name: str) -> SurfaceElements:
        front = self.fronts[name]
        cached = self._surface_cache.get(name)
        if cached is not None and cached.version == front.version:
            return cached.to_surface_elements()
        cache = LocalCorrectedInterfaceCache(
            front.solid_mask, self.voxel_size_m, self.fitting_radius_m, front.version
        )
        self.reconstruction_count += 1
        self._surface_cache[name] = cache
        return cache.to_surface_elements()

    def _active_surface(self, name: str) -> SurfaceElements:
        front = self.fronts[name]
        elements = self._full_surface(name)
        alignment = elements.normals @ front.reference_normal
        offsets = elements.centroids_m - front.centre_m
        axial = offsets @ front.reference_normal
        transverse = offsets - axial[:, None] * front.reference_normal
        keep = elements.normal_valid & (alignment >= front.minimum_alignment)
        if front.radius_m is not None:
            keep &= np.linalg.norm(transverse, axis=1) <= front.radius_m
        if not np.any(keep):
            return SurfaceElements(
                areas_m2=np.empty(0), centroids_m=np.empty((0, 3)),
                normals=np.empty((0, 3)), candidate_coordinates=np.empty((0, 3), dtype=np.int64),
                normal_valid=np.empty(0, dtype=bool), method=elements.method + "_empty_window",
                runtime_s=elements.runtime_s,
            )
        return _subset_elements(elements, keep)

    def _candidate_fluxes(self) -> tuple[dict[Owner, float], dict[str, float]]:
        weights: dict[Owner, float] = {}
        areas = {name: 0.0 for name in self.fronts}
        for name in self.fronts:
            elements = self._active_surface(name)
            if not len(elements.areas_m2):
                continue
            kappa = cubic_attachment_weight(elements.normals, self.grain_to_lab, self.alpha)
            face_weights = elements.areas_m2 * kappa
            areas[name] = math.fsum(float(v) for v in elements.areas_m2)
            for coordinate, weight in zip(elements.candidate_coordinates, face_weights, strict=True):
                owner = (name, tuple(int(v) for v in coordinate))
                weights[owner] = weights.get(owner, 0.0) + float(weight)
        if any(not math.isfinite(value) or value <= 0.0 for value in weights.values()):
            raise RuntimeError("candidate flux weights must be finite and positive")
        return weights, areas

    def _front_geometry(self, name: str) -> tuple[float, float, int]:
        elements = self._active_surface(name)
        if not len(elements.areas_m2):
            return 0.0, float("nan"), 0
        area = math.fsum(float(v) for v in elements.areas_m2)
        position = float(
            np.sum(elements.areas_m2 * (elements.centroids_m @ self.fronts[name].reference_normal))
            / area
        )
        candidates = len({tuple(int(v) for v in row) for row in elements.candidate_coordinates})
        return area, position, candidates

    def active_owners(self) -> set[Owner]:
        weights, _ = self._candidate_fluxes()
        return set(weights)

    def closure_error_m3(self) -> float:
        return self.cumulative_input_volume_m3 - (
            self.cumulative_realized_voxels * self.voxel_volume_m3
            + math.fsum(self.local_reservoirs_m3.values())
            + self.cumulative_deferred_volume_m3
        )

    def assert_invariants(self) -> None:
        if any(
            not math.isfinite(value) or value < -64.0 * math.ulp(self.voxel_volume_m3)
            or value >= self.voxel_volume_m3 * (1.0 + 1.0e-10)
            for value in self.local_reservoirs_m3.values()
        ):
            raise RuntimeError("local reservoir is outside [0, voxel volume)")
        realized = sum(self.realized_by_front.values())
        if realized != self.cumulative_realized_voxels:
            raise RuntimeError("per-front and global realized inventories disagree")
        scale = max(self.cumulative_input_volume_m3, self.voxel_volume_m3)
        if abs(self.closure_error_m3()) > 256.0 * math.ulp(scale):
            raise RuntimeError("dynamic volume ledger does not close")
        if not math.isclose(
            self.cumulative_interface_allocated_volume_m3
            + self.cumulative_deferred_volume_m3,
            self.cumulative_input_volume_m3,
            rel_tol=0.0,
            abs_tol=256.0 * math.ulp(scale),
        ):
            raise RuntimeError("allocated plus deferred volume differs from input")

    def front_snapshot(self, name: str) -> FrontSnapshot:
        area, position, candidates = self._front_geometry(name)
        active = self.active_owners()
        active_volume = math.fsum(
            volume for owner, volume in self.local_reservoirs_m3.items()
            if owner[0] == name and owner in active
        )
        inactive_volume = math.fsum(
            volume for owner, volume in self.local_reservoirs_m3.items()
            if owner[0] == name and owner not in active
        )
        return FrontSnapshot(
            name=name,
            reconstructed_area_m2=area,
            area_weighted_position_m=position,
            active_candidate_count=candidates,
            active_fractional_volume_m3=active_volume,
            inactive_fractional_volume_m3=inactive_volume,
            realized_voxels=self.realized_by_front[name],
            cumulative_allocated_volume_m3=self.allocated_by_front_m3[name],
        )

    def local_cache_consistency(self, name: str) -> dict[str, float | int | bool]:
        """Compare current local patch data with a fresh full reconstruction."""

        self._full_surface(name)
        local = self._surface_cache[name]
        front = self.fronts[name]
        full = LocalCorrectedInterfaceCache(
            front.solid_mask,
            self.voxel_size_m,
            self.fitting_radius_m,
            front.version,
        )
        keys_exact = set(local.geometry) == set(full.geometry)
        if not keys_exact:
            return {
                "keys_exact": False,
                "local_face_count": len(local.geometry),
                "full_face_count": len(full.geometry),
                "maximum_area_abs_error_m2": float("inf"),
                "maximum_normal_vector_error": float("inf"),
            }
        area_error = max(
            (abs(local.areas_m2[key] - full.areas_m2[key]) for key in full.geometry),
            default=0.0,
        )
        normal_error = max(
            (
                float(np.linalg.norm(local.normals[key] - full.normals[key]))
                for key in full.geometry
                if local.valid[key] and full.valid[key]
            ),
            default=0.0,
        )
        valid_exact = all(local.valid[key] == full.valid[key] for key in full.geometry)
        return {
            "keys_exact": True,
            "valid_flags_exact": valid_exact,
            "local_face_count": len(local.geometry),
            "full_face_count": len(full.geometry),
            "maximum_area_abs_error_m2": area_error,
            "maximum_normal_vector_error": normal_error,
        }

    def advance(
        self,
        incoming_mass_kg: float,
        *,
        inject_failure_after_event_batches: int | None = None,
        maximum_event_batches: int = 100_000,
    ) -> DynamicUpdate:
        """Allocate and realize one input transaction with event substepping."""

        if not math.isfinite(incoming_mass_kg) or incoming_mass_kg < 0.0:
            raise ValueError("incoming mass must be finite and nonnegative")
        before = self.snapshot_state()
        started = perf_counter()
        reconstructions_before = self.reconstruction_count
        local_updates_before = self.local_patch_update_count
        local_elements_before = self.local_recomputed_elements
        realized_before = self.cumulative_realized_voxels
        allocated_before = self.cumulative_interface_allocated_volume_m3
        deferred_before = self.cumulative_deferred_volume_m3
        incoming_volume = incoming_mass_kg / self.solid_density_kg_m3
        remaining = incoming_volume
        event_batches = 0
        try:
            self.transaction_id += 1
            self.cumulative_input_volume_m3 += incoming_volume
            tolerance = 128.0 * math.ulp(self.voxel_volume_m3)
            while remaining > 128.0 * math.ulp(max(incoming_volume, self.voxel_volume_m3)):
                weights, _ = self._candidate_fluxes()
                if not weights:
                    self.cumulative_deferred_volume_m3 += remaining
                    remaining = 0.0
                    break
                total_weight = math.fsum(weights.values())
                owners = sorted(weights)
                shares = np.asarray([weights[owner] / total_weight for owner in owners])
                gaps = np.asarray([
                    self.voxel_volume_m3 - self.local_reservoirs_m3.get(owner, 0.0)
                    for owner in owners
                ])
                if np.any(gaps <= tolerance):
                    critical_indices = np.flatnonzero(gaps <= tolerance)
                    allocation = 0.0
                else:
                    thresholds = gaps / shares
                    allocation = min(remaining, float(np.min(thresholds)))
                    additions = allocation * shares
                    residual = allocation - math.fsum(float(v) for v in additions)
                    additions[int(np.argmax(shares))] += residual
                    for owner, addition in zip(owners, additions, strict=True):
                        self.local_reservoirs_m3[owner] = (
                            self.local_reservoirs_m3.get(owner, 0.0) + float(addition)
                        )
                        self.cumulative_local_allocated_m3[owner] = (
                            self.cumulative_local_allocated_m3.get(owner, 0.0) + float(addition)
                        )
                        self.allocated_by_front_m3[owner[0]] += float(addition)
                    self.cumulative_interface_allocated_volume_m3 += allocation
                    remaining -= allocation
                    critical_indices = np.flatnonzero(
                        np.asarray([self.local_reservoirs_m3[owner] for owner in owners])
                        >= self.voxel_volume_m3 - tolerance
                    )
                if critical_indices.size:
                    event_batches += 1
                    if event_batches > maximum_event_batches:
                        raise RuntimeError("maximum dynamic event batches exceeded")
                    changed_fronts: dict[str, list[Coordinate]] = {}
                    for index in critical_indices:
                        owner = owners[int(index)]
                        name, coordinate = owner
                        front = self.fronts[name]
                        if front.solid_mask[coordinate]:
                            raise RuntimeError("candidate voxel was already solid")
                        reservoir = self.local_reservoirs_m3[owner]
                        excess = reservoir - self.voxel_volume_m3
                        if abs(excess) > tolerance:
                            raise RuntimeError("event substep overshot a local voxel reservoir")
                        self.local_reservoirs_m3[owner] = max(0.0, excess)
                        front.solid_mask[coordinate] = True
                        front.version += 1
                        changed_fronts.setdefault(name, []).append(coordinate)
                        self.cumulative_realized_voxels += 1
                        self.realized_by_front[name] += 1
                    for name, coordinates in changed_fronts.items():
                        cache = self._surface_cache.get(name)
                        if cache is None:
                            continue
                        recomputed = cache.update_after_insertions(
                            self.fronts[name].solid_mask,
                            coordinates,
                            self.fronts[name].version,
                        )
                        self.local_patch_update_count += 1
                        self.local_recomputed_elements += recomputed
                    if inject_failure_after_event_batches == event_batches:
                        raise RuntimeError("injected dynamic interface failure")
                elif allocation >= remaining + allocation - tolerance:
                    # No event occurred and the remaining input was fully stored.
                    remaining = 0.0
                    break
            self.assert_invariants()
            active = self.active_owners()
            active_fractional = math.fsum(
                volume for owner, volume in self.local_reservoirs_m3.items() if owner in active
            )
            inactive_fractional = math.fsum(
                volume for owner, volume in self.local_reservoirs_m3.items() if owner not in active
            )
            result = DynamicUpdate(
                transaction_id=self.transaction_id,
                incoming_mass_kg=float(incoming_mass_kg),
                incoming_volume_m3=incoming_volume,
                interface_allocated_volume_m3=(
                    self.cumulative_interface_allocated_volume_m3 - allocated_before
                ),
                deferred_volume_m3=self.cumulative_deferred_volume_m3 - deferred_before,
                realized_voxels=self.cumulative_realized_voxels - realized_before,
                event_batches=event_batches,
                reconstruction_count=self.reconstruction_count - reconstructions_before,
                local_patch_update_count=self.local_patch_update_count - local_updates_before,
                local_recomputed_elements=self.local_recomputed_elements - local_elements_before,
                active_fractional_volume_m3=active_fractional,
                inactive_fractional_volume_m3=inactive_fractional,
                closure_error_m3=self.closure_error_m3(),
                runtime_s=perf_counter() - started,
            )
            self.history.append(result)
            return result
        except Exception:
            self.restore_state(before)
            self.assert_invariants()
            raise
