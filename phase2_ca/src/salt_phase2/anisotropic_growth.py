"""Candidate-normal mapping and experimental cubic weighting for Phase 2A.3b-0.

This module deliberately separates reconstructed geometric normals from
crystallographic interpretation.  Edge/corner labels describe the digital
candidate geometry only; they are not facet labels.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import math
from time import perf_counter
from typing import Literal, Sequence

import numpy as np
from numpy.typing import NDArray

from .state import Coordinate
from .surface_orientation import (
    cubic_attachment_weight,
    estimate_plane_fit_normals,
)


FloatArray = NDArray[np.float64]
EdgePolicy = Literal["smooth_average", "isotropic_ambiguous"]
UpdateStrategy = Literal["full", "local", "cached"]

_FACE_OFFSETS: tuple[Coordinate, ...] = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)


@dataclass(frozen=True)
class CandidateNormalDiagnostic:
    """Auditable geometric and anisotropy information for one candidate."""

    coordinate: Coordinate
    site_class: str
    adjacent_salt_faces: int
    valid_face_normals: int
    normal_set: tuple[tuple[float, float, float], ...]
    combined_normal: tuple[float, float, float]
    normal_quality: float
    anisotropy_eligible: bool
    anisotropy_weight: float
    edge_policy: str


def _validate_rotation(rotation: NDArray[np.float64]) -> FloatArray:
    matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("grain_to_lab must have shape (3, 3)")
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1.0e-12):
        raise ValueError("grain_to_lab must be orthogonal")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1.0e-12):
        raise ValueError("grain_to_lab must be a proper rotation")
    return matrix


def _site_class(offsets: Sequence[Coordinate]) -> str:
    if len(offsets) == 1:
        return "face_like"
    if len(offsets) == 2 and int(np.dot(offsets[0], offsets[1])) == 0:
        return "edge_like"
    return "corner_like"


def candidate_interface_faces(
    salt_mask: NDArray[np.bool_],
    coordinates: Sequence[Coordinate],
    voxel_size_m: float,
) -> tuple[FloatArray, list[tuple[int, int]], list[str]]:
    """Map empty candidate voxels to adjacent salt-interface face centres.

    Returns one point per distinct candidate/salt face, slices into that point
    array for each candidate, and digital site classes.
    """

    salt = np.asarray(salt_mask, dtype=bool)
    if salt.ndim != 3 or min(salt.shape) < 3:
        raise ValueError("salt_mask must be a three-dimensional array of size >= 3")
    if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
        raise ValueError("voxel_size_m must be positive and finite")
    origin = -0.5 * np.asarray(salt.shape, dtype=float) * voxel_size_m
    points: list[FloatArray] = []
    slices: list[tuple[int, int]] = []
    classes: list[str] = []
    for coordinate in coordinates:
        if salt[coordinate]:
            raise ValueError("candidate coordinate is already salt")
        offsets: list[Coordinate] = []
        centre = origin + (np.asarray(coordinate, dtype=float) + 0.5) * voxel_size_m
        start = len(points)
        for offset in _FACE_OFFSETS:
            neighbour = tuple(coordinate[i] + offset[i] for i in range(3))
            if not all(0 <= neighbour[i] < salt.shape[i] for i in range(3)):
                continue
            if salt[neighbour]:
                offsets.append(offset)
                points.append(centre + 0.5 * np.asarray(offset) * voxel_size_m)
        if not offsets:
            raise ValueError("candidate has no adjacent salt face")
        slices.append((start, len(points)))
        classes.append(_site_class(offsets))
    return np.asarray(points, dtype=float), slices, classes


def map_candidate_normals(
    salt_mask: NDArray[np.bool_],
    coordinates: Sequence[Coordinate],
    *,
    voxel_size_m: float,
    fitting_radius_m: float,
    grain_to_lab: NDArray[np.float64],
    alpha: float,
    edge_policy: EdgePolicy,
    quality_minimum: float = 0.75,
) -> tuple[CandidateNormalDiagnostic, ...]:
    """Reconstruct candidate normals and apply one predeclared edge policy."""

    if edge_policy not in {"smooth_average", "isotropic_ambiguous"}:
        raise ValueError("unknown edge_policy")
    if not math.isfinite(alpha):
        raise ValueError("alpha must be finite")
    if not 0.0 <= quality_minimum <= 1.0:
        raise ValueError("quality_minimum must lie in [0, 1]")
    rotation = _validate_rotation(grain_to_lab)
    checked = [tuple(int(v) for v in coordinate) for coordinate in coordinates]
    if not checked:
        return ()
    points, slices, classes = candidate_interface_faces(
        salt_mask, checked, voxel_size_m
    )
    estimate = estimate_plane_fit_normals(
        salt_mask,
        voxel_size_m,
        fitting_radius_m,
        points_m=points,
    )
    records: list[CandidateNormalDiagnostic] = []
    for coordinate, (start, stop), site_class in zip(
        checked, slices, classes, strict=True
    ):
        local_valid = estimate.valid[start:stop]
        local = estimate.normals[start:stop][local_valid]
        normal_set = tuple(tuple(float(value) for value in row) for row in local)
        if local.shape[0]:
            resultant = np.mean(local, axis=0)
            resultant_length = float(np.linalg.norm(resultant))
            if resultant_length > 1.0e-14:
                combined = resultant / resultant_length
            else:
                combined = np.full(3, np.nan)
            valid_fraction = local.shape[0] / (stop - start)
            quality = float(valid_fraction * resultant_length)
        else:
            combined = np.full(3, np.nan)
            quality = 0.0
        eligible = bool(
            np.all(np.isfinite(combined))
            and quality >= quality_minimum
            and (
                edge_policy == "smooth_average"
                or site_class == "face_like"
            )
        )
        weight = (
            float(cubic_attachment_weight(combined[None, :], rotation, alpha)[0])
            if eligible
            else 1.0
        )
        records.append(
            CandidateNormalDiagnostic(
                coordinate=coordinate,
                site_class=site_class,
                adjacent_salt_faces=stop - start,
                valid_face_normals=int(local.shape[0]),
                normal_set=normal_set,
                combined_normal=tuple(float(value) for value in combined),
                normal_quality=quality,
                anisotropy_eligible=eligible,
                anisotropy_weight=weight,
                edge_policy=edge_policy,
            )
        )
    return tuple(records)


class OrientationCandidateWeightModel:
    """Stateful, rollback-safe candidate weighting with explicit invalidation."""

    def __init__(
        self,
        *,
        voxel_size_m: float,
        fitting_radius_m: float,
        grain_to_lab: NDArray[np.float64],
        alpha: float,
        edge_policy: EdgePolicy,
        update_strategy: UpdateStrategy = "full",
        quality_minimum: float = 0.75,
    ) -> None:
        if update_strategy not in {"full", "local", "cached"}:
            raise ValueError("unknown update_strategy")
        self.voxel_size_m = float(voxel_size_m)
        self.fitting_radius_m = float(fitting_radius_m)
        self.grain_to_lab = _validate_rotation(grain_to_lab).copy()
        self.alpha = float(alpha)
        self.edge_policy = edge_policy
        self.update_strategy = update_strategy
        self.quality_minimum = float(quality_minimum)
        self._cache: dict[Coordinate, CandidateNormalDiagnostic] = {}
        self._last_records: dict[Coordinate, CandidateNormalDiagnostic] = {}
        self.selection_history: list[CandidateNormalDiagnostic] = []
        self.weight_calls = 0
        self.recomputed_candidates = 0
        self.invalidated_candidates = 0
        self.normal_runtime_s = 0.0
        self.candidate_evaluations = 0
        self.class_counts = {"face_like": 0, "edge_like": 0, "corner_like": 0}
        self.anisotropy_eligible_evaluations = 0

    def _compute(
        self,
        salt_mask: NDArray[np.bool_],
        coordinates: list[Coordinate],
    ) -> tuple[CandidateNormalDiagnostic, ...]:
        if not coordinates:
            return ()
        started = perf_counter()
        records = map_candidate_normals(
            salt_mask,
            coordinates,
            voxel_size_m=self.voxel_size_m,
            fitting_radius_m=self.fitting_radius_m,
            grain_to_lab=self.grain_to_lab,
            alpha=self.alpha,
            edge_policy=self.edge_policy,
            quality_minimum=self.quality_minimum,
        )
        self.normal_runtime_s += perf_counter() - started
        self.recomputed_candidates += len(coordinates)
        return records

    def weights_for_candidates(
        self,
        salt_mask: NDArray[np.bool_],
        coordinates: list[Coordinate],
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
        for record in records:
            self.class_counts[record.site_class] += 1
            self.anisotropy_eligible_evaluations += int(record.anisotropy_eligible)
        weights = np.asarray([record.anisotropy_weight for record in records], dtype=float)
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise RuntimeError("candidate anisotropy weights must be finite and positive")
        return weights

    def record_selection(self, coordinate: Coordinate) -> None:
        try:
            self.selection_history.append(self._last_records[coordinate])
        except KeyError as error:
            raise RuntimeError("selected candidate lacks a normal diagnostic") from error

    def diagnostics_for_candidates(
        self,
        salt_mask: NDArray[np.bool_],
        coordinates: list[Coordinate],
    ) -> tuple[CandidateNormalDiagnostic, ...]:
        """Return current ordered diagnostics through the configured cache policy."""

        self.weights_for_candidates(salt_mask, coordinates)
        return tuple(self._last_records[coordinate] for coordinate in coordinates)

    def notify_insertion(self, coordinate: Coordinate) -> None:
        if self.update_strategy == "full":
            invalid = set(self._cache)
        elif self.update_strategy == "local":
            radius_cells = self.fitting_radius_m / self.voxel_size_m + math.sqrt(3.0)
            centre = np.asarray(coordinate, dtype=float)
            invalid = {
                candidate
                for candidate in self._cache
                if np.linalg.norm(np.asarray(candidate, dtype=float) - centre)
                <= radius_cells
            }
        else:
            invalid = {coordinate}
            invalid.update(
                tuple(coordinate[i] + offset[i] for i in range(3))
                for offset in _FACE_OFFSETS
            )
        for candidate in invalid:
            self._cache.pop(candidate, None)
        self.invalidated_candidates += len(invalid)

    def snapshot_state(self) -> object:
        return copy.deepcopy(
            {
                "cache": self._cache,
                "last_records": self._last_records,
                "selection_history": self.selection_history,
                "weight_calls": self.weight_calls,
                "recomputed_candidates": self.recomputed_candidates,
                "invalidated_candidates": self.invalidated_candidates,
                "normal_runtime_s": self.normal_runtime_s,
                "candidate_evaluations": self.candidate_evaluations,
                "class_counts": self.class_counts,
                "anisotropy_eligible_evaluations": self.anisotropy_eligible_evaluations,
            }
        )

    def restore_state(self, snapshot: object) -> None:
        if not isinstance(snapshot, dict):
            raise TypeError("invalid candidate-weight snapshot")
        restored = copy.deepcopy(snapshot)
        self._cache = restored["cache"]
        self._last_records = restored["last_records"]
        self.selection_history = restored["selection_history"]
        self.weight_calls = restored["weight_calls"]
        self.recomputed_candidates = restored["recomputed_candidates"]
        self.invalidated_candidates = restored["invalidated_candidates"]
        self.normal_runtime_s = restored["normal_runtime_s"]
        self.candidate_evaluations = restored["candidate_evaluations"]
        self.class_counts = restored["class_counts"]
        self.anisotropy_eligible_evaluations = restored[
            "anisotropy_eligible_evaluations"
        ]

    def performance_summary(self) -> dict[str, float | int | str]:
        return {
            "update_strategy": self.update_strategy,
            "weight_calls": self.weight_calls,
            "recomputed_candidates": self.recomputed_candidates,
            "invalidated_candidates": self.invalidated_candidates,
            "normal_runtime_s": self.normal_runtime_s,
            "selected_events": len(self.selection_history),
            "candidate_evaluations": self.candidate_evaluations,
            "face_like_evaluations": self.class_counts["face_like"],
            "edge_like_evaluations": self.class_counts["edge_like"],
            "corner_like_evaluations": self.class_counts["corner_like"],
            "anisotropy_eligible_evaluations": self.anisotropy_eligible_evaluations,
        }
