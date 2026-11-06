"""Growth-only diagnostics for the Phase 2A.3b-1 mechanism verification."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .anisotropic_growth import CandidateNormalDiagnostic, map_candidate_normals
from .state import Coordinate
from .surface_orientation import cubic_q4


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class CrystallographicCandidate:
    coordinate: Coordinate
    site_class: str
    adjacent_salt_faces: int
    normal_lab: tuple[float, float, float]
    normal_crystal: tuple[float, float, float]
    q4: float
    orientation_family: str
    family_angular_error_deg: float
    represented_area_m2: float
    normal_quality: float
    anisotropy_eligible: bool
    anisotropy_weight: float


def _unit_rows(values: NDArray[np.float64]) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    norms = np.linalg.norm(array, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0.0):
        raise ValueError("direction rows must be finite and nonzero")
    return array / norms[:, None]


def ideal_cubic_direction_families() -> dict[str, FloatArray]:
    axes = []
    face_diagonals = []
    body_diagonals = []
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for z in (-1, 0, 1):
                nonzero = sum(value != 0 for value in (x, y, z))
                if nonzero == 1:
                    axes.append((x, y, z))
                elif nonzero == 2:
                    face_diagonals.append((x, y, z))
                elif nonzero == 3:
                    body_diagonals.append((x, y, z))
    return {
        "<100>": _unit_rows(np.asarray(axes, dtype=float)),
        "<110>": _unit_rows(np.asarray(face_diagonals, dtype=float)),
        "<111>": _unit_rows(np.asarray(body_diagonals, dtype=float)),
    }


def classify_crystal_directions(
    normals_crystal: NDArray[np.float64],
    *,
    near_tolerance_deg: float = 15.0,
) -> tuple[NDArray[np.str_], FloatArray]:
    """Classify by the nearest ideal cubic direction, retaining an `other` bin."""

    values = np.asarray(normals_crystal, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("normals_crystal must have shape (n, 3)")
    if not 0.0 < near_tolerance_deg < 45.0:
        raise ValueError("near_tolerance_deg must lie in (0, 45)")
    norms = np.linalg.norm(values, axis=1)
    valid = np.all(np.isfinite(values), axis=1) & np.isfinite(norms) & (norms > 0.0)
    labels = np.full(values.shape[0], "other", dtype="<U7")
    minimum = np.full(values.shape[0], np.nan, dtype=float)
    if not np.any(valid):
        return labels, minimum
    normals = values[valid] / norms[valid, None]
    families = ideal_cubic_direction_families()
    names = list(families)
    angles = np.column_stack(
        [
            np.degrees(
                np.arccos(
                    np.clip(np.max(normals @ directions.T, axis=1), -1.0, 1.0)
                )
            )
            for directions in families.values()
        ]
    )
    nearest = np.argmin(angles, axis=1)
    valid_minimum = angles[np.arange(normals.shape[0]), nearest]
    valid_labels = np.asarray([names[index] for index in nearest], dtype="<U7")
    valid_labels[valid_minimum > near_tolerance_deg] = "other"
    labels[valid] = valid_labels
    minimum[valid] = valid_minimum
    return labels, minimum


def crystallographic_candidate_diagnostics(
    salt_mask: NDArray[np.bool_],
    coordinates: Sequence[Coordinate],
    *,
    voxel_size_m: float,
    fitting_radius_m: float,
    grain_to_lab: NDArray[np.float64],
    alpha: float,
    edge_policy: str,
    quality_minimum: float,
    near_tolerance_deg: float = 15.0,
) -> tuple[CrystallographicCandidate, ...]:
    geometric = map_candidate_normals(
        salt_mask,
        coordinates,
        voxel_size_m=voxel_size_m,
        fitting_radius_m=fitting_radius_m,
        grain_to_lab=grain_to_lab,
        alpha=alpha,
        edge_policy=edge_policy,
        quality_minimum=quality_minimum,
    )
    lab = np.asarray([record.combined_normal for record in geometric], dtype=float)
    rotation = np.asarray(grain_to_lab, dtype=float)
    crystal = lab @ rotation
    labels, errors = classify_crystal_directions(
        crystal, near_tolerance_deg=near_tolerance_deg
    )
    q4_values = cubic_q4(crystal)
    return tuple(
        CrystallographicCandidate(
            coordinate=record.coordinate,
            site_class=record.site_class,
            adjacent_salt_faces=record.adjacent_salt_faces,
            normal_lab=record.combined_normal,
            normal_crystal=tuple(float(value) for value in crystal[index]),
            q4=float(q4_values[index]),
            orientation_family=str(labels[index]),
            family_angular_error_deg=float(errors[index]),
            represented_area_m2=record.adjacent_salt_faces * voxel_size_m**2,
            normal_quality=record.normal_quality,
            anisotropy_eligible=record.anisotropy_eligible,
            anisotropy_weight=record.anisotropy_weight,
        )
        for index, record in enumerate(geometric)
    )


def fixed_front_probabilities(
    records: Sequence[CrystallographicCandidate],
) -> FloatArray:
    weights = np.asarray([record.anisotropy_weight for record in records], dtype=float)
    if weights.size == 0 or np.any(~np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("records must contain finite positive weights")
    return weights / np.sum(weights)


def sample_fixed_front(
    probabilities: NDArray[np.float64],
    *,
    sample_count: int,
    rng_seed: int,
) -> NDArray[np.int64]:
    if isinstance(sample_count, bool) or int(sample_count) != sample_count or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")
    probabilities = np.asarray(probabilities, dtype=float)
    if not np.isclose(np.sum(probabilities), 1.0, atol=1.0e-12):
        raise ValueError("probabilities must sum to one")
    return np.random.default_rng(int(rng_seed)).multinomial(
        int(sample_count), probabilities
    )


def air_exposed_face_area_m2(
    salt_mask: NDArray[np.bool_],
    empty_mask: NDArray[np.bool_],
    voxel_size_m: float,
) -> float:
    """Return internal salt/empty face area; substrate and domain exterior are excluded."""

    salt = np.asarray(salt_mask, dtype=bool)
    empty = np.asarray(empty_mask, dtype=bool)
    if salt.shape != empty.shape or salt.ndim != 3:
        raise ValueError("salt and empty masks must be matching 3D arrays")
    faces = 0
    for axis in range(3):
        left = [slice(None)] * 3
        right = [slice(None)] * 3
        left[axis] = slice(0, -1)
        right[axis] = slice(1, None)
        faces += int(np.count_nonzero(salt[tuple(left)] & empty[tuple(right)]))
        faces += int(np.count_nonzero(empty[tuple(left)] & salt[tuple(right)]))
    return faces * voxel_size_m**2


def directional_growth_extents(
    initial_mask: NDArray[np.bool_],
    final_mask: NDArray[np.bool_],
    *,
    voxel_size_m: float,
    grain_to_lab: NDArray[np.float64],
) -> dict[str, float]:
    if initial_mask.shape != final_mask.shape:
        raise ValueError("initial and final masks must match")
    if np.any(initial_mask & ~final_mask):
        raise ValueError("growth metrics require a persistent initial seed")
    initial = np.argwhere(initial_mask).astype(float) * voxel_size_m
    final = np.argwhere(final_mask).astype(float) * voxel_size_m
    rotation = np.asarray(grain_to_lab, dtype=float)
    initial_crystal = initial @ rotation
    final_crystal = final @ rotation
    result: dict[str, float] = {}
    for axis, name in enumerate(("x", "y", "z")):
        negative = float(np.min(initial_crystal[:, axis]) - np.min(final_crystal[:, axis]))
        positive = float(np.max(final_crystal[:, axis]) - np.max(initial_crystal[:, axis]))
        result[f"growth_extent_{name}_negative_m"] = max(negative, 0.0)
        result[f"growth_extent_{name}_positive_m"] = max(positive, 0.0)
        result[f"growth_extent_{name}_total_m"] = max(negative, 0.0) + max(positive, 0.0)
    return result


def symmetric_chamfer_distance(
    first_points: NDArray[np.float64],
    second_points: NDArray[np.float64],
) -> float:
    """Mean symmetric nearest-neighbour distance for two non-empty point sets."""

    from scipy.spatial import cKDTree

    first = np.asarray(first_points, dtype=float)
    second = np.asarray(second_points, dtype=float)
    if first.ndim != 2 or second.ndim != 2 or first.shape[1] != second.shape[1]:
        raise ValueError("point sets must be two-dimensional with matching columns")
    if first.shape[0] == 0 or second.shape[0] == 0:
        raise ValueError("point sets must be non-empty")
    first_to_second = cKDTree(second).query(first, k=1)[0]
    second_to_first = cKDTree(first).query(second, k=1)[0]
    return 0.5 * (float(np.mean(first_to_second)) + float(np.mean(second_to_first)))


def closure_tolerance_kg(input_mass_kg: float, voxel_mass_kg: float) -> float:
    return 32.0 * math.ulp(max(abs(input_mass_kg), voxel_mass_kg))
