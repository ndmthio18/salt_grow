"""Mechanism-diagnosis metrics for Phase 2A.3b-2."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


def uniform_cube_facet_displacement_m(
    initial_side_m: float,
    added_volume_fraction: float,
) -> float:
    """Return the one-sided facet displacement for uniform cubic expansion."""

    if not math.isfinite(initial_side_m) or initial_side_m <= 0.0:
        raise ValueError("initial_side_m must be finite and positive")
    if not math.isfinite(added_volume_fraction) or added_volume_fraction < 0.0:
        raise ValueError("added_volume_fraction must be finite and nonnegative")
    return 0.5 * initial_side_m * (
        (1.0 + added_volume_fraction) ** (1.0 / 3.0) - 1.0
    )


def volume_fraction_for_facet_displacement(
    initial_side_m: float,
    displacement_m: float,
) -> float:
    """Invert :func:`uniform_cube_facet_displacement_m`."""

    if not math.isfinite(initial_side_m) or initial_side_m <= 0.0:
        raise ValueError("initial_side_m must be finite and positive")
    if not math.isfinite(displacement_m) or displacement_m < 0.0:
        raise ValueError("displacement_m must be finite and nonnegative")
    return (1.0 + 2.0 * displacement_m / initial_side_m) ** 3 - 1.0


def equivalent_displacement_m(
    event_count: float,
    voxel_size_m: float,
    physical_area_m2: float,
) -> float:
    """Convert fixed-volume voxel events to an area-averaged normal advance."""

    if not math.isfinite(event_count) or event_count < 0.0:
        raise ValueError("event_count must be finite and nonnegative")
    if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
        raise ValueError("voxel_size_m must be finite and positive")
    if not math.isfinite(physical_area_m2) or physical_area_m2 <= 0.0:
        raise ValueError("physical_area_m2 must be finite and positive")
    return event_count * voxel_size_m**3 / physical_area_m2


def group_allocation(
    candidate_counts: FloatArray,
    physical_areas_m2: FloatArray,
    mean_weights: FloatArray,
    *,
    allocation: str,
) -> FloatArray:
    """Return group probability for site- or physical-area-weighted allocation."""

    counts = np.asarray(candidate_counts, dtype=float)
    areas = np.asarray(physical_areas_m2, dtype=float)
    weights = np.asarray(mean_weights, dtype=float)
    if counts.ndim != 1 or areas.shape != counts.shape or weights.shape != counts.shape:
        raise ValueError("counts, areas, and weights must be matching vectors")
    if (
        np.any(~np.isfinite(counts))
        or np.any(~np.isfinite(areas))
        or np.any(~np.isfinite(weights))
        or np.any(counts <= 0.0)
        or np.any(areas <= 0.0)
        or np.any(weights <= 0.0)
    ):
        raise ValueError("allocation inputs must be finite and positive")
    if allocation == "site_weighted":
        attractiveness = counts * weights
    elif allocation == "area_weighted":
        attractiveness = areas * weights
    else:
        raise ValueError("allocation must be site_weighted or area_weighted")
    return attractiveness / np.sum(attractiveness)


def support_advances_m(
    initial_mask: NDArray[np.bool_],
    final_mask: NDArray[np.bool_],
    *,
    voxel_size_m: float,
    grain_to_lab: FloatArray,
    directions_crystal: FloatArray,
) -> FloatArray:
    """Return support-function advances along prescribed crystal-frame directions."""

    initial = np.asarray(initial_mask, dtype=bool)
    final = np.asarray(final_mask, dtype=bool)
    if initial.shape != final.shape or initial.ndim != 3:
        raise ValueError("initial and final masks must be matching 3D arrays")
    if not np.any(initial) or np.any(initial & ~final):
        raise ValueError("final mask must contain a non-empty initial mask")
    if not math.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
        raise ValueError("voxel_size_m must be finite and positive")
    rotation = np.asarray(grain_to_lab, dtype=float)
    if rotation.shape != (3, 3) or not np.all(np.isfinite(rotation)):
        raise ValueError("grain_to_lab must be a finite 3 by 3 matrix")
    directions = np.asarray(directions_crystal, dtype=float)
    if directions.ndim != 2 or directions.shape[1] != 3:
        raise ValueError("directions_crystal must have shape (n, 3)")
    norms = np.linalg.norm(directions, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= 0.0):
        raise ValueError("directions must be finite and nonzero")
    directions = directions / norms[:, None]
    initial_crystal = np.argwhere(initial).astype(float) * voxel_size_m @ rotation
    final_crystal = np.argwhere(final).astype(float) * voxel_size_m @ rotation
    initial_support = np.max(initial_crystal @ directions.T, axis=0)
    final_support = np.max(final_crystal @ directions.T, axis=0)
    advances = final_support - initial_support
    return np.maximum(advances, 0.0)
