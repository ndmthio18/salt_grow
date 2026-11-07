"""Measurement and attribution diagnostics for Phase 2A.3b-4R.

The functions in this module are read-only with respect to the dynamic system.
They deliberately separate deposited fractional volume, binary geometry, and
surface-based front measurements.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

from .dynamic_interface_flux import DynamicInterfaceFluxSystem
from .surface_orientation import cubic_attachment_weight


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


def cell_centres(shape: tuple[int, int, int], voxel_size_m: float) -> FloatArray:
    """Return cell-centre coordinates with the cache's domain convention."""

    axes = [
        (-0.5 * count + 0.5 + np.arange(count, dtype=float)) * voxel_size_m
        for count in shape
    ]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)


def fixed_cylinder_mask(
    shape: tuple[int, int, int],
    voxel_size_m: float,
    normal: FloatArray,
    radius_m: float,
    axial_half_width_m: float,
) -> BoolArray:
    """Cell-centre approximation to a fixed physical cylindrical region."""

    points = cell_centres(shape, voxel_size_m)
    axial = points @ normal
    transverse = points - axial[..., None] * normal
    return (np.linalg.norm(transverse, axis=-1) <= radius_m) & (
        np.abs(axial) <= axial_half_width_m
    )


def weighted_quantile(values: FloatArray, weights: FloatArray, quantile: float) -> float:
    """Deterministic weighted quantile with midpoint cumulative positions."""

    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must lie in [0, 1]")
    if len(values) == 0 or len(values) != len(weights):
        raise ValueError("values and weights must be non-empty and aligned")
    order = np.argsort(values, kind="stable")
    ordered_values = values[order]
    ordered_weights = weights[order]
    total = float(np.sum(ordered_weights))
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("weights must have a finite positive sum")
    positions = (np.cumsum(ordered_weights) - 0.5 * ordered_weights) / total
    return float(np.interp(quantile, positions, ordered_values))


def front_height_statistics(
    system: DynamicInterfaceFluxSystem,
    name: str,
    radius_m: float,
) -> dict[str, float | int]:
    """Area-weighted reconstructed height distribution in a fixed window."""

    front = system.fronts[name]
    elements = system._full_surface(name)
    axial = elements.centroids_m @ front.reference_normal
    offsets = elements.centroids_m - front.centre_m
    transverse = offsets - (offsets @ front.reference_normal)[:, None] * front.reference_normal
    alignment = elements.normals @ front.reference_normal
    keep = (
        elements.normal_valid
        & (alignment >= front.minimum_alignment)
        & (np.linalg.norm(transverse, axis=1) <= radius_m)
    )
    if not np.any(keep):
        return {
            "face_count": 0,
            "area_m2": 0.0,
            "mean_m": float("nan"),
            "median_m": float("nan"),
            "std_m": float("nan"),
            "iqr_m": float("nan"),
            "p05_m": float("nan"),
            "p95_m": float("nan"),
        }
    heights = axial[keep]
    areas = elements.areas_m2[keep]
    total = float(np.sum(areas))
    mean = float(np.sum(areas * heights) / total)
    variance = float(np.sum(areas * (heights - mean) ** 2) / total)
    p05 = weighted_quantile(heights, areas, 0.05)
    q25 = weighted_quantile(heights, areas, 0.25)
    median = weighted_quantile(heights, areas, 0.50)
    q75 = weighted_quantile(heights, areas, 0.75)
    p95 = weighted_quantile(heights, areas, 0.95)
    return {
        "face_count": int(np.count_nonzero(keep)),
        "area_m2": total,
        "mean_m": mean,
        "median_m": median,
        "std_m": math.sqrt(max(0.0, variance)),
        "iqr_m": q75 - q25,
        "p05_m": p05,
        "p95_m": p95,
    }


def fixed_region_inventory(
    system: DynamicInterfaceFluxSystem,
    name: str,
    initial_mask: BoolArray,
    region_mask: BoolArray,
) -> dict[str, float | int]:
    """Realized and deposited inventories in one immutable observation region."""

    front = system.fronts[name]
    if initial_mask.shape != front.solid_mask.shape or region_mask.shape != front.solid_mask.shape:
        raise ValueError("measurement masks must match the front geometry")
    gained = int(np.count_nonzero(region_mask & front.solid_mask & ~initial_mask))
    lost = int(np.count_nonzero(region_mask & initial_mask & ~front.solid_mask))
    realized_volume = (gained - lost) * system.voxel_volume_m3
    deposited_volume = 0.0
    residual_volume = 0.0
    for owner, cumulative in system.cumulative_local_allocated_m3.items():
        if owner[0] == name and region_mask[owner[1]]:
            deposited_volume += float(cumulative)
    for owner, residual in system.local_reservoirs_m3.items():
        if owner[0] == name and region_mask[owner[1]]:
            residual_volume += float(residual)
    return {
        "gained_voxels": gained,
        "lost_voxels": lost,
        "realized_volume_m3": realized_volume,
        "deposited_volume_m3": deposited_volume,
        "residual_fractional_volume_m3": residual_volume,
        "binary_lag_volume_m3": deposited_volume - realized_volume,
    }


def interface_filter_audit(
    system: DynamicInterfaceFluxSystem,
    name: str,
    initial_face_keys: Iterable[tuple[int, int, int, int]],
) -> dict[str, float | int]:
    """Exclusive accounting of active and excluded reconstructed interface area."""

    front = system.fronts[name]
    elements = system._full_surface(name)
    cache = system._surface_cache[name]
    keys = sorted(cache.geometry)
    initial = set(initial_face_keys)
    offsets = elements.centroids_m - front.centre_m
    axial = offsets @ front.reference_normal
    transverse = offsets - axial[:, None] * front.reference_normal
    radius_ok = np.ones(len(keys), dtype=bool)
    if front.radius_m is not None:
        radius_ok = np.linalg.norm(transverse, axis=1) <= front.radius_m
    alignment = elements.normals @ front.reference_normal
    valid = elements.normal_valid
    active = valid & radius_ok & (alignment >= front.minimum_alignment)
    invalid = ~valid
    boundary = valid & ~radius_ok
    misaligned = valid & radius_ok & (alignment < front.minimum_alignment)
    created = np.asarray([key not in initial for key in keys], dtype=bool)
    digital_alignment = np.asarray([
        float(np.dot(cache.geometry[key][2], front.reference_normal)) for key in keys
    ])
    areas = elements.areas_m2

    def area(mask: BoolArray) -> float:
        return float(np.sum(areas[mask]))

    finite_alignment = alignment[valid]
    return {
        "total_face_count": len(keys),
        "total_corrected_area_m2": area(np.ones(len(keys), dtype=bool)),
        "active_area_m2": area(active),
        "invalid_reconstruction_area_m2": area(invalid),
        "physical_window_excluded_area_m2": area(boundary),
        "normal_misalignment_excluded_area_m2": area(misaligned),
        "growth_created_total_area_m2": area(created),
        "growth_created_active_area_m2": area(created & active),
        "digital_staircase_raw_area_m2": len(keys) * system.voxel_size_m**2,
        "digital_alignment_below_threshold_raw_area_m2": float(
            np.count_nonzero(digital_alignment < front.minimum_alignment)
            * system.voxel_size_m**2
        ),
        "mean_reconstructed_alignment": (
            float(np.mean(finite_alignment)) if len(finite_alignment) else float("nan")
        ),
        "p05_reconstructed_alignment": (
            float(np.percentile(finite_alignment, 5)) if len(finite_alignment) else float("nan")
        ),
        "exclusive_area_closure_m2": area(np.ones(len(keys), dtype=bool))
        - area(invalid) - area(boundary) - area(misaligned) - area(active),
    }


def controlled_variant_weights(
    system: DynamicInterfaceFluxSystem,
    names: Iterable[str],
    analytical_area_m2: float,
    alpha: float,
) -> list[dict[str, float | str]]:
    """A/B/C frozen-geometry weights for analytical and reconstructed controls."""

    rows: list[dict[str, float | str]] = []
    rotation = system.grain_to_lab
    for name in names:
        front = system.fronts[name]
        elements = system._active_surface(name)
        analytical_kappa = float(
            cubic_attachment_weight(front.reference_normal[None, :], rotation, alpha)[0]
        )
        weight_a = analytical_area_m2 * analytical_kappa
        if len(elements.areas_m2):
            reconstructed_kappa = cubic_attachment_weight(elements.normals, rotation, alpha)
            scale = analytical_area_m2 / float(np.sum(elements.areas_m2))
            weight_b = float(np.sum(scale * elements.areas_m2 * reconstructed_kappa))
            weight_c = float(np.sum(elements.areas_m2 * reconstructed_kappa))
            reconstructed_area = float(np.sum(elements.areas_m2))
        else:
            weight_b = 0.0
            weight_c = 0.0
            reconstructed_area = 0.0
        rows.append({
            "front": name,
            "alpha": float(alpha),
            "variant_a_weight_m2": weight_a,
            "variant_b_weight_m2": weight_b,
            "variant_c_weight_m2": weight_c,
            "reconstructed_active_area_m2": reconstructed_area,
        })
    totals = {
        key: math.fsum(float(row[key]) for row in rows)
        for key in ("variant_a_weight_m2", "variant_b_weight_m2", "variant_c_weight_m2")
    }
    for row in rows:
        for variant in "abc":
            key = f"variant_{variant}_weight_m2"
            row[f"variant_{variant}_share"] = float(row[key]) / totals[key] if totals[key] else float("nan")
    return rows


def reservoir_rows(
    system: DynamicInterfaceFluxSystem,
    names: Iterable[str],
) -> list[dict[str, float | int | str | bool]]:
    """Return every current or historical owner, including zero-filled active owners."""

    active = system.active_owners()
    owners = active | set(system.local_reservoirs_m3) | set(system.cumulative_local_allocated_m3)
    rows: list[dict[str, float | int | str | bool]] = []
    shape = next(iter(system.fronts.values())).solid_mask.shape
    origin = -0.5 * np.asarray(shape, dtype=float) * system.voxel_size_m
    allowed = set(names)
    for name, coordinate in sorted(owners):
        if name not in allowed:
            continue
        front = system.fronts[name]
        centre = origin + (np.asarray(coordinate, dtype=float) + 0.5) * system.voxel_size_m
        offset = centre - front.centre_m
        axial = float(np.dot(offset, front.reference_normal))
        transverse = float(np.linalg.norm(offset - axial * front.reference_normal))
        residual = float(system.local_reservoirs_m3.get((name, coordinate), 0.0))
        cumulative = float(system.cumulative_local_allocated_m3.get((name, coordinate), 0.0))
        rows.append({
            "front": name,
            "i": coordinate[0],
            "j": coordinate[1],
            "k": coordinate[2],
            "axial_position_m": axial,
            "transverse_radius_m": transverse,
            "active": (name, coordinate) in active,
            "binary_solid": bool(front.solid_mask[coordinate]),
            "residual_fractional_volume_m3": residual,
            "partial_fill_fraction": residual / system.voxel_volume_m3,
            "cumulative_deposited_volume_m3": cumulative,
            "cumulative_deposited_voxel_equivalents": cumulative / system.voxel_volume_m3,
        })
    return rows
