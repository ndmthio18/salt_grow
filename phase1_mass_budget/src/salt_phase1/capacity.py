"""Domain-capacity diagnostics for the future spatial salt model.

The functions in this module compare a scalar solid inventory with a prescribed
voxel domain.  They do not inspect geometry, connectivity, or local site
eligibility and therefore cannot replace Phase 2 occupancy checks.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


def _positive_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")


def _nonnegative_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")


def _validate_grid_shape(grid_shape: Sequence[int]) -> tuple[int, int, int]:
    if len(grid_shape) != 3:
        raise ValueError("grid_shape must contain exactly three dimensions")
    shape: list[int] = []
    for dimension in grid_shape:
        if isinstance(dimension, bool) or int(dimension) != dimension or dimension <= 0:
            raise ValueError("grid_shape dimensions must be positive integers")
        shape.append(int(dimension))
    return shape[0], shape[1], shape[2]


@dataclass(frozen=True)
class CapacityDiagnostic:
    """Scalar solid-volume capacity result for one prescribed domain."""

    grid_shape: tuple[int, int, int]
    voxel_size_m: float
    domain_volume_m3: float
    initial_seed_voxel_count: int
    initial_seed_solid_volume_m3: float
    cumulative_precipitated_salt_volume_m3: float
    total_solid_volume_m3: float
    solid_volume_fraction: float
    maximum_permitted_solid_volume_fraction: float
    permitted_solid_volume_m3: float
    remaining_permitted_solid_volume_m3: float
    overflow_solid_volume_m3: float
    fits_permitted_capacity: bool
    capacity_reach_time_s: float | None
    status: str
    note: str


def diagnose_domain_capacity(
    *,
    grid_shape: Sequence[int],
    voxel_size_m: float,
    cumulative_precipitated_salt_volume_m3: float,
    initial_seed_voxel_count: int = 0,
    maximum_permitted_solid_volume_fraction: float = 1.0,
    precipitation_volume_rate_m3_s: float | None = None,
) -> CapacityDiagnostic:
    """Compare seed plus precipitated solid volume with a domain-wide limit.

    ``capacity_reach_time_s`` is measured from the start of precipitation and
    assumes a constant supplied solid-volume rate.  It is ``None`` when the
    rate is omitted or zero before capacity is reached.  This time is an
    inventory extrapolation, not a prediction that a feasible geometry exists.
    """

    shape = _validate_grid_shape(grid_shape)
    _positive_finite("voxel_size_m", voxel_size_m)
    _nonnegative_finite(
        "cumulative_precipitated_salt_volume_m3",
        cumulative_precipitated_salt_volume_m3,
    )
    if isinstance(initial_seed_voxel_count, bool):
        raise ValueError("initial_seed_voxel_count must be a non-negative integer")
    if int(initial_seed_voxel_count) != initial_seed_voxel_count or initial_seed_voxel_count < 0:
        raise ValueError("initial_seed_voxel_count must be a non-negative integer")
    if not math.isfinite(maximum_permitted_solid_volume_fraction) or not (
        0.0 < maximum_permitted_solid_volume_fraction <= 1.0
    ):
        raise ValueError(
            "maximum_permitted_solid_volume_fraction must be in (0, 1]"
        )
    if precipitation_volume_rate_m3_s is not None:
        _nonnegative_finite(
            "precipitation_volume_rate_m3_s", precipitation_volume_rate_m3_s
        )

    voxel_volume = voxel_size_m**3
    domain_volume = math.prod(shape) * voxel_volume
    seed_volume = int(initial_seed_voxel_count) * voxel_volume
    total_solid_volume = seed_volume + cumulative_precipitated_salt_volume_m3
    permitted_volume = maximum_permitted_solid_volume_fraction * domain_volume
    tolerance = 16.0 * math.ulp(max(domain_volume, permitted_volume, total_solid_volume))
    overflow = max(total_solid_volume - permitted_volume, 0.0)
    remaining = max(permitted_volume - total_solid_volume, 0.0)
    fits = total_solid_volume <= permitted_volume + tolerance

    if total_solid_volume > permitted_volume + tolerance:
        status = "capacity_exceeded"
        note = (
            "The scalar solid inventory exceeds the permitted domain capacity; "
            "no salt mass may be silently discarded or converted to fictitious voxels."
        )
    elif math.isclose(total_solid_volume, permitted_volume, rel_tol=0.0, abs_tol=tolerance):
        status = "at_capacity"
        note = "The scalar solid inventory has reached the permitted domain capacity."
    else:
        status = "within_capacity"
        note = (
            "The scalar inventory fits the permitted capacity, but local geometric "
            "occupancy and connectivity remain untested."
        )

    capacity_available_after_seed = max(permitted_volume - seed_volume, 0.0)
    if seed_volume >= permitted_volume - tolerance:
        capacity_reach_time: float | None = 0.0
    elif precipitation_volume_rate_m3_s is None or precipitation_volume_rate_m3_s == 0.0:
        capacity_reach_time = None
    else:
        capacity_reach_time = (
            capacity_available_after_seed / precipitation_volume_rate_m3_s
        )

    return CapacityDiagnostic(
        grid_shape=shape,
        voxel_size_m=voxel_size_m,
        domain_volume_m3=domain_volume,
        initial_seed_voxel_count=int(initial_seed_voxel_count),
        initial_seed_solid_volume_m3=seed_volume,
        cumulative_precipitated_salt_volume_m3=cumulative_precipitated_salt_volume_m3,
        total_solid_volume_m3=total_solid_volume,
        solid_volume_fraction=total_solid_volume / domain_volume,
        maximum_permitted_solid_volume_fraction=maximum_permitted_solid_volume_fraction,
        permitted_solid_volume_m3=permitted_volume,
        remaining_permitted_solid_volume_m3=remaining,
        overflow_solid_volume_m3=overflow,
        fits_permitted_capacity=fits,
        capacity_reach_time_s=capacity_reach_time,
        status=status,
        note=note,
    )


def required_domain_height_m(
    *,
    total_solid_volume_m3: float,
    base_grid_shape: Sequence[int],
    voxel_size_m: float,
    target_solid_volume_fraction: float,
) -> float:
    """Return the continuous height required for a target bulk solid fraction."""

    if len(base_grid_shape) != 2:
        raise ValueError("base_grid_shape must contain exactly two dimensions")
    nx, ny, _ = _validate_grid_shape((base_grid_shape[0], base_grid_shape[1], 1))
    _positive_finite("voxel_size_m", voxel_size_m)
    _nonnegative_finite("total_solid_volume_m3", total_solid_volume_m3)
    if not math.isfinite(target_solid_volume_fraction) or not (
        0.0 < target_solid_volume_fraction <= 1.0
    ):
        raise ValueError("target_solid_volume_fraction must be in (0, 1]")
    base_area = nx * ny * voxel_size_m**2
    return total_solid_volume_m3 / (target_solid_volume_fraction * base_area)
