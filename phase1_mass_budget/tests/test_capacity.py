import math

import pytest

from salt_phase1.capacity import diagnose_domain_capacity, required_domain_height_m


DX = 10.0e-6
GRID = (64, 64, 64)
DOMAIN_VOLUME = (64 * DX) ** 3


def test_capacity_overflow_is_explicit():
    result = diagnose_domain_capacity(
        grid_shape=GRID,
        voxel_size_m=DX,
        cumulative_precipitated_salt_volume_m3=1.2 * DOMAIN_VOLUME,
    )
    assert result.status == "capacity_exceeded"
    assert result.fits_permitted_capacity is False
    assert result.solid_volume_fraction == pytest.approx(1.2)
    assert result.overflow_solid_volume_m3 == pytest.approx(0.2 * DOMAIN_VOLUME)


@pytest.mark.parametrize(
    ("height_voxels", "expected_fraction", "expected_status"),
    [(32, 1.5, "capacity_exceeded"), (64, 0.75, "within_capacity"), (128, 0.375, "within_capacity")],
)
def test_multiple_domain_heights(height_voxels, expected_fraction, expected_status):
    inventory = 0.75 * DOMAIN_VOLUME
    result = diagnose_domain_capacity(
        grid_shape=(64, 64, height_voxels),
        voxel_size_m=DX,
        cumulative_precipitated_salt_volume_m3=inventory,
    )
    assert result.solid_volume_fraction == pytest.approx(expected_fraction)
    assert result.status == expected_status


@pytest.mark.parametrize(
    ("permitted_fraction", "fits"), [(0.3, False), (0.5, True), (0.7, True)]
)
def test_configurable_permitted_solid_fraction(permitted_fraction, fits):
    result = diagnose_domain_capacity(
        grid_shape=GRID,
        voxel_size_m=DX,
        cumulative_precipitated_salt_volume_m3=0.4 * DOMAIN_VOLUME,
        maximum_permitted_solid_volume_fraction=permitted_fraction,
    )
    assert result.fits_permitted_capacity is fits


def test_existing_seed_volume_is_included_in_capacity_and_reach_time():
    seed_voxels = 100
    rate = 2.0 * DX**3
    result = diagnose_domain_capacity(
        grid_shape=(10, 10, 10),
        voxel_size_m=DX,
        cumulative_precipitated_salt_volume_m3=200 * DX**3,
        initial_seed_voxel_count=seed_voxels,
        maximum_permitted_solid_volume_fraction=0.5,
        precipitation_volume_rate_m3_s=rate,
    )
    assert result.initial_seed_solid_volume_m3 == pytest.approx(seed_voxels * DX**3)
    assert result.total_solid_volume_m3 == pytest.approx(300 * DX**3)
    assert result.capacity_reach_time_s == pytest.approx(200.0)


def test_required_height_obeys_packing_fraction():
    inventory = 0.6 * DOMAIN_VOLUME
    height = required_domain_height_m(
        total_solid_volume_m3=inventory,
        base_grid_shape=(64, 64),
        voxel_size_m=DX,
        target_solid_volume_fraction=0.3,
    )
    assert height == pytest.approx(2.0 * 64 * DX)


def test_zero_rate_has_no_finite_future_capacity_time():
    result = diagnose_domain_capacity(
        grid_shape=GRID,
        voxel_size_m=DX,
        cumulative_precipitated_salt_volume_m3=0.0,
        precipitation_volume_rate_m3_s=0.0,
    )
    assert result.capacity_reach_time_s is None
    assert math.isfinite(result.domain_volume_m3)
