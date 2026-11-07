import math

import numpy as np
import pytest

from salt_phase2.dynamic_interface_flux import DynamicInterfaceFluxSystem, FrontRegion
from salt_phase2.front_error_attribution import (
    controlled_variant_weights,
    fixed_cylinder_mask,
    fixed_region_inventory,
    front_height_statistics,
    interface_filter_audit,
    weighted_quantile,
)


def make_system(cells=16, dx=5e-6):
    mask = np.zeros((cells, cells, cells), dtype=bool)
    mask[: cells // 2] = True
    front = FrontRegion(
        name="axis",
        solid_mask=mask,
        reference_normal=np.asarray([1.0, 0.0, 0.0]),
        centre_m=np.zeros(3),
        radius_m=25e-6,
        minimum_alignment=0.65,
    )
    return DynamicInterfaceFluxSystem(
        [front],
        voxel_size_m=dx,
        fitting_radius_m=15e-6,
        solid_density_kg_m3=2160.0,
        grain_to_lab=np.eye(3),
        alpha=0.0,
        rng_seed=13,
    )


def test_weighted_quantile_and_height_statistics_are_well_defined():
    assert weighted_quantile(
        np.asarray([0.0, 1.0, 2.0]), np.ones(3), 0.5
    ) == pytest.approx(1.0)
    stats = front_height_statistics(make_system(), "axis", 20e-6)
    assert stats["face_count"] > 0
    assert stats["mean_m"] == pytest.approx(0.0, abs=1e-15)
    assert stats["median_m"] == pytest.approx(0.0, abs=1e-15)
    assert stats["std_m"] == pytest.approx(0.0, abs=1e-15)


def test_fixed_region_separates_fractional_deposition_from_binary_geometry():
    system = make_system()
    initial = system.fronts["axis"].solid_mask.copy()
    region = fixed_cylinder_mask(
        initial.shape, system.voxel_size_m, np.asarray([1.0, 0.0, 0.0]), 20e-6, 25e-6
    )
    owner = next(iter(system.active_owners()))
    assert region[owner[1]]
    system.local_reservoirs_m3[owner] = 0.4 * system.voxel_volume_m3
    system.cumulative_local_allocated_m3[owner] = 0.4 * system.voxel_volume_m3
    record = fixed_region_inventory(system, "axis", initial, region)
    assert record["realized_volume_m3"] == 0.0
    assert record["deposited_volume_m3"] == pytest.approx(0.4 * system.voxel_volume_m3)
    assert record["binary_lag_volume_m3"] == pytest.approx(0.4 * system.voxel_volume_m3)


def test_filter_area_partition_closes_and_alpha_zero_controls_match():
    system = make_system()
    system.front_snapshot("axis")
    initial_keys = set(system._surface_cache["axis"].geometry)
    record = interface_filter_audit(system, "axis", initial_keys)
    assert record["exclusive_area_closure_m2"] == pytest.approx(0.0, abs=1e-24)
    assert record["growth_created_total_area_m2"] == 0.0
    area = math.pi * (25e-6) ** 2
    controls = controlled_variant_weights(system, ["axis"], area, alpha=0.0)
    assert controls[0]["variant_a_share"] == pytest.approx(1.0)
    assert controls[0]["variant_b_share"] == pytest.approx(1.0)
    assert controls[0]["variant_c_share"] == pytest.approx(1.0)
