import copy
import math

import numpy as np
import pytest

from salt_phase2.interface_flux import (
    ConservativePatchLedger,
    allocate_interface_volume,
    bundle_elements_by_candidate,
    corrected_digital_interface,
    frozen_front_realization,
    triangulated_interface,
)


def central_plane(shape, normal):
    axes = [np.arange(count, dtype=float) - 0.5 * (count - 1) for count in shape]
    x, y, z = np.meshgrid(*axes, indexing="ij", sparse=True)
    normal = np.asarray(normal, dtype=float)
    normal /= np.linalg.norm(normal)
    return np.asarray(normal[0] * x + normal[1] * y + normal[2] * z <= 0.0)


def test_continuum_allocation_is_exact_and_candidate_count_independent():
    areas = np.asarray([1.0, math.sqrt(2.0), 3.0 * math.sqrt(3.0) / 4.0]) * 1.0e-8
    kappa = np.asarray([2.0, math.sqrt(2.0), 1.0])
    result = allocate_interface_volume(7.3e-15, areas, kappa)
    assert math.fsum(result.patch_volumes_m3) == 7.3e-15
    assert result.normal_advances_m / result.normal_advances_m[-1] == pytest.approx(
        kappa / kappa[-1], rel=2.0e-15
    )
    repeated = allocate_interface_volume(7.3e-15, areas, kappa)
    assert repeated.patch_volumes_m3 == pytest.approx(result.patch_volumes_m3, rel=0, abs=0)


def test_alpha_zero_equal_advance_for_unequal_areas():
    areas = np.asarray([1.0, 2.0, 5.0])
    result = allocate_interface_volume(8.0, areas, np.ones(3))
    assert result.normal_advances_m == pytest.approx(np.ones(3))


def test_surface_elements_are_nonoverlapping_and_candidate_bundling_closes():
    dx = 10.0e-6
    mask = central_plane((18, 18, 18), (1.0, 1.0, 0.0))
    digital = corrected_digital_interface(mask, dx, 35.0e-6)
    bundles = bundle_elements_by_candidate(digital)
    assert all(not mask[tuple(coordinate)] for coordinate in digital.candidate_coordinates)
    assert math.fsum(p.area_m2 for p in bundles) == pytest.approx(
        digital.total_area_m2, rel=2.0e-15
    )
    assert sum(p.element_count for p in bundles) == len(digital.areas_m2)
    assert any(p.element_count > 1 for p in bundles)


def test_triangulated_axis_plane_area_is_exact_on_cell_centre_support():
    dx = 10.0e-6
    mask = central_plane((18, 18, 18), (1.0, 0.0, 0.0))
    surface = triangulated_interface(mask, dx)
    expected = ((mask.shape[1] - 1) * dx) * ((mask.shape[2] - 1) * dx)
    assert surface.total_area_m2 == pytest.approx(expected, rel=2.0e-14)
    assert np.all(surface.normals[:, 0] > 0.99)


def test_signed_distance_triangulation_tracks_rotated_plane_without_grid_normals():
    dx = 10.0e-6
    truth = np.asarray([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    mask = central_plane((24, 24, 24), truth)
    surface = triangulated_interface(mask, dx)
    angles = np.degrees(
        np.arccos(np.clip(surface.normals @ truth, -1.0, 1.0))
    )
    assert np.percentile(angles, 90) < 35.0
    assert np.count_nonzero(np.abs(surface.normals) > 0.1, axis=1).mean() > 2.0


def test_patch_ledger_closes_prevents_double_consumption_and_rolls_back():
    ledger = ConservativePatchLedger(2160.0, 10.0e-6)
    voxel_mass = ledger.voxel_mass_kg
    update = ledger.update(
        2.5 * voxel_mass,
        ("a", "b"),
        (1.0, 1.0),
        (1.0, 1.0),
        (("shared", "a2"), ("shared", "b2")),
    )
    assert update.realized_voxels == 2
    assert len(ledger.consumed_candidates) == 2
    assert update.closure_error_kg == pytest.approx(0.0, abs=64 * math.ulp(voxel_mass))
    snapshot = copy.deepcopy(ledger.snapshot_state())
    with pytest.raises(RuntimeError, match="injected"):
        ledger.update(
            2.0 * voxel_mass,
            ("a",),
            (1.0,),
            (1.0,),
            (("a3", "a4"),),
            inject_failure_after_voxels=1,
        )
    assert ledger.snapshot_state() == snapshot


def test_frozen_front_integer_error_refines_and_mass_closes():
    areas = np.asarray([1.0, math.sqrt(2.0), 3.0 * math.sqrt(3.0) / 4.0]) * (320e-6) ** 2
    total = math.fsum(areas) * 1.0e-6
    errors = []
    for dx in (20e-6, 10e-6, 5e-6):
        capacities = np.ceil(areas / dx**2 * 2.0).astype(int)
        result = frozen_front_realization(total, areas, np.ones(3), dx, capacities)
        errors.append(float(np.max(np.abs(result["voxel_only_advances_m"] / 1.0e-6 - 1.0))))
        assert result["closure_error_m3"] == pytest.approx(0.0, abs=128 * math.ulp(total))
        assert result["capacity_limited_patches"] == 0
    assert errors[2] < errors[1] < errors[0]
