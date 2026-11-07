import copy
import math

import numpy as np
import pytest

from salt_phase2.dynamic_interface_flux import (
    DynamicInterfaceFluxSystem,
    FrontRegion,
    candidate_face_flux_audit,
    perpendicular_face_counterexample,
)
from salt_phase2.interface_flux import corrected_digital_interface


def axis_plane(cells):
    mask = np.zeros((cells, cells, cells), dtype=bool)
    mask[: cells // 2] = True
    return mask


def make_system(cells=14, dx=10e-6, radius=35e-6):
    front = FrontRegion(
        name="axis",
        solid_mask=axis_plane(cells),
        reference_normal=np.asarray([1.0, 0.0, 0.0]),
        centre_m=np.zeros(3),
        radius_m=radius,
        minimum_alignment=0.8,
    )
    return DynamicInterfaceFluxSystem(
        [front], voxel_size_m=dx, fitting_radius_m=25e-6,
        solid_density_kg_m3=2160.0, grain_to_lab=np.eye(3), alpha=0.0,
        rng_seed=17,
    )


def test_perpendicular_face_counterexample_is_exact():
    for alpha in (-1.5 * math.log(2.0), 0.0, 1.5 * math.log(2.0)):
        row = perpendicular_face_counterexample(alpha)
        assert row["adapter_to_face_ratio"] == pytest.approx(
            math.exp(-0.5 * alpha), rel=2e-15
        )


def test_face_flux_audit_detects_multiface_bundle_difference():
    mask = np.ones((12, 12, 12), dtype=bool)
    mask[6:, 6:, :] = False
    elements = corrected_digital_interface(mask, 10e-6, 25e-6)
    rows = candidate_face_flux_audit(elements, np.eye(3), 1.5 * math.log(2.0))
    multi = [row for row in rows if row["element_count"] > 1]
    assert multi
    assert max(float(row["absolute_relative_discrepancy"]) for row in multi) > 0.1


def test_fractional_inventory_is_timestep_partition_invariant():
    cumulative = make_system()
    split = make_system()
    initial_owners = cumulative.active_owners()
    initial_faces = set(cumulative._surface_cache["axis"].geometry)
    count = cumulative.front_snapshot("axis").active_candidate_count
    total_volume = count * 1.25 * cumulative.voxel_volume_m3
    total_mass = total_volume * cumulative.solid_density_kg_m3
    cumulative.advance(total_mass)
    for _ in range(5):
        split.advance(total_mass / 5.0)
    assert np.array_equal(
        cumulative.fronts["axis"].solid_mask,
        split.fronts["axis"].solid_mask,
    )
    assert cumulative.local_reservoirs_m3.keys() == split.local_reservoirs_m3.keys()
    assert np.asarray(list(cumulative.local_reservoirs_m3.values())) == pytest.approx(
        np.asarray(list(split.local_reservoirs_m3.values())), abs=1e-28
    )
    assert cumulative.closure_error_m3() == pytest.approx(0.0, abs=1e-28)
    assert split.closure_error_m3() == pytest.approx(0.0, abs=1e-28)
    final_owners = cumulative.active_owners()
    final_faces = set(cumulative._surface_cache["axis"].geometry)
    assert final_owners - initial_owners
    assert initial_faces - final_faces
    assert final_faces - initial_faces
    assert set(cumulative.cumulative_local_allocated_m3).issuperset(initial_owners)
    consistency = cumulative.local_cache_consistency("axis")
    assert consistency["keys_exact"]
    assert consistency["valid_flags_exact"]
    assert consistency["maximum_area_abs_error_m2"] == pytest.approx(0.0, abs=1e-28)
    assert consistency["maximum_normal_vector_error"] == pytest.approx(0.0, abs=1e-12)


def test_dynamic_failure_restores_geometry_inventory_rng_and_cache():
    system = make_system()
    count = system.front_snapshot("axis").active_candidate_count
    total_mass = count * 1.1 * system.voxel_volume_m3 * system.solid_density_kg_m3
    before = copy.deepcopy(system.snapshot_state())
    with pytest.raises(RuntimeError, match="injected"):
        system.advance(total_mass, inject_failure_after_event_batches=1)
    after = system.snapshot_state()
    assert np.array_equal(
        before["fronts"]["axis"].solid_mask,
        after["fronts"]["axis"].solid_mask,
    )
    assert before["local_reservoirs_m3"] == after["local_reservoirs_m3"]
    assert before["cumulative_input_volume_m3"] == after["cumulative_input_volume_m3"]
    assert before["rng_state"] == after["rng_state"]
    assert before["surface_cache"].keys() == after["surface_cache"].keys()


def test_domain_boundary_is_not_an_external_growth_candidate():
    system = make_system(cells=12, radius=1e-3)
    active = system.active_owners()
    shape = system.fronts["axis"].solid_mask.shape
    assert active
    assert all(0 <= coordinate[i] < shape[i] for _, coordinate in active for i in range(3))
    assert all(coordinate[0] != 0 for _, coordinate in active)
