import numpy as np
import pytest

from salt_phase2 import Material, SaltGrowthCA, initialize_physical_state


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def _ca(*, wet=True, enforce_wetting=True):
    state = initialize_physical_state(
        shape=(9, 9, 8),
        substrate_layers=1,
        seed_coordinates=((2, 2, 1), (6, 6, 1)),
        wetted_seed_indices=(0, 1) if wet else (),
    )
    return SaltGrowthCA(
        state,
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        enforce_wetting=enforce_wetting,
    )


def test_growth_commits_only_realized_voxels_and_separates_seed_mass():
    ca = _ca()
    result = ca.grow(2.5 * VOXEL_MASS)
    assert result.transaction.requested_solid_voxels == 2
    assert result.transaction.realized_solid_voxels == 2
    assert result.morphology.seed_salt_voxels == 2
    assert result.morphology.new_salt_voxels == 2
    assert result.morphology.total_salt_voxels == 4
    assert result.morphology.seed_salt_mass_kg == pytest.approx(2 * VOXEL_MASS)
    assert result.morphology.new_salt_mass_kg == pytest.approx(2 * VOXEL_MASS)
    assert result.morphology.mass_closure_error_kg == pytest.approx(0.0, abs=1.0e-26)


def test_request_above_geometric_capacity_retains_all_unrealized_mass():
    ca = _ca()
    candidate_count = int(np.count_nonzero(ca.eligible_growth_mask()))
    result = ca.grow(100.5 * VOXEL_MASS)
    assert result.candidate_voxel_count == candidate_count
    assert result.transaction.realized_solid_voxels == candidate_count
    assert result.transaction.requested_solid_voxels > candidate_count
    assert result.transaction.status == "partial_spatial_capacity"
    assert result.morphology.new_salt_voxels == candidate_count
    assert result.morphology.retained_unallocated_mass_kg > 0.0
    assert result.morphology.mass_closure_error_kg == pytest.approx(0.0, abs=1.0e-24)
    assert np.all(ca.state.material[:, :, 0] == Material.SUBSTRATE)


def test_no_wetted_seed_gives_no_growth_and_explicit_retained_mass():
    ca = _ca(wet=False)
    before = ca.state.material.copy()
    result = ca.grow(3.5 * VOXEL_MASS)
    assert result.transaction.status == "no_spatial_capacity"
    assert result.transaction.realized_solid_voxels == 0
    assert np.array_equal(ca.state.material, before)
    assert result.morphology.retained_unallocated_mass_kg == pytest.approx(
        3.5 * VOXEL_MASS
    )


def test_supply_connectivity_is_recomputed_after_growth():
    ca = _ca()
    first = ca.grow(4.5 * VOXEL_MASS)
    connected = ca.supply_connected_salt_mask()
    assert np.count_nonzero(connected) == first.morphology.total_salt_voxels
    second_candidates = ca.eligible_growth_mask()
    assert np.count_nonzero(second_candidates) > 0


def test_candidate_selection_is_deterministic():
    first = _ca()
    second = _ca()
    first_result = first.grow(5.5 * VOXEL_MASS)
    second_result = second.grow(5.5 * VOXEL_MASS)
    assert first_result.selected_coordinates == second_result.selected_coordinates
    assert np.array_equal(first.state.material, second.state.material)


def test_external_unaccounted_salt_is_detected():
    ca = _ca()
    ca.state.material[4, 4, 2] = Material.SALT
    with pytest.raises(RuntimeError, match="ledger"):
        ca.assert_invariants()
