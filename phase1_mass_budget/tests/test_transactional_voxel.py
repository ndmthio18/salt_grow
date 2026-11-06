from dataclasses import replace
import math

import pytest

from salt_phase1.evaporation import EvaporationParameters, compute_mass_budget
from salt_phase1.voxel import TransactionalMassToVoxelAllocator


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def test_prepare_does_not_commit_before_geometry_confirmation():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    proposal = allocator.prepare(3.25 * VOXEL_MASS)
    assert proposal.requested_solid_voxels == 3
    assert allocator.cumulative_committed_solid_mass_kg == 0.0
    assert allocator.cumulative_realized_solid_voxels == 0


def test_mass_request_larger_than_spatial_capacity_is_retained():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    proposal = allocator.prepare(5.5 * VOXEL_MASS)
    result = allocator.commit(proposal, realized_solid_voxels=2)
    assert result.requested_solid_voxels == 5
    assert result.realized_solid_voxels == 2
    assert result.status == "partial_spatial_capacity"
    assert result.unrealized_requested_mass_kg == pytest.approx(3.0 * VOXEL_MASS)
    assert result.unallocated_mass_kg == pytest.approx(3.5 * VOXEL_MASS)
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=1.0e-26)


def test_zero_spatial_capacity_returns_explicit_status_and_no_fictitious_solid():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    proposal = allocator.prepare(4.0 * VOXEL_MASS)
    result = allocator.commit(proposal, realized_solid_voxels=0, unallocated_state="solution")
    assert result.status == "no_spatial_capacity"
    assert result.committed_solid_mass_kg == 0.0
    assert result.mass_remaining_in_solution_kg == pytest.approx(4.0 * VOXEL_MASS)


def test_zero_growth_below_one_voxel_is_explicit_and_carried_forward():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    first = allocator.prepare(0.4 * VOXEL_MASS)
    first_result = allocator.commit(first, realized_solid_voxels=0)
    assert first_result.status == "zero_growth_insufficient_mass"
    second = allocator.prepare(0.7 * VOXEL_MASS)
    assert second.requested_solid_voxels == 1
    second_result = allocator.commit(second, realized_solid_voxels=1)
    assert second_result.subvoxel_mass_kg == pytest.approx(0.1 * VOXEL_MASS)
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=2.0e-26)


def test_solution_mass_requires_explicit_release_before_reconsideration():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    proposal = allocator.prepare(2.0 * VOXEL_MASS)
    allocator.commit(proposal, realized_solid_voxels=0, unallocated_state="solution")
    zero_request = allocator.prepare(0.0)
    assert zero_request.requested_solid_voxels == 0
    allocator.commit(zero_request, realized_solid_voxels=0)
    released = allocator.prepare(0.0, released_solution_mass_kg=2.0 * VOXEL_MASS)
    result = allocator.commit(released, realized_solid_voxels=2)
    assert result.status == "committed"
    assert allocator.mass_remaining_in_solution_kg == pytest.approx(0.0)
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=1.0e-26)


def test_variable_rh_mass_budget_closes_cumulatively():
    base = EvaporationParameters(
        temperature_k=298.15,
        water_activity=0.753,
        relative_humidity=0.30,
        mass_transfer_coefficient_m_s=0.005,
        evaporation_area_m2=(64 * DX) ** 2,
        timestep_s=60.0,
        salt_mass_fraction=0.264,
        solid_density_kg_m3=RHO,
    )
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    masses = []
    for rh in (0.30, 0.753, 0.50, 0.80, 0.70, 0.30):
        mass = compute_mass_budget(replace(base, relative_humidity=rh)).precipitated_salt_mass_kg
        masses.append(mass)
        proposal = allocator.prepare(mass)
        allocator.commit(
            proposal, realized_solid_voxels=proposal.requested_solid_voxels
        )
    assert allocator.cumulative_input_mass_kg == pytest.approx(math.fsum(masses))
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=5.0e-25)


def test_stale_or_over_realized_transaction_is_rejected():
    allocator = TransactionalMassToVoxelAllocator(
        solid_density_kg_m3=RHO, voxel_size_m=DX
    )
    proposal = allocator.prepare(VOXEL_MASS)
    with pytest.raises(ValueError):
        allocator.commit(proposal, realized_solid_voxels=2)
