import math

import pytest

from salt_phase1.voxel import MassToVoxelAllocator, allocate_salt_mass


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def test_voxel_mass_and_integer_allocation():
    result = allocate_salt_mass(2.5 * VOXEL_MASS, 0.0, RHO, DX)
    assert result.voxel_mass_kg == pytest.approx(VOXEL_MASS)
    assert result.new_solid_voxels == 2
    assert result.allocated_mass_kg == pytest.approx(2.0 * VOXEL_MASS)
    assert result.residual_mass_kg == pytest.approx(0.5 * VOXEL_MASS)


def test_subvoxel_mass_is_carried_forward_without_loss():
    allocator = MassToVoxelAllocator(solid_density_kg_m3=RHO, voxel_size_m=DX)
    first = allocator.allocate(0.4 * VOXEL_MASS)
    second = allocator.allocate(0.7 * VOXEL_MASS)
    assert first.new_solid_voxels == 0
    assert second.new_solid_voxels == 1
    assert allocator.residual_mass_kg == pytest.approx(0.1 * VOXEL_MASS)
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=1.0e-26)


def test_cumulative_mass_identity_over_many_steps():
    allocator = MassToVoxelAllocator(solid_density_kg_m3=RHO, voxel_size_m=DX)
    increments = [0.03, 1.2, 0.001, 4.7, 0.4, 3.0]
    for multiplier in increments:
        allocator.allocate(multiplier * VOXEL_MASS)
    assert allocator.cumulative_allocated_mass_kg + allocator.residual_mass_kg == pytest.approx(
        math.fsum(multiplier * VOXEL_MASS for multiplier in increments), rel=2.0e-15
    )
    assert allocator.closure_error_kg() == pytest.approx(0.0, abs=2.0e-26)


@pytest.mark.parametrize("voxel_size", [1.0e-6, 10.0e-6, 50.0e-6])
def test_multiple_voxel_sizes(voxel_size):
    voxel_mass = RHO * voxel_size**3
    result = allocate_salt_mass(7.25 * voxel_mass, 0.0, RHO, voxel_size)
    assert result.new_solid_voxels == 7
    assert result.residual_mass_kg == pytest.approx(0.25 * voxel_mass)


def test_initial_seeds_are_not_new_precipitated_voxels():
    allocator = MassToVoxelAllocator(
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        initial_seed_voxel_count=64,
    )
    allocator.allocate(2.2 * VOXEL_MASS)
    assert allocator.initial_seed_voxel_count == 64
    assert allocator.cumulative_new_voxels == 2
    assert allocator.total_morphology_voxel_count == 66
    assert allocator.cumulative_allocated_mass_kg == pytest.approx(2.0 * VOXEL_MASS)


def test_zero_mass_produces_no_negative_or_invalid_count():
    result = allocate_salt_mass(0.0, 0.0, RHO, DX)
    assert result.new_solid_voxels == 0
    assert result.allocated_mass_kg == 0.0
    assert result.residual_mass_kg == 0.0


@pytest.mark.parametrize(
    "args",
    [
        (-1.0, 0.0, RHO, DX),
        (0.0, -1.0, RHO, DX),
        (0.0, 0.0, 0.0, DX),
        (0.0, 0.0, RHO, 0.0),
        (math.nan, 0.0, RHO, DX),
    ],
)
def test_invalid_allocation_inputs_are_rejected(args):
    with pytest.raises(ValueError):
        allocate_salt_mass(*args)


def test_invalid_seed_count_is_rejected():
    with pytest.raises(ValueError):
        MassToVoxelAllocator(
            solid_density_kg_m3=RHO,
            voxel_size_m=DX,
            initial_seed_voxel_count=1.5,
        )
