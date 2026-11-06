import numpy as np
import pytest

from salt_phase2 import AttachmentWeights, SaltGrowthCA, initialize_physical_state


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def _make_ca(
    *,
    rng_seed=123,
    weights=AttachmentWeights(),
    update_mode="synchronous",
):
    state = initialize_physical_state(
        shape=(14, 14, 12),
        substrate_layers=1,
        seed_coordinates=((4, 4, 1), (9, 9, 1)),
    )
    return SaltGrowthCA(
        state,
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        growth_mode="stochastic",
        update_mode=update_mode,
        attachment_weights=weights,
        rng_seed=rng_seed,
    )


def _run_partition(update_mode, increments, *, seed=123, weights=AttachmentWeights()):
    ca = _make_ca(rng_seed=seed, weights=weights, update_mode=update_mode)
    for increment in increments:
        ca.grow(increment)
    return ca


def test_stochastic_growth_is_reproducible_for_fixed_seed():
    first = _run_partition("synchronous", [30.25 * VOXEL_MASS] * 3)
    second = _run_partition("synchronous", [30.25 * VOXEL_MASS] * 3)
    assert np.array_equal(first.state.material, second.state.material)
    assert first.rng.bit_generator.state == second.rng.bit_generator.state


def test_weight_sets_change_geometry_without_changing_committed_mass():
    uniform = _run_partition(
        "sequential", [100.5 * VOXEL_MASS], weights=AttachmentWeights(1, 1, 1, 1)
    )
    preferential = _run_partition(
        "sequential", [100.5 * VOXEL_MASS], weights=AttachmentWeights(1, 4, 16, 64)
    )
    assert uniform.allocator.cumulative_realized_solid_voxels == 100
    assert preferential.allocator.cumulative_realized_solid_voxels == 100
    assert uniform.allocator.cumulative_committed_solid_mass_kg == pytest.approx(
        preferential.allocator.cumulative_committed_solid_mass_kg
    )
    assert not np.array_equal(uniform.state.material, preferential.state.material)


def test_sequential_result_depends_on_event_sequence_not_timestep_partition():
    total_mass = 50.5 * VOXEL_MASS
    partitions = {
        30: [total_mass / 4] * 4,
        60: [total_mass / 2] * 2,
        120: [total_mass],
    }
    results = {
        timestep: _run_partition("sequential", increments, seed=991)
        for timestep, increments in partitions.items()
    }
    assert all(
        ca.allocator.cumulative_realized_solid_voxels == 50
        for ca in results.values()
    )
    assert np.array_equal(results[30].state.material, results[60].state.material)
    assert np.array_equal(results[30].state.material, results[120].state.material)


def test_synchronous_batch_geometry_can_depend_on_timestep_partition():
    total_mass = 50.5 * VOXEL_MASS
    fine = _run_partition("synchronous", [total_mass / 4] * 4, seed=991)
    coarse = _run_partition("synchronous", [total_mass], seed=991)
    assert fine.allocator.cumulative_realized_solid_voxels == 50
    assert coarse.allocator.cumulative_realized_solid_voxels == 10
    assert coarse.allocator.deferred_precipitation_mass_kg > 0.0
    assert not np.array_equal(fine.state.material, coarse.state.material)


@pytest.mark.parametrize("weights", [(0, 1, 1, 1), (1, -1, 1, 1), (1, 1, np.nan, 1)])
def test_attachment_weights_must_be_positive_and_finite(weights):
    with pytest.raises(ValueError):
        AttachmentWeights(*weights)
