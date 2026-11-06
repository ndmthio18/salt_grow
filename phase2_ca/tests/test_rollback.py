import numpy as np
import pytest
import copy

from salt_phase2 import SaltGrowthCA, initialize_physical_state


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def _ca(*, update_mode="synchronous"):
    state = initialize_physical_state(
        shape=(8, 8, 8),
        substrate_layers=1,
        seed_coordinates=((4, 4, 1),),
    )
    return SaltGrowthCA(
        state,
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        growth_mode="stochastic",
        update_mode=update_mode,
        rng_seed=5,
    )


def test_candidate_selection_failure_aborts_pending_transaction(monkeypatch):
    ca = _ca()
    before = ca.state.material.copy()
    rng_before = copy.deepcopy(ca.rng.bit_generator.state)

    def fail(*args, **kwargs):
        raise RuntimeError("injected selection failure")

    monkeypatch.setattr(ca, "_synchronous_selection", fail)
    with pytest.raises(RuntimeError, match="selection"):
        ca.grow(3.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, before)
    assert not ca.allocator.has_pending_transaction
    assert ca.allocator.cumulative_input_mass_kg == 0.0
    assert ca.allocator.closure_error_kg() == 0.0
    assert ca.rng.bit_generator.state == rng_before


@pytest.mark.parametrize("update_mode", ["synchronous", "sequential"])
def test_commit_failure_rolls_back_geometry_and_aborts(monkeypatch, update_mode):
    ca = _ca(update_mode=update_mode)
    before = ca.state.material.copy()
    rng_before = copy.deepcopy(ca.rng.bit_generator.state)

    def fail(*args, **kwargs):
        raise RuntimeError("injected commit failure")

    monkeypatch.setattr(ca.allocator, "commit", fail)
    with pytest.raises(RuntimeError, match="commit"):
        ca.grow(3.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, before)
    assert not ca.allocator.has_pending_transaction
    assert ca.allocator.cumulative_input_mass_kg == 0.0
    assert ca.allocator.closure_error_kg() == 0.0
    assert ca.rng.bit_generator.state == rng_before


def test_sequential_stochastic_selection_failure_after_insertions_rolls_back(
    monkeypatch,
):
    ca = _ca(update_mode="sequential")
    before = ca.state.material.copy()
    rng_before = copy.deepcopy(ca.rng.bit_generator.state)
    original = ca._stochastic_candidate
    selection_count = 0

    def fail_on_fourth_selection(coordinates, coordination):
        nonlocal selection_count
        selection_count += 1
        if selection_count == 4:
            raise RuntimeError("injected stochastic selection failure")
        return original(coordinates, coordination)

    monkeypatch.setattr(ca, "_stochastic_candidate", fail_on_fourth_selection)
    with pytest.raises(RuntimeError, match="stochastic selection"):
        ca.grow(8.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, before)
    assert ca.allocator.cumulative_input_mass_kg == 0.0
    assert ca.allocator.cumulative_realized_solid_voxels == 0
    assert not ca.allocator.has_pending_transaction
    assert ca.rng.bit_generator.state == rng_before


@pytest.mark.parametrize("failure_insertion", [1, 4])
def test_sequential_insertion_failure_restores_complete_step(
    monkeypatch, failure_insertion
):
    ca = _ca(update_mode="sequential")
    before = ca.state.material.copy()
    rng_before = copy.deepcopy(ca.rng.bit_generator.state)
    allocator_before = ca.allocator.snapshot_state()
    original = ca._insert_tentative_voxel
    insertion_count = 0

    def insert_then_fail(coordinate):
        nonlocal insertion_count
        insertion_count += 1
        original(coordinate)
        if insertion_count == failure_insertion:
            raise RuntimeError(f"injected insertion {failure_insertion} failure")

    monkeypatch.setattr(ca, "_insert_tentative_voxel", insert_then_fail)
    with pytest.raises(RuntimeError, match="injected insertion"):
        ca.grow(8.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, before)
    assert ca.allocator.snapshot_state() == allocator_before
    assert not ca.allocator.has_pending_transaction
    assert ca.rng.bit_generator.state == rng_before
    ca.assert_invariants()


def test_sequential_candidate_generation_failure_after_insertions_rolls_back(
    monkeypatch,
):
    ca = _ca(update_mode="sequential")
    before = ca.state.material.copy()
    rng_before = copy.deepcopy(ca.rng.bit_generator.state)
    original = ca.eligible_growth_mask
    calls = 0

    def fail_on_fourth_generation():
        nonlocal calls
        calls += 1
        if calls == 4:
            raise RuntimeError("injected candidate generation failure")
        return original()

    monkeypatch.setattr(ca, "eligible_growth_mask", fail_on_fourth_generation)
    with pytest.raises(RuntimeError, match="candidate generation"):
        ca.grow(8.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, before)
    assert ca.allocator.cumulative_input_mass_kg == 0.0
    assert ca.allocator.cumulative_realized_solid_voxels == 0
    assert not ca.allocator.has_pending_transaction
    assert ca.rng.bit_generator.state == rng_before
