import numpy as np
import pytest
import copy

from salt_phase2 import SaltGrowthCA, initialize_physical_state


RHO = 2160.0
DX = 10.0e-6
VOXEL_MASS = RHO * DX**3


def _ca():
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


def test_commit_failure_rolls_back_geometry_and_aborts(monkeypatch):
    ca = _ca()
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
