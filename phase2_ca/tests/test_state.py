import numpy as np
import pytest

from salt_phase2 import Material, initialize_physical_state


def test_material_liquid_and_seed_state_are_separate_arrays():
    state = initialize_physical_state(
        shape=(8, 7, 6),
        substrate_layers=2,
        seed_coordinates=((1, 2, 2), (6, 4, 2)),
        wetted_seed_indices=(0,),
    )
    assert state.material.shape == (8, 7, 6)
    assert state.material.dtype == np.uint8
    assert state.liquid_accessible.dtype == np.bool_
    assert state.seed_mask.dtype == np.bool_
    assert not np.shares_memory(state.material, state.liquid_accessible)
    assert np.all(state.material[:, :, :2] == Material.SUBSTRATE)
    assert state.initial_seed_voxel_count == 2
    assert state.liquid_accessible[1, 2, 2]
    assert not state.liquid_accessible[6, 4, 2]


@pytest.mark.parametrize(
    "coordinates",
    [(), ((1, 1, 2), (1, 1, 2)), ((1, 1, 3),), ((9, 1, 2),)],
)
def test_invalid_seed_definitions_are_rejected(coordinates):
    with pytest.raises(ValueError):
        initialize_physical_state(
            shape=(8, 8, 6),
            substrate_layers=2,
            seed_coordinates=coordinates,
        )


def test_seed_must_remain_salt():
    state = initialize_physical_state(
        shape=(5, 5, 5),
        substrate_layers=1,
        seed_coordinates=((2, 2, 1),),
    )
    state.material[2, 2, 1] = Material.EMPTY
    with pytest.raises(ValueError, match="seed"):
        state.validate()
