import numpy as np

from salt_phase2 import SaltGrowthCA, initialize_physical_state


RHO = 2160.0
DX = 10.0e-6


def _two_seed_ca(*, enforce_wetting=True):
    state = initialize_physical_state(
        shape=(7, 7, 6),
        substrate_layers=1,
        seed_coordinates=((1, 1, 1), (5, 5, 1)),
        wetted_seed_indices=(0,),
    )
    return SaltGrowthCA(
        state,
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        enforce_wetting=enforce_wetting,
    )


def test_only_component_of_wetted_seed_is_supply_connected():
    ca = _two_seed_ca()
    connected = ca.supply_connected_salt_mask()
    assert connected[1, 1, 1]
    assert not connected[5, 5, 1]
    assert np.count_nonzero(connected) == 1


def test_air_target_adjacent_to_wetted_salt_is_eligible():
    ca = _two_seed_ca()
    target = (1, 1, 2)
    assert not ca.state.liquid_accessible[target]
    assert ca.eligible_growth_mask()[target]


def test_target_near_unwetted_disconnected_seed_is_ineligible():
    ca = _two_seed_ca()
    assert not ca.eligible_growth_mask()[5, 5, 2]


def test_disabling_wetting_filter_is_an_unphysical_control():
    ca = _two_seed_ca(enforce_wetting=False)
    assert ca.eligible_growth_mask()[1, 1, 2]
    assert ca.eligible_growth_mask()[5, 5, 2]


def test_face_connectivity_is_not_periodic():
    state = initialize_physical_state(
        shape=(5, 5, 5),
        substrate_layers=1,
        seed_coordinates=((0, 2, 1),),
    )
    ca = SaltGrowthCA(state, solid_density_kg_m3=RHO, voxel_size_m=DX)
    assert ca.eligible_growth_mask()[1, 2, 1]
    assert not ca.eligible_growth_mask()[4, 2, 1]
