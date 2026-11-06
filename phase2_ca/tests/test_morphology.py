import numpy as np

from salt_phase2 import (
    Material,
    analyze_morphology,
    column_envelope,
    face_shell_envelope,
    initialize_physical_state,
)


def test_salt_cluster_count_and_sizes_distinguish_disconnected_seeds():
    state = initialize_physical_state(
        shape=(9, 9, 7),
        substrate_layers=1,
        seed_coordinates=((2, 4, 1), (6, 4, 1)),
    )
    disconnected = analyze_morphology(state)
    assert disconnected.salt_cluster_count == 2
    assert disconnected.salt_cluster_sizes == (1, 1)
    assert disconnected.seed_occupied_cluster_count == 2

    state.material[3:6, 4, 1] = Material.SALT
    connected = analyze_morphology(state)
    assert connected.salt_cluster_count == 1
    assert connected.salt_cluster_sizes == (5,)
    assert connected.seed_occupied_cluster_count == 1


def test_enclosed_void_is_not_confused_with_external_air():
    state = initialize_physical_state(
        shape=(7, 7, 7),
        substrate_layers=1,
        seed_coordinates=((3, 3, 1),),
    )
    for x in range(2, 5):
        for y in range(2, 5):
            for z in range(2, 5):
                if x in {2, 4} or y in {2, 4} or z in {2, 4}:
                    state.material[x, y, z] = Material.SALT
    analysis = analyze_morphology(state)
    assert analysis.bounding_box.enclosed_void_voxels == 1
    assert analysis.bounding_box.enclosed_void_component_sizes == (1,)
    assert analysis.bounding_box.external_air_connected_void_voxels > 0
    assert analysis.external_air_voxels > analysis.bounding_box.void_voxels


def test_envelope_definitions_report_sensitivity_not_one_intrinsic_porosity():
    state = initialize_physical_state(
        shape=(8, 8, 8),
        substrate_layers=1,
        seed_coordinates=((2, 2, 1),),
    )
    state.material[2, 2, 2:5] = Material.SALT
    state.material[5, 5, 1:3] = Material.SALT
    analysis = analyze_morphology(state)
    assert analysis.bounding_box.envelope_voxels > analysis.column_envelope.envelope_voxels
    assert (
        analysis.bounding_box.salt_volume_fraction
        != analysis.column_envelope.salt_volume_fraction
    )
    assert not analysis.intrinsic_porosity_identifiable
    assert analysis.tested_envelope_void_fraction_range > 0.0


def test_computational_percolation_is_separate_from_salt_envelope_pathway():
    state = initialize_physical_state(
        shape=(7, 7, 7),
        substrate_layers=1,
        seed_coordinates=((3, 3, 1),),
    )
    analysis = analyze_morphology(state, voxel_size_m=10.0e-6)
    assert analysis.computational_box.bottom_to_top_void_percolates
    assert analysis.computational_box.bottom_to_lateral_void_percolates
    assert not analysis.bounding_box.bottom_to_top_void_percolates
    assert analysis.bounding_box.accessible_pore_volume_m3 == 0.0


def test_accessible_and_enclosed_pore_volumes_use_explicit_voxel_scale():
    state = initialize_physical_state(
        shape=(7, 7, 7),
        substrate_layers=1,
        seed_coordinates=((3, 3, 1),),
    )
    for x in range(2, 5):
        for y in range(2, 5):
            for z in range(2, 5):
                if x in {2, 4} or y in {2, 4} or z in {2, 4}:
                    state.material[x, y, z] = Material.SALT
    voxel_size = 5.0e-6
    analysis = analyze_morphology(state, voxel_size_m=voxel_size)
    assert analysis.bounding_box.enclosed_void_voxels == 1
    assert analysis.bounding_box.enclosed_void_volume_m3 == voxel_size**3
    assert analysis.salt_surface_area_m2 is not None
    assert analysis.surface_to_volume_ratio_m_inv is not None
    assert analysis.characteristic_solid_size_m is not None


def test_column_envelope_implicit_global_top_has_no_empty_voxels_by_construction():
    state = initialize_physical_state(
        shape=(8, 8, 8),
        substrate_layers=1,
        seed_coordinates=((2, 2, 1), (5, 5, 1)),
    )
    state.material[2, 2, 2:6] = Material.SALT
    state.material[5, 5, 2:4] = Material.SALT
    envelope = column_envelope(state)
    global_top = int(np.max(np.argwhere(envelope)[:, 2]))
    empty = state.material == Material.EMPTY
    assert not np.any(envelope[:, :, global_top] & empty[:, :, global_top])


def test_physical_shell_thickness_is_grid_consistent():
    extents_um = []
    for voxel_size_um in (5, 10, 20):
        side = 20 // voxel_size_um
        shape = tuple(200 // voxel_size_um for _ in range(3))
        substrate_layers = 20 // voxel_size_um
        lower = 80 // voxel_size_um
        seed_coordinates = tuple(
            (x, y, z)
            for x in range(lower, lower + side)
            for y in range(lower, lower + side)
            for z in range(substrate_layers, substrate_layers + side)
        )
        state = initialize_physical_state(
            shape=shape,
            substrate_layers=substrate_layers,
            seed_coordinates=seed_coordinates,
        )
        shell = face_shell_envelope(
            state,
            voxel_size_m=voxel_size_um * 1.0e-6,
            shell_thickness_m=20.0e-6,
        )
        coordinates = np.argwhere(shell)
        extents_um.append(
            tuple(
                (int(np.ptp(coordinates[:, axis])) + 1) * voxel_size_um
                for axis in range(3)
            )
        )
    assert extents_um == [(60, 60, 40)] * 3


def test_intrinsic_porosity_requires_explicit_material_boundary():
    state = initialize_physical_state(
        shape=(7, 7, 7),
        substrate_layers=1,
        seed_coordinates=((3, 3, 1),),
    )
    unresolved = analyze_morphology(state)
    assert not unresolved.intrinsic_porosity_identifiable
    assert unresolved.intrinsic_porosity is None

    material_mask = np.zeros(state.shape, dtype=bool)
    material_mask[2:5, 2:5, 1:4] = True
    defined = analyze_morphology(
        state,
        intrinsic_material_mask=material_mask,
        intrinsic_boundary_source="synthetic reference cube",
    )
    assert defined.intrinsic_porosity_identifiable
    assert defined.intrinsic_porosity == 26 / 27
