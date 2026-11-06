import copy

import numpy as np
import pytest
from scipy import ndimage

from salt_phase2 import (
    Material,
    OrientationCandidateWeightModel,
    PhysicalState,
    SaltGrowthCA,
    map_candidate_normals,
    proper_cubic_rotations,
    rotation_matrix_xyz_degrees,
)
from salt_phase2.orientation_verification import geometry_catalogue


DX = 10.0e-6
RHO = 2160.0
VOXEL_MASS = RHO * DX**3
ALPHA = 1.5 * np.log(2.0)


def _cube_state() -> PhysicalState:
    shape = (18, 18, 18)
    material = np.full(shape, Material.EMPTY, dtype=np.uint8)
    material[:, :, 0] = Material.SUBSTRATE
    seed = np.zeros(shape, dtype=bool)
    seed[6:12, 6:12, 5:11] = True
    material[seed] = Material.SALT
    liquid = seed.copy()
    state = PhysicalState(material, liquid, seed.copy(), 1)
    state.validate()
    return state


def test_rotated_plane_candidate_normal_is_reconstructed_not_lattice_direction():
    geometry = geometry_catalogue()["plane_33deg"]
    salt = geometry.voxelize(20.0e-6)
    candidates = ndimage.binary_dilation(salt) & ~salt
    coordinates = [tuple(int(v) for v in row) for row in np.argwhere(candidates)[::200]][:8]
    records = map_candidate_normals(
        salt,
        coordinates,
        voxel_size_m=20.0e-6,
        fitting_radius_m=120.0e-6,
        grain_to_lab=np.eye(3),
        alpha=ALPHA,
        edge_policy="smooth_average",
    )
    expected = np.array([np.cos(np.radians(33)), np.sin(np.radians(33)), 0.0])
    for record in records:
        normal = np.asarray(record.combined_normal)
        assert abs(float(np.dot(normal, expected))) > 0.98
        assert np.count_nonzero(np.abs(normal) > 0.1) >= 2


def test_edge_policies_are_explicit_and_preserve_alpha_zero_weights():
    salt = np.zeros((11, 11, 11), dtype=bool)
    salt[:5, :, :] = True
    salt[:, :5, :] = True
    coordinate = (5, 5, 5)
    smooth = map_candidate_normals(
        salt,
        [coordinate],
        voxel_size_m=DX,
        fitting_radius_m=35.0e-6,
        grain_to_lab=np.eye(3),
        alpha=0.0,
        edge_policy="smooth_average",
        quality_minimum=0.0,
    )[0]
    fallback = map_candidate_normals(
        salt,
        [coordinate],
        voxel_size_m=DX,
        fitting_radius_m=35.0e-6,
        grain_to_lab=np.eye(3),
        alpha=0.0,
        edge_policy="isotropic_ambiguous",
        quality_minimum=0.0,
    )[0]
    assert smooth.site_class == fallback.site_class == "edge_like"
    assert smooth.anisotropy_eligible
    assert not fallback.anisotropy_eligible
    assert smooth.anisotropy_weight == fallback.anisotropy_weight == 1.0


def test_cubic_equivalent_orientations_give_same_candidate_weights():
    state = _cube_state()
    salt = state.material == Material.SALT
    candidates = ndimage.binary_dilation(salt) & (state.material == Material.EMPTY)
    coordinates = [tuple(int(v) for v in row) for row in np.argwhere(candidates)]
    base = map_candidate_normals(
        salt,
        coordinates,
        voxel_size_m=DX,
        fitting_radius_m=40.0e-6,
        grain_to_lab=np.eye(3),
        alpha=ALPHA,
        edge_policy="smooth_average",
    )
    symmetry = proper_cubic_rotations()[7]
    equivalent = map_candidate_normals(
        salt,
        coordinates,
        voxel_size_m=DX,
        fitting_radius_m=40.0e-6,
        grain_to_lab=symmetry,
        alpha=ALPHA,
        edge_policy="smooth_average",
    )
    assert [record.anisotropy_weight for record in equivalent] == pytest.approx(
        [record.anisotropy_weight for record in base]
    )


def test_alpha_zero_orientation_mode_reproduces_isotropic_sequential_ca():
    baseline = SaltGrowthCA(
        _cube_state(),
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=712,
    )
    model = OrientationCandidateWeightModel(
        voxel_size_m=DX,
        fitting_radius_m=40.0e-6,
        grain_to_lab=rotation_matrix_xyz_degrees(17.0, 23.0, 11.0),
        alpha=0.0,
        edge_policy="smooth_average",
        update_strategy="local",
    )
    oriented = SaltGrowthCA(
        _cube_state(),
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=712,
        candidate_weight_model=model,
    )
    baseline_result = baseline.grow(6.5 * VOXEL_MASS)
    oriented_result = oriented.grow(6.5 * VOXEL_MASS)
    assert baseline_result.selected_coordinates == oriented_result.selected_coordinates
    assert np.array_equal(baseline.state.material, oriented.state.material)
    assert baseline.allocator.snapshot_state() == oriented.allocator.snapshot_state()


def test_orientation_cache_rolls_back_with_failed_transaction(monkeypatch):
    model = OrientationCandidateWeightModel(
        voxel_size_m=DX,
        fitting_radius_m=40.0e-6,
        grain_to_lab=np.eye(3),
        alpha=ALPHA,
        edge_policy="smooth_average",
        update_strategy="local",
    )
    ca = SaltGrowthCA(
        _cube_state(),
        solid_density_kg_m3=RHO,
        voxel_size_m=DX,
        growth_mode="stochastic",
        update_mode="sequential",
        rng_seed=9,
        candidate_weight_model=model,
    )
    material_before = ca.state.material.copy()
    model_before = copy.deepcopy(model.snapshot_state())

    def fail(*args, **kwargs):
        raise RuntimeError("injected orientation commit failure")

    monkeypatch.setattr(ca.allocator, "commit", fail)
    with pytest.raises(RuntimeError, match="orientation commit"):
        ca.grow(4.5 * VOXEL_MASS)
    assert np.array_equal(ca.state.material, material_before)
    assert model.snapshot_state() == model_before
    ca.assert_invariants()
