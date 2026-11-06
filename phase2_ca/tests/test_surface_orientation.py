import numpy as np
import pytest

from salt_phase2 import (
    angular_errors_deg,
    crystal_frame_normals,
    cubic_attachment_weight,
    cubic_q4,
    cube_reference_normals,
    estimate_plane_fit_normals,
    estimate_sdf_normals,
    interface_face_centres,
    proper_cubic_rotations,
    rotation_matrix_xyz_degrees,
    signed_distance_from_mask,
)


def _axis_plane(n=32, dx=10.0e-6):
    centres = (-0.5 * n + 0.5 + np.arange(n)) * dx
    solid = np.broadcast_to(centres[:, None, None] < 0.0, (n, n, n)).copy()
    return solid, dx


def test_signed_distance_sign_and_axis_plane_outward_normal():
    solid, dx = _axis_plane()
    phi = signed_distance_from_mask(solid, dx)
    assert np.all(phi[solid] < 0.0)
    assert np.all(phi[~solid] > 0.0)
    estimate = estimate_sdf_normals(solid, dx, 20.0e-6, max_points=300)
    truth = np.broadcast_to([1.0, 0.0, 0.0], estimate.normals.shape)
    errors = angular_errors_deg(estimate.normals, truth)
    assert np.nanmax(errors) < 1.0e-5


def test_weighted_plane_fit_uses_physical_radius_and_recovers_plane():
    solid, dx = _axis_plane()
    estimate = estimate_plane_fit_normals(solid, dx, 35.0e-6, max_points=200)
    truth = np.broadcast_to([1.0, 0.0, 0.0], estimate.normals.shape)
    errors = angular_errors_deg(estimate.normals, truth)
    assert np.count_nonzero(estimate.valid) == estimate.points_m.shape[0]
    assert np.nanmax(errors) < 1.0e-5
    with pytest.raises(ValueError, match="exceed one voxel"):
        estimate_plane_fit_normals(solid, dx, dx)


def test_cubic_q4_and_predeclared_alpha_invariants():
    directions = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )
    q4 = cubic_q4(directions)
    assert q4 == pytest.approx([1.0, 0.5, 1.0 / 3.0])
    identity = np.eye(3)
    assert cubic_attachment_weight(directions, identity, 0.0) == pytest.approx(1.0)
    alpha = 1.5 * np.log(2.0)
    positive = cubic_attachment_weight(directions, identity, alpha)
    negative = cubic_attachment_weight(directions, identity, -alpha)
    assert positive[0] > positive[1] > positive[2]
    assert negative[0] < negative[1] < negative[2]


def test_all_24_proper_cubic_rotations_preserve_weight():
    rotations = proper_cubic_rotations()
    assert len(rotations) == 24
    normals = np.array([[0.2, -0.4, 0.7], [1.0, 2.0, -3.0]])
    alpha = 1.5 * np.log(2.0)
    baseline = cubic_attachment_weight(normals, np.eye(3), alpha)
    for symmetry in rotations:
        equivalent_grain_to_lab = symmetry
        assert cubic_attachment_weight(
            normals, equivalent_grain_to_lab, alpha
        ) == pytest.approx(baseline)


def test_joint_rotation_preserves_crystal_normal_and_weight():
    grain = rotation_matrix_xyz_degrees(13.0, -21.0, 7.0)
    joint = rotation_matrix_xyz_degrees(-17.0, 11.0, 29.0)
    normals = np.array([[1.0, 2.0, 3.0], [-2.0, 1.0, 0.5]])
    rotated_normals = normals @ joint.T
    assert crystal_frame_normals(rotated_normals, joint @ grain) == pytest.approx(
        crystal_frame_normals(normals, grain)
    )
    alpha = 1.5 * np.log(2.0)
    assert cubic_attachment_weight(rotated_normals, joint @ grain, alpha) == pytest.approx(
        cubic_attachment_weight(normals, grain, alpha)
    )


def test_cube_edges_and_corners_are_not_assigned_unique_normals():
    half_width = 1.0
    points = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )
    normals, labels = cube_reference_normals(
        points, np.zeros(3), half_width, np.eye(3), edge_exclusion_m=0.1
    )
    assert labels.tolist() == ["face", "edge", "corner"]
    assert normals[0] == pytest.approx([1.0, 0.0, 0.0])
    assert np.all(np.isnan(normals[1:]))


def test_interface_points_are_physical_face_centres():
    solid, dx = _axis_plane(n=8, dx=2.0)
    points = interface_face_centres(solid, dx)
    assert points.shape == (64, 3)
    assert np.all(points[:, 0] == 0.0)
