"""Surface-normal and cubic-orientation utilities for voxel geometries.

The signed-distance convention is ``phi < 0`` in the solid and ``phi > 0``
outside.  Consequently, ``grad(phi) / |grad(phi)|`` is the outward normal.
All smoothing and fitting lengths are supplied in physical units.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations, product

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import ndimage
from scipy.spatial import cKDTree


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


@dataclass(frozen=True)
class NormalEstimate:
    """Normals evaluated at voxel-interface face centres."""

    points_m: FloatArray
    normals: FloatArray
    valid: BoolArray
    method: str
    physical_length_m: float


def _validate_mask_and_spacing(mask: ArrayLike, voxel_size_m: float) -> BoolArray:
    solid = np.asarray(mask, dtype=bool)
    if solid.ndim != 3 or min(solid.shape) < 3:
        raise ValueError("solid_mask must be a three-dimensional array of size >= 3")
    if not np.isfinite(voxel_size_m) or voxel_size_m <= 0.0:
        raise ValueError("voxel_size_m must be positive and finite")
    if not np.any(solid) or np.all(solid):
        raise ValueError("solid_mask must contain both solid and exterior cells")
    return solid


def signed_distance_from_mask(mask: ArrayLike, voxel_size_m: float) -> FloatArray:
    """Return a cell-centred Euclidean signed-distance approximation.

    The interface lies between unlike cell centres.  The EDT magnitude is
    therefore shifted by half a voxel so that the reconstructed zero level is
    at the intervening voxel face.
    """

    solid = _validate_mask_and_spacing(mask, voxel_size_m)
    outside = ndimage.distance_transform_edt(~solid, sampling=voxel_size_m)
    inside = ndimage.distance_transform_edt(solid, sampling=voxel_size_m)
    phi = outside - inside
    phi[~solid] -= 0.5 * voxel_size_m
    phi[solid] += 0.5 * voxel_size_m
    return np.asarray(phi, dtype=float)


def interface_face_centres(
    mask: ArrayLike,
    voxel_size_m: float,
    *,
    max_points: int | None = None,
) -> FloatArray:
    """Return deterministic centres of all internal solid/exterior voxel faces."""

    solid = _validate_mask_and_spacing(mask, voxel_size_m)
    origin = -0.5 * np.asarray(solid.shape, dtype=float) * voxel_size_m
    pieces: list[FloatArray] = []
    for axis in range(3):
        left_slice = [slice(None)] * 3
        right_slice = [slice(None)] * 3
        left_slice[axis] = slice(0, -1)
        right_slice[axis] = slice(1, None)
        transitions = solid[tuple(left_slice)] != solid[tuple(right_slice)]
        indices = np.argwhere(transitions).astype(float)
        if indices.size == 0:
            continue
        points = origin + (indices + 0.5) * voxel_size_m
        points[:, axis] += 0.5 * voxel_size_m
        pieces.append(points)
    if not pieces:
        raise ValueError("mask has no internal solid/exterior interface")
    points = np.concatenate(pieces, axis=0)
    order = np.lexsort((points[:, 2], points[:, 1], points[:, 0]))
    points = points[order]
    if max_points is not None:
        if max_points < 1:
            raise ValueError("max_points must be positive")
        if points.shape[0] > max_points:
            selected = np.linspace(0, points.shape[0] - 1, max_points, dtype=int)
            points = points[selected]
    return points


def _grid_coordinates(points_m: FloatArray, shape: tuple[int, ...], dx: float) -> FloatArray:
    origin = -0.5 * np.asarray(shape, dtype=float) * dx
    return ((points_m - origin) / dx - 0.5).T


def _sample_vector_field(
    field: tuple[FloatArray, FloatArray, FloatArray],
    points_m: FloatArray,
    dx: float,
) -> FloatArray:
    coords = _grid_coordinates(points_m, field[0].shape, dx)
    return np.column_stack(
        [ndimage.map_coordinates(component, coords, order=1, mode="nearest") for component in field]
    )


def _normalise_rows(vectors: ArrayLike) -> tuple[FloatArray, BoolArray]:
    values = np.asarray(vectors, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("vectors must have shape (n, 3)")
    magnitudes = np.linalg.norm(values, axis=1)
    valid = np.isfinite(magnitudes) & (magnitudes > 1.0e-14)
    result = np.full_like(values, np.nan, dtype=float)
    result[valid] = values[valid] / magnitudes[valid, None]
    return result, valid


def estimate_sdf_normals(
    mask: ArrayLike,
    voxel_size_m: float,
    smoothing_length_m: float,
    *,
    points_m: ArrayLike | None = None,
    max_points: int | None = None,
) -> NormalEstimate:
    """Estimate outward normals from a smoothed signed-distance gradient."""

    solid = _validate_mask_and_spacing(mask, voxel_size_m)
    if not np.isfinite(smoothing_length_m) or smoothing_length_m <= 0.0:
        raise ValueError("smoothing_length_m must be positive and finite")
    points = (
        interface_face_centres(solid, voxel_size_m, max_points=max_points)
        if points_m is None
        else np.asarray(points_m, dtype=float)
    )
    phi = signed_distance_from_mask(solid, voxel_size_m)
    sigma_voxels = smoothing_length_m / voxel_size_m
    smoothed = ndimage.gaussian_filter(phi, sigma=sigma_voxels, mode="nearest")
    gradient = tuple(np.gradient(smoothed, voxel_size_m, edge_order=2))
    sampled = _sample_vector_field(gradient, points, voxel_size_m)
    normals, valid = _normalise_rows(sampled)
    return NormalEstimate(points, normals, valid, "signed_distance_gradient", smoothing_length_m)


def estimate_plane_fit_normals(
    mask: ArrayLike,
    voxel_size_m: float,
    fitting_radius_m: float,
    *,
    points_m: ArrayLike | None = None,
    max_points: int | None = None,
) -> NormalEstimate:
    """Estimate outward normals with a Gaussian-weighted local PCA plane fit.

    The fit uses all voxel-interface face centres within ``fitting_radius_m``.
    A raw signed-distance gradient is used only to resolve the PCA sign.
    """

    solid = _validate_mask_and_spacing(mask, voxel_size_m)
    if not np.isfinite(fitting_radius_m) or fitting_radius_m <= voxel_size_m:
        raise ValueError("fitting_radius_m must be finite and exceed one voxel")
    full_points = interface_face_centres(solid, voxel_size_m)
    query_points = (
        interface_face_centres(solid, voxel_size_m, max_points=max_points)
        if points_m is None
        else np.asarray(points_m, dtype=float)
    )
    tree = cKDTree(full_points)
    neighbour_lists = tree.query_ball_point(query_points, fitting_radius_m)
    phi = signed_distance_from_mask(solid, voxel_size_m)
    raw_gradient = tuple(np.gradient(phi, voxel_size_m, edge_order=2))
    orientation = _sample_vector_field(raw_gradient, query_points, voxel_size_m)
    orientation, orientation_valid = _normalise_rows(orientation)
    normals = np.full_like(query_points, np.nan, dtype=float)
    valid = np.zeros(query_points.shape[0], dtype=bool)
    sigma = 0.5 * fitting_radius_m
    for index, neighbours in enumerate(neighbour_lists):
        if len(neighbours) < 6 or not orientation_valid[index]:
            continue
        local = full_points[np.asarray(neighbours, dtype=int)]
        offsets = local - query_points[index]
        distances = np.linalg.norm(offsets, axis=1)
        weights = np.exp(-0.5 * (distances / sigma) ** 2)
        weight_sum = float(np.sum(weights))
        if weight_sum <= 0.0:
            continue
        centroid = np.sum(local * weights[:, None], axis=0) / weight_sum
        centred = local - centroid
        covariance = (centred * weights[:, None]).T @ centred / weight_sum
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        if not np.all(np.isfinite(eigenvalues)) or eigenvalues[1] <= 1.0e-24:
            continue
        normal = eigenvectors[:, 0]
        if np.dot(normal, orientation[index]) < 0.0:
            normal = -normal
        normals[index] = normal / np.linalg.norm(normal)
        valid[index] = True
    return NormalEstimate(query_points, normals, valid, "weighted_plane_fit", fitting_radius_m)


def angular_errors_deg(estimated: ArrayLike, reference: ArrayLike) -> FloatArray:
    """Return signed-direction angular errors in degrees (0 to 180)."""

    estimate, valid_estimate = _normalise_rows(estimated)
    truth, valid_truth = _normalise_rows(reference)
    valid = valid_estimate & valid_truth
    result = np.full(estimate.shape[0], np.nan, dtype=float)
    dot = np.einsum("ij,ij->i", estimate[valid], truth[valid])
    result[valid] = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))
    return result


def rotation_matrix_xyz_degrees(x_deg: float, y_deg: float, z_deg: float) -> FloatArray:
    """Return active grain-to-laboratory rotation ``Rz @ Ry @ Rx``."""

    x, y, z = np.radians([x_deg, y_deg, z_deg])
    rx = np.array([[1, 0, 0], [0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]])
    ry = np.array([[np.cos(y), 0, np.sin(y)], [0, 1, 0], [-np.sin(y), 0, np.cos(y)]])
    rz = np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0], [0, 0, 1]])
    return np.asarray(rz @ ry @ rx, dtype=float)


def _validate_rotation(rotation: ArrayLike) -> FloatArray:
    matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3)")
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1.0e-12):
        raise ValueError("rotation must be orthogonal")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1.0e-12):
        raise ValueError("rotation must be proper (determinant +1)")
    return matrix


def crystal_frame_normals(normals_lab: ArrayLike, grain_to_lab: ArrayLike) -> FloatArray:
    """Transform lab normals using ``n_c = R_g.T @ n_lab``."""

    normals, valid = _normalise_rows(normals_lab)
    rotation = _validate_rotation(grain_to_lab)
    result = np.full_like(normals, np.nan)
    result[valid] = normals[valid] @ rotation
    return result


def cubic_q4(normals_crystal: ArrayLike) -> FloatArray:
    normals, valid = _normalise_rows(normals_crystal)
    result = np.full(normals.shape[0], np.nan, dtype=float)
    result[valid] = np.sum(normals[valid] ** 4, axis=1)
    return result


def cubic_attachment_weight(
    normals_lab: ArrayLike,
    grain_to_lab: ArrayLike,
    alpha: float,
) -> FloatArray:
    """Evaluate ``exp(alpha * (q4 - 1/3))`` without kinetic interpretation."""

    if not np.isfinite(alpha):
        raise ValueError("alpha must be finite")
    q4 = cubic_q4(crystal_frame_normals(normals_lab, grain_to_lab))
    return np.exp(alpha * (q4 - 1.0 / 3.0))


def proper_cubic_rotations() -> tuple[FloatArray, ...]:
    """Return the 24 signed-permutation matrices of the proper cubic group."""

    rotations: list[FloatArray] = []
    for permutation in permutations(range(3)):
        base = np.eye(3)[:, permutation]
        for signs in product((-1.0, 1.0), repeat=3):
            candidate = base @ np.diag(signs)
            if np.isclose(np.linalg.det(candidate), 1.0):
                rotations.append(candidate)
    return tuple(rotations)


def cube_reference_normals(
    points_lab_m: ArrayLike,
    centre_m: ArrayLike,
    half_width_m: float,
    cube_to_lab: ArrayLike,
    edge_exclusion_m: float,
) -> tuple[FloatArray, NDArray[np.str_]]:
    """Return unique face normals and face/edge/corner classifications.

    Normals at mathematical edges and corners are deliberately returned as NaN.
    The classification is based on proximity of the two or three largest local
    absolute coordinates, with a physical exclusion width.
    """

    points = np.asarray(points_lab_m, dtype=float)
    centre = np.asarray(centre_m, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or centre.shape != (3,):
        raise ValueError("invalid point or centre shape")
    if half_width_m <= 0.0 or edge_exclusion_m < 0.0:
        raise ValueError("invalid cube length")
    rotation = _validate_rotation(cube_to_lab)
    local = (points - centre) @ rotation
    absolute = np.abs(local)
    order = np.argsort(absolute, axis=1)
    largest = absolute[np.arange(points.shape[0]), order[:, 2]]
    second = absolute[np.arange(points.shape[0]), order[:, 1]]
    third = absolute[np.arange(points.shape[0]), order[:, 0]]
    near_second = largest - second <= edge_exclusion_m
    near_third = second - third <= edge_exclusion_m
    labels = np.full(points.shape[0], "face", dtype="<U6")
    labels[near_second] = "edge"
    labels[near_second & near_third] = "corner"
    normals = np.full_like(points, np.nan, dtype=float)
    face_indices = np.flatnonzero(labels == "face")
    axes = order[face_indices, 2]
    signs = np.sign(local[face_indices, axes])
    signs[signs == 0.0] = 1.0
    local_normals = np.zeros((face_indices.size, 3), dtype=float)
    local_normals[np.arange(face_indices.size), axes] = signs
    normals[face_indices] = local_normals @ rotation.T
    return normals, labels


def orientation_corrected_surface_area(
    estimated_normals: ArrayLike,
    voxel_size_m: float,
) -> float:
    """Correct digital face area by the local L1 projection factor."""

    normals, valid = _normalise_rows(estimated_normals)
    if not np.any(valid):
        return float("nan")
    factors = np.sum(np.abs(normals[valid]), axis=1)
    return float(np.sum(voxel_size_m**2 / factors))
