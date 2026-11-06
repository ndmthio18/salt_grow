"""Analytical geometries and metrics for Phase 2A.3a verification."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from numpy.typing import NDArray

from .surface_orientation import (
    NormalEstimate,
    angular_errors_deg,
    cube_reference_normals,
    estimate_plane_fit_normals,
    estimate_sdf_normals,
    interface_face_centres,
    rotation_matrix_xyz_degrees,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ReferenceGeometry:
    name: str
    kind: str
    domain_length_m: float
    normal: tuple[float, float, float] | None = None
    centre_m: tuple[float, float, float] = (0.0, 0.0, 0.0)
    radius_m: float | None = None
    half_width_m: float | None = None
    rotation_xyz_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @property
    def rotation(self) -> FloatArray:
        return rotation_matrix_xyz_degrees(*self.rotation_xyz_deg)

    @property
    def exact_volume_m3(self) -> float:
        if self.kind == "plane":
            return 0.5 * self.domain_length_m**3
        if self.kind == "sphere":
            assert self.radius_m is not None
            return 4.0 * np.pi * self.radius_m**3 / 3.0
        if self.kind == "cube":
            assert self.half_width_m is not None
            return (2.0 * self.half_width_m) ** 3
        raise ValueError(self.kind)

    @property
    def exact_area_m2(self) -> float:
        if self.kind == "plane":
            assert self.normal is not None
            return central_plane_box_area(np.asarray(self.normal), self.domain_length_m)
        if self.kind == "sphere":
            assert self.radius_m is not None
            return 4.0 * np.pi * self.radius_m**2
        if self.kind == "cube":
            assert self.half_width_m is not None
            return 24.0 * self.half_width_m**2
        raise ValueError(self.kind)

    def voxelize(self, voxel_size_m: float) -> NDArray[np.bool_]:
        cells = int(round(self.domain_length_m / voxel_size_m))
        if not np.isclose(cells * voxel_size_m, self.domain_length_m, atol=1.0e-15):
            raise ValueError("domain length must be an integer number of voxels")
        coordinates = (-0.5 * cells + 0.5 + np.arange(cells)) * voxel_size_m
        x, y, z = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij", sparse=True)
        centre = np.asarray(self.centre_m, dtype=float)
        if self.kind == "plane":
            normal = np.asarray(self.normal, dtype=float)
            normal /= np.linalg.norm(normal)
            phi = normal[0] * x + normal[1] * y + normal[2] * z
            return np.asarray(phi <= 0.0, dtype=bool)
        if self.kind == "sphere":
            assert self.radius_m is not None
            radius_squared = (
                (x - centre[0]) ** 2 + (y - centre[1]) ** 2 + (z - centre[2]) ** 2
            )
            return np.asarray(radius_squared <= self.radius_m**2, dtype=bool)
        if self.kind == "cube":
            assert self.half_width_m is not None
            rotation = self.rotation
            offsets = (x - centre[0], y - centre[1], z - centre[2])
            local = [sum(offsets[j] * rotation[j, i] for j in range(3)) for i in range(3)]
            return np.asarray(
                np.maximum.reduce([np.abs(component) for component in local])
                <= self.half_width_m,
                dtype=bool,
            )
        raise ValueError(self.kind)

    def reference_normals(
        self, points_m: FloatArray, edge_exclusion_m: float
    ) -> tuple[FloatArray, NDArray[np.str_]]:
        if self.kind == "plane":
            normal = np.asarray(self.normal, dtype=float)
            normal /= np.linalg.norm(normal)
            return np.broadcast_to(normal, points_m.shape).copy(), np.full(
                points_m.shape[0], "smooth", dtype="<U6"
            )
        if self.kind == "sphere":
            offsets = points_m - np.asarray(self.centre_m, dtype=float)
            offsets /= np.linalg.norm(offsets, axis=1)[:, None]
            return offsets, np.full(points_m.shape[0], "smooth", dtype="<U6")
        if self.kind == "cube":
            assert self.half_width_m is not None
            return cube_reference_normals(
                points_m,
                np.asarray(self.centre_m),
                self.half_width_m,
                self.rotation,
                edge_exclusion_m,
            )
        raise ValueError(self.kind)


def central_plane_box_area(normal: FloatArray, side_length_m: float) -> float:
    """Area of a central plane section through an axis-aligned cube."""

    normal = np.asarray(normal, dtype=float)
    normal /= np.linalg.norm(normal)
    half = 0.5 * side_length_m
    vertices = np.array(
        [[x, y, z] for x in (-half, half) for y in (-half, half) for z in (-half, half)]
    )
    edges = []
    for index, first in enumerate(vertices):
        for second in vertices[index + 1 :]:
            if np.count_nonzero(first != second) == 1:
                edges.append((first, second))
    intersections: list[FloatArray] = []
    for first, second in edges:
        d_first = float(np.dot(normal, first))
        d_second = float(np.dot(normal, second))
        if d_first * d_second > 0.0:
            continue
        denominator = d_first - d_second
        if abs(denominator) < 1.0e-15:
            continue
        fraction = d_first / denominator
        if -1.0e-12 <= fraction <= 1.0 + 1.0e-12:
            intersections.append(first + fraction * (second - first))
    points = np.unique(np.round(np.asarray(intersections), decimals=15), axis=0)
    if points.shape[0] < 3:
        raise RuntimeError("plane/cube intersection is degenerate")
    helper = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(helper, normal)) > 0.8:
        helper = np.array([0.0, 1.0, 0.0])
    axis_one = np.cross(normal, helper)
    axis_one /= np.linalg.norm(axis_one)
    axis_two = np.cross(normal, axis_one)
    projected = np.column_stack((points @ axis_one, points @ axis_two))
    centre = np.mean(projected, axis=0)
    angles = np.arctan2(projected[:, 1] - centre[1], projected[:, 0] - centre[0])
    polygon = projected[np.argsort(angles)]
    area = 0.5 * abs(
        np.dot(polygon[:, 0], np.roll(polygon[:, 1], -1))
        - np.dot(polygon[:, 1], np.roll(polygon[:, 0], -1))
    )
    return float(area)


def geometry_catalogue() -> dict[str, ReferenceGeometry]:
    length = 1.2e-3
    return {
        "plane_axis": ReferenceGeometry("plane_axis", "plane", length, normal=(1, 0, 0)),
        "plane_17deg": ReferenceGeometry(
            "plane_17deg", "plane", length, normal=(np.cos(np.radians(17)), np.sin(np.radians(17)), 0)
        ),
        "plane_33deg": ReferenceGeometry(
            "plane_33deg", "plane", length, normal=(np.cos(np.radians(33)), np.sin(np.radians(33)), 0)
        ),
        "sphere_dev": ReferenceGeometry(
            "sphere_dev", "sphere", length, radius_m=0.28e-3
        ),
        "sphere_locked": ReferenceGeometry(
            "sphere_locked",
            "sphere",
            length,
            centre_m=(37e-6, -21e-6, 15e-6),
            radius_m=0.27e-3,
        ),
        "cube_axis": ReferenceGeometry(
            "cube_axis", "cube", length, half_width_m=0.25e-3
        ),
        "cube_rot_a": ReferenceGeometry(
            "cube_rot_a",
            "cube",
            length,
            half_width_m=0.25e-3,
            rotation_xyz_deg=(17.0, 23.0, 11.0),
        ),
        "cube_rot_b": ReferenceGeometry(
            "cube_rot_b",
            "cube",
            length,
            half_width_m=0.25e-3,
            rotation_xyz_deg=(-19.0, 31.0, 7.0),
        ),
    }


def evaluate_normal_case(
    geometry: ReferenceGeometry,
    voxel_size_m: float,
    method: str,
    physical_length_m: float,
    *,
    max_points: int = 4000,
) -> tuple[dict[str, float | int | str], NormalEstimate, NDArray[np.str_]]:
    started = perf_counter()
    mask = geometry.voxelize(voxel_size_m)
    full_point_count = interface_face_centres(mask, voxel_size_m).shape[0]
    if method == "signed_distance_gradient":
        estimate = estimate_sdf_normals(
            mask,
            voxel_size_m,
            physical_length_m,
            max_points=max_points,
        )
    elif method == "weighted_plane_fit":
        estimate = estimate_plane_fit_normals(
            mask,
            voxel_size_m,
            physical_length_m,
            max_points=max_points,
        )
    else:
        raise ValueError(f"unknown method {method}")
    edge_exclusion = max(2.0 * voxel_size_m, physical_length_m)
    reference, classes = geometry.reference_normals(estimate.points_m, edge_exclusion)
    errors = angular_errors_deg(estimate.normals, reference)
    unique = np.isin(classes, ["smooth", "face"]) & estimate.valid & np.isfinite(errors)
    if not np.any(unique):
        raise RuntimeError(f"no unique-normal samples for {geometry.name}")
    inverse_l1 = 1.0 / np.sum(np.abs(estimate.normals[estimate.valid]), axis=1)
    corrected_area = full_point_count * voxel_size_m**2 * float(np.mean(inverse_l1))
    digital_area = full_point_count * voxel_size_m**2
    volume = float(np.count_nonzero(mask) * voxel_size_m**3)
    values = errors[unique]
    record: dict[str, float | int | str] = {
        "case": geometry.name,
        "kind": geometry.kind,
        "method": method,
        "physical_length_um": physical_length_m * 1.0e6,
        "voxel_size_um": voxel_size_m * 1.0e6,
        "cells_per_domain": mask.shape[0],
        "surface_face_count": full_point_count,
        "sample_count": estimate.points_m.shape[0],
        "valid_sample_count": int(np.count_nonzero(estimate.valid)),
        "unique_normal_count": int(np.count_nonzero(unique)),
        "edge_count": int(np.count_nonzero(classes == "edge")),
        "corner_count": int(np.count_nonzero(classes == "corner")),
        "median_angular_error_deg": float(np.median(values)),
        "p90_angular_error_deg": float(np.percentile(values, 90)),
        "max_angular_error_deg": float(np.max(values)),
        "estimated_surface_area_m2": corrected_area,
        "digital_face_area_m2": digital_area,
        "reference_surface_area_m2": geometry.exact_area_m2,
        "surface_area_relative_error": abs(corrected_area / geometry.exact_area_m2 - 1.0),
        "digital_area_relative_error": abs(digital_area / geometry.exact_area_m2 - 1.0),
        "estimated_volume_m3": volume,
        "reference_volume_m3": geometry.exact_volume_m3,
        "volume_relative_error": abs(volume / geometry.exact_volume_m3 - 1.0),
        "runtime_s": perf_counter() - started,
    }
    return record, estimate, classes
