import numpy as np
import pytest

from salt_phase2.orientation_verification import (
    central_plane_box_area,
    geometry_catalogue,
)


def test_central_axis_plane_area_equals_square_cross_section():
    assert central_plane_box_area(np.array([1.0, 0.0, 0.0]), 2.0) == pytest.approx(4.0)


def test_geometry_catalogue_volumes_and_areas_are_analytical():
    cases = geometry_catalogue()
    assert cases["plane_axis"].exact_volume_m3 == pytest.approx(0.5 * 1.2e-3**3)
    assert cases["sphere_dev"].exact_area_m2 == pytest.approx(4 * np.pi * 0.28e-3**2)
    assert cases["cube_axis"].exact_volume_m3 == pytest.approx((0.5e-3) ** 3)


def test_rotated_cube_is_fully_contained_and_volume_refines():
    cube = geometry_catalogue()["cube_rot_a"]
    errors = []
    for voxel_um in (40.0, 10.0):
        mask = cube.voxelize(voxel_um * 1.0e-6)
        volume = np.count_nonzero(mask) * (voxel_um * 1.0e-6) ** 3
        errors.append(abs(volume / cube.exact_volume_m3 - 1.0))
    assert errors[1] < errors[0]
