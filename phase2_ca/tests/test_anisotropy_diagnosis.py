import numpy as np
import pytest

from salt_phase2 import (
    equivalent_displacement_m,
    group_allocation,
    support_advances_m,
    uniform_cube_facet_displacement_m,
    volume_fraction_for_facet_displacement,
)


def test_uniform_cube_displacement_and_inverse_for_one_voxel_layer():
    fraction = volume_fraction_for_facet_displacement(160.0e-6, 10.0e-6)
    assert fraction == pytest.approx(0.423828125)
    assert uniform_cube_facet_displacement_m(160.0e-6, fraction) == pytest.approx(
        10.0e-6
    )


def test_equivalent_displacement_is_deposited_volume_over_area():
    assert equivalent_displacement_m(100, 2.0, 50.0) == pytest.approx(16.0)


def test_site_and_area_allocation_expose_candidate_density_dependence():
    counts = np.array([100.0, 200.0])
    areas = np.array([1.0, 1.0])
    weights = np.array([2.0, 1.0])
    assert group_allocation(counts, areas, weights, allocation="site_weighted") == pytest.approx(
        [0.5, 0.5]
    )
    assert group_allocation(counts, areas, weights, allocation="area_weighted") == pytest.approx(
        [2.0 / 3.0, 1.0 / 3.0]
    )


def test_support_advance_uses_growth_only_and_crystal_directions():
    initial = np.zeros((7, 7, 7), dtype=bool)
    initial[2:5, 2:5, 2:5] = True
    final = initial.copy()
    final[5, 3, 3] = True
    advances = support_advances_m(
        initial,
        final,
        voxel_size_m=2.0,
        grain_to_lab=np.eye(3),
        directions_crystal=np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]),
    )
    assert advances == pytest.approx([2.0, 0.0])
