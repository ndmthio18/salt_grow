import numpy as np
import pytest

from salt_phase2 import (
    CrystallographicCandidate,
    air_exposed_face_area_m2,
    classify_crystal_directions,
    directional_growth_extents,
    fixed_front_probabilities,
    sample_fixed_front,
    symmetric_chamfer_distance,
)


def _record(weight):
    return CrystallographicCandidate(
        coordinate=(0, 0, 0),
        site_class="face_like",
        adjacent_salt_faces=1,
        normal_lab=(1.0, 0.0, 0.0),
        normal_crystal=(1.0, 0.0, 0.0),
        q4=1.0,
        orientation_family="<100>",
        family_angular_error_deg=0.0,
        represented_area_m2=1.0,
        normal_quality=1.0,
        anisotropy_eligible=True,
        anisotropy_weight=weight,
    )


def test_ideal_direction_classification_and_other_bin():
    directions = np.array(
        [
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.35, 0.2],
        ]
    )
    labels, errors = classify_crystal_directions(directions, near_tolerance_deg=10.0)
    assert labels.tolist() == ["<100>", "<110>", "<111>", "other"]
    assert errors[:3] == pytest.approx(0.0, abs=2.0e-6)
    assert errors[3] > 10.0


def test_invalid_direction_is_retained_as_unresolved_other():
    directions = np.array([[np.nan, np.nan, np.nan], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    labels, errors = classify_crystal_directions(directions, near_tolerance_deg=10.0)
    assert labels.tolist() == ["other", "other", "<100>"]
    assert np.isnan(errors[:2]).all()
    assert errors[2] == pytest.approx(0.0)


def test_fixed_front_multinomial_matches_normalized_weights():
    probabilities = fixed_front_probabilities([_record(1.0), _record(2.0), _record(3.0)])
    assert probabilities == pytest.approx([1 / 6, 2 / 6, 3 / 6])
    counts = sample_fixed_front(probabilities, sample_count=200_000, rng_seed=91)
    empirical = counts / np.sum(counts)
    standard_error = np.sqrt(probabilities * (1.0 - probabilities) / np.sum(counts))
    assert np.all(np.abs(empirical - probabilities) <= 4.0 * standard_error)


def test_air_surface_area_excludes_domain_exterior_and_nonempty_neighbors():
    salt = np.zeros((5, 5, 5), dtype=bool)
    salt[2, 2, 2] = True
    empty = ~salt
    assert air_exposed_face_area_m2(salt, empty, 2.0) == pytest.approx(6 * 4.0)
    empty[2, 2, 3] = False
    assert air_exposed_face_area_m2(salt, empty, 2.0) == pytest.approx(5 * 4.0)


def test_directional_extents_use_growth_increment_in_crystal_frame():
    initial = np.zeros((7, 7, 7), dtype=bool)
    initial[2:5, 2:5, 2:5] = True
    final = initial.copy()
    final[1, 3, 3] = True
    final[5, 3, 3] = True
    result = directional_growth_extents(
        initial,
        final,
        voxel_size_m=2.0,
        grain_to_lab=np.eye(3),
    )
    assert result["growth_extent_x_negative_m"] == pytest.approx(2.0)
    assert result["growth_extent_x_positive_m"] == pytest.approx(2.0)
    assert result["growth_extent_y_total_m"] == pytest.approx(0.0)


def test_symmetric_chamfer_is_zero_for_permuted_identical_points():
    first = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    second = first[[2, 0, 1]]
    assert symmetric_chamfer_distance(first, second) == pytest.approx(0.0)
