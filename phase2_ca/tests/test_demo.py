import pytest

from run_demo import generate_demo


def test_demo_is_reproducible_spatial_and_mass_conservative(tmp_path):
    first_path = tmp_path / "first" / "demo_summary.json"
    second_path = tmp_path / "second" / "demo_summary.json"
    first = generate_demo(first_path)
    second = generate_demo(second_path)
    assert first == second
    assert first_path.read_bytes() == second_path.read_bytes()
    assert (first_path.parent / "SHA256SUMS.txt").read_bytes() == (
        second_path.parent / "SHA256SUMS.txt"
    ).read_bytes()
    final = first["final_morphology"]
    assert final["new_salt_voxels"] > 0
    assert final["maximum_height_voxels_above_substrate"] > 1
    assert final["mass_closure_error_kg"] == pytest.approx(0.0, abs=1.0e-24)
