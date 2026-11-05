import numpy as np

from prat_growth.grid import top_refined_grid
from prat_growth.material_tracking import reconstruct, synthetic_profile


def test_synthetic_fixed_material_birth_and_reaction_integral():
    h0, velocity, eps_birth = 1e-3, 3e-8, .43
    amplitude, length, birth = .19, .5e-6, 3600.
    age = np.linspace(0, 300, 30001)
    result = synthetic_profile(h0, velocity, eps_birth, amplitude, length, birth, age)
    assert result["z0"] == h0+velocity*birth
    np.testing.assert_allclose(result["distance"], velocity*age)
    np.testing.assert_allclose(result["epsilon_reaction"], result["epsilon"], atol=2e-10)
    # A fixed ALE index x follows z=x*h(t), not the fixed material z0.
    x = result["z0"]/(h0+velocity*birth)
    later_z_at_fixed_x = x*(h0+velocity*(birth+age[-1]))
    assert later_z_at_fixed_x != result["z0"]


def test_reconstruction_is_bounded_and_top_uses_birth_state():
    faces, centers, _ = top_refined_grid(8, 2)
    values = np.linspace(.25, .40, 8)
    h = .0012
    assert reconstruct(h, h, faces, centers, values, .43) == .43
    for z in np.linspace(0, h, 100):
        value = reconstruct(z, h, faces, centers, values, .43, "linear")
        assert .25 <= value <= .43
    assert reconstruct(.5*h, h, faces, centers, values, .43, "piecewise_constant") in values


def test_archived_m4_diagnostic_targets_and_explicit_limitation():
    import csv
    from pathlib import Path
    out = Path(__file__).resolve().parents[1]/"outputs/m4"
    if not (out/"cohort_summary.csv").exists():
        import pytest
        pytest.skip("run python -m prat_growth.m4_analysis")
    with (out/"cohort_summary.csv").open() as stream:
        cohorts = list(csv.DictReader(stream))
    newborn = [r for r in cohorts if r["marker_type"] == "newborn"]
    assert newborn
    assert all(float(r["cross_grid_final_epsilon_difference"]) <= .001 for r in newborn)
    assert all(float(r["output_halving_reaction_difference_relative"]) <= .01 for r in newborn)
    # Preserve the reported limitation: pointwise reconstructed balance exceeds
    # the 1% target transiently, although its final error is below 1%.
    assert any(float(r["relative_balance_error_to_drop"]) > .01 for r in newborn)
    assert all(float(r["relative_balance_error_to_drop"]) < .022 for r in newborn)
    assert all(float(r["final_balance_error_relative_to_drop"]) < .01 for r in newborn)
    for closure in ("B1", "B2"):
        fine = [float(r["final_epsilon"]) for r in newborn
                if r["closure"] == closure and r["n_cells"] == "640"]
        assert max(fine)-min(fine) < .001
