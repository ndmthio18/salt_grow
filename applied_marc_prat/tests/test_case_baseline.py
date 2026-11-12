"""The self-consistent baseline case (saltv2/cases/baseline.py).

Regime: evaporation-wicking.  The case must state everything that distinguishes it from the
archived crust chain -- the salt structure, c_k = 1, the temperature, the far field re-derived at
that temperature -- and a run of it must not freeze, must not seal and must close its books.
"""
import numpy as np
import pytest

from saltv2.cases import baseline, vtp2014_column
from saltv2.coupling import driver3d
from saltv2.coupling.driver3d import col
from saltv2.growth.crust import HIERARCHICAL, LEGACY_MIXED, SINGLE_SCALE
from saltv2.materials import air


def test_the_far_field_is_the_old_calibration_at_the_temperature_of_the_experiment():
    """Same four measured rates, same procedure: with the 20 C constants it returns the 18 mm the
    archived cases use, at 22 C it is longer by exactly the ratio of rho_g D_v W_sat."""
    legacy = air.LEGACY.rho_D * air.LEGACY.W_sat
    gas = air.moist_air(baseline.T, baseline.RH)
    at_22C = baseline.calibrated_delta_eff(baseline.T) * 1e3
    assert at_22C * legacy / (gas.rho_D * gas.W_sat) == pytest.approx(18.11, abs=0.01)   # plan/14 s.12
    assert baseline.DELTA_EFF_MM == round(at_22C, 1) == 20.6
    assert vtp2014_column.DELTA_EFF_MM == 18.0


def test_the_case_states_what_distinguishes_it(tmp_path):
    dgf = vtp2014_column.write_dgf("1mm", tmp_path)
    kw = baseline.case("1mm", dgf, SINGLE_SCALE)
    assert kw["c_k"] == 1.0 and kw["T"] == 295.15 and kw["salt_structure"] is SINGLE_SCALE
    assert kw["d_sub"] == baseline.BEADS["1mm"] and kw["delta_eff_mm"] == baseline.DELTA_EFF_MM
    with pytest.raises(TypeError):
        baseline.case("1mm", dgf)                       # the structure has no default
    assert baseline.dumux_column(HIERARCHICAL)["salt_structure"] is HIERARCHICAL
    assert LEGACY_MIXED not in (kw["salt_structure"], baseline.dumux_column(SINGLE_SCALE)["salt_structure"])


@pytest.mark.parametrize("structure", [SINGLE_SCALE, HIERARCHICAL], ids=lambda s: s.name)
def test_a_baseline_run_grows_without_freezing_or_sealing_and_closes_its_books(tmp_path, structure):
    dgf = vtp2014_column.write_dgf("1mm", tmp_path)
    kw = dict(baseline.case("1mm", dgf, structure), t_end=6 * 3600.0, verbose=False)
    res = driver3d.run(**kw)
    h = res["hist"]
    assert np.isfinite(res["t_nuc"]).all() and col(h, "h_max")[-1] > 0.0
    assert col(h, "n_frozen").max() == 0 and col(h, "n_sealed_cells").max() == 0
    led = res["ledger"]
    for balance in ("salt_pores", "salt_interface", "salt_network", "salt_precipitation",
                    "liquid_network"):
        assert abs(led[f"relative_residual_{balance}"]) < 1e-10
    assert np.abs(col(h, "res_salt")).max() < 1e-10 and np.abs(col(h, "res_water")).max() < 1e-10
    assert res["closure"]["salt_structure"] == structure.name and res["closure"]["T"] == baseline.T
