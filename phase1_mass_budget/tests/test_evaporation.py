from dataclasses import replace
import math

import pytest

from salt_phase1.evaporation import (
    EvaporationParameters,
    compute_mass_budget,
    saturation_vapor_concentration_kg_m3,
)


AREA_64_10UM = (64 * 10.0e-6) ** 2


@pytest.fixture
def reference_parameters() -> EvaporationParameters:
    return EvaporationParameters(
        temperature_k=298.15,
        water_activity=0.753,
        relative_humidity=0.30,
        mass_transfer_coefficient_m_s=0.005,
        evaporation_area_m2=AREA_64_10UM,
        timestep_s=60.0,
        salt_mass_fraction=0.264,
        solid_density_kg_m3=2160.0,
        saturation_vapor_concentration_input_kg_m3=0.0230,
    )


def test_saturation_concentration_is_approximately_23_g_per_m3_at_25c():
    concentration = saturation_vapor_concentration_kg_m3(298.15)
    assert concentration == pytest.approx(0.02302, rel=2.0e-3)


def test_positive_evaporation_below_equilibrium(reference_parameters):
    result = compute_mass_budget(reference_parameters)
    assert result.regime == "evaporation"
    assert result.evaporation_flux_kg_m2_s > 0.0
    assert result.precipitated_salt_mass_kg > 0.0


def test_zero_evaporation_at_equilibrium(reference_parameters):
    result = compute_mass_budget(
        replace(reference_parameters, relative_humidity=0.753 + 5.0e-14)
    )
    assert result.regime == "equilibrium"
    assert result.evaporation_flux_kg_m2_s == 0.0
    assert result.precipitated_salt_mass_kg == 0.0


def test_above_equilibrium_is_zero_and_flagged(reference_parameters):
    result = compute_mass_budget(replace(reference_parameters, relative_humidity=0.80))
    assert result.regime == "above_equilibrium_unmodeled"
    assert result.evaporation_flux_kg_m2_s == 0.0
    assert result.condensation_or_deliquescence_modeled is False
    assert "does not model" in result.note


def test_evaporation_monotonically_increases_as_rh_decreases(reference_parameters):
    rhs = [0.70, 0.50, 0.30]
    fluxes = [
        compute_mass_budget(replace(reference_parameters, relative_humidity=rh))
        .evaporation_flux_kg_m2_s
        for rh in rhs
    ]
    assert fluxes[0] < fluxes[1] < fluxes[2]


@pytest.mark.parametrize(
    ("field", "value"),
    [("evaporation_area_m2", 0.0), ("mass_transfer_coefficient_m_s", 0.0)],
)
def test_zero_area_or_transfer_coefficient_gives_zero_mass(
    reference_parameters, field, value
):
    result = compute_mass_budget(replace(reference_parameters, **{field: value}))
    assert result.evaporated_water_mass_kg == 0.0
    assert result.precipitated_salt_mass_kg == 0.0


def test_analytical_mass_budget_with_independent_input(reference_parameters):
    result = compute_mass_budget(reference_parameters)
    expected_flux = 0.005 * (0.753 - 0.30) * 0.0230
    expected_water = expected_flux * AREA_64_10UM * 60.0
    expected_ratio = 0.264 / (1.0 - 0.264)
    expected_salt = expected_ratio * expected_water
    assert result.evaporation_flux_kg_m2_s == pytest.approx(expected_flux, rel=1.0e-14)
    assert result.evaporated_water_mass_kg == pytest.approx(expected_water, rel=1.0e-14)
    assert result.salt_to_water_mass_ratio == pytest.approx(expected_ratio, rel=1.0e-14)
    assert result.precipitated_salt_mass_kg == pytest.approx(expected_salt, rel=1.0e-14)
    assert result.precipitated_salt_volume_m3 == pytest.approx(
        expected_salt / 2160.0, rel=1.0e-14
    )


def test_timestep_partition_does_not_change_integrated_mass(reference_parameters):
    total_s = 3600.0
    coarse = compute_mass_budget(replace(reference_parameters, timestep_s=total_s))
    fine_step = compute_mass_budget(replace(reference_parameters, timestep_s=15.0))
    fine_total = math.fsum(
        fine_step.precipitated_salt_mass_kg for _ in range(int(total_s / 15.0))
    )
    assert fine_total == pytest.approx(coarse.precipitated_salt_mass_kg, rel=2.0e-15)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("temperature_k", 260.0),
        ("water_activity", 1.01),
        ("relative_humidity", -0.01),
        ("mass_transfer_coefficient_m_s", -1.0),
        ("evaporation_area_m2", -1.0),
        ("timestep_s", -1.0),
        ("salt_mass_fraction", 0.0),
        ("salt_mass_fraction", 1.0),
        ("solid_density_kg_m3", 0.0),
        ("saturation_vapor_concentration_input_kg_m3", -1.0),
    ],
)
def test_invalid_parameters_are_rejected(reference_parameters, field, value):
    with pytest.raises(ValueError):
        changed = replace(reference_parameters, **{field: value})
        compute_mass_budget(changed)


def test_nonfinite_parameters_are_rejected(reference_parameters):
    with pytest.raises(ValueError):
        replace(reference_parameters, relative_humidity=math.nan)
