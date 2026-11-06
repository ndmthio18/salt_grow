"""Evaporation-driven NaCl precipitation budget for Phase 1.

All quantities use SI units.  The module describes a continuously replenished,
saturated brine with a stationary liquid inventory and no dissolved-salt export.
It intentionally contains no cellular-automaton or morphology logic.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


MOLAR_MASS_WATER_KG_MOL = 0.01801528
UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.31446261815324


def _require_finite(name: str, value: float) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")


def saturation_vapor_pressure_pa(temperature_k: float) -> float:
    """Return saturation vapor pressure over liquid water in Pa.

    Uses the Buck (1981) relation

    ``p_sat = 0.61121 exp((18.678 - T_c/234.5) T_c/(257.14 + T_c))`` kPa.

    This implementation is deliberately limited to 0--50 degC, which covers
    the Phase 1 isothermal reference case.  It is not a brine activity model;
    solution activity is applied separately in :func:`compute_mass_budget`.
    """

    _require_finite("temperature_k", temperature_k)
    temperature_c = temperature_k - 273.15
    if not 0.0 <= temperature_c <= 50.0:
        raise ValueError("Buck relation is restricted here to 273.15--323.15 K")
    return 611.21 * math.exp(
        (18.678 - temperature_c / 234.5)
        * (temperature_c / (257.14 + temperature_c))
    )


def saturation_vapor_concentration_kg_m3(temperature_k: float) -> float:
    """Return ideal-gas saturation water-vapor concentration in kg/m3."""

    pressure_pa = saturation_vapor_pressure_pa(temperature_k)
    return (
        pressure_pa
        * MOLAR_MASS_WATER_KG_MOL
        / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature_k)
    )


@dataclass(frozen=True)
class EvaporationParameters:
    """Configurable inputs to one Phase 1 physical time increment.

    Parameters
    ----------
    temperature_k:
        Absolute temperature [K].
    water_activity:
        Saturated-brine water activity [-], interpreted as equilibrium RH.
    relative_humidity:
        Environmental relative humidity as a fraction [-].
    mass_transfer_coefficient_m_s:
        Effective vapor mass-transfer coefficient [m/s].  It is an imposed,
        uncalibrated coefficient in Phase 1.
    evaporation_area_m2:
        Fixed effective evaporation area [m2].
    timestep_s:
        Physical accounting interval [s].
    salt_mass_fraction:
        NaCl mass fraction in supplied saturated brine [kg/kg solution].
    solid_density_kg_m3:
        Density of precipitated solid NaCl [kg/m3].
    saturation_vapor_concentration_input_kg_m3:
        Optional validated saturation concentration [kg/m3].  If omitted, the
        Buck pressure relation plus the ideal-gas law is used.
    """

    temperature_k: float
    water_activity: float
    relative_humidity: float
    mass_transfer_coefficient_m_s: float
    evaporation_area_m2: float
    timestep_s: float
    salt_mass_fraction: float
    solid_density_kg_m3: float
    saturation_vapor_concentration_input_kg_m3: float | None = None

    def __post_init__(self) -> None:
        values = {
            "temperature_k": self.temperature_k,
            "water_activity": self.water_activity,
            "relative_humidity": self.relative_humidity,
            "mass_transfer_coefficient_m_s": self.mass_transfer_coefficient_m_s,
            "evaporation_area_m2": self.evaporation_area_m2,
            "timestep_s": self.timestep_s,
            "salt_mass_fraction": self.salt_mass_fraction,
            "solid_density_kg_m3": self.solid_density_kg_m3,
        }
        for name, value in values.items():
            _require_finite(name, value)
        if not 273.15 <= self.temperature_k <= 323.15:
            raise ValueError("temperature_k must be in the Phase 1 range 273.15--323.15 K")
        if not 0.0 <= self.water_activity <= 1.0:
            raise ValueError("water_activity must be in [0, 1]")
        if not 0.0 <= self.relative_humidity <= 1.0:
            raise ValueError("relative_humidity must be in [0, 1]")
        if self.mass_transfer_coefficient_m_s < 0.0:
            raise ValueError("mass_transfer_coefficient_m_s must be non-negative")
        if self.evaporation_area_m2 < 0.0:
            raise ValueError("evaporation_area_m2 must be non-negative")
        if self.timestep_s < 0.0:
            raise ValueError("timestep_s must be non-negative")
        if not 0.0 < self.salt_mass_fraction < 1.0:
            raise ValueError("salt_mass_fraction must be strictly between 0 and 1")
        if self.solid_density_kg_m3 <= 0.0:
            raise ValueError("solid_density_kg_m3 must be positive")
        if self.saturation_vapor_concentration_input_kg_m3 is not None:
            _require_finite(
                "saturation_vapor_concentration_input_kg_m3",
                self.saturation_vapor_concentration_input_kg_m3,
            )
            if self.saturation_vapor_concentration_input_kg_m3 <= 0.0:
                raise ValueError(
                    "saturation_vapor_concentration_input_kg_m3 must be positive"
                )


@dataclass(frozen=True)
class MassBudgetResult:
    """Complete water and salt budget for one interval, in SI units."""

    saturation_vapor_concentration_kg_m3: float
    interface_vapor_concentration_kg_m3: float
    ambient_vapor_concentration_kg_m3: float
    positive_driving_concentration_kg_m3: float
    evaporation_flux_kg_m2_s: float
    water_evaporation_rate_kg_s: float
    evaporated_water_mass_kg: float
    salt_to_water_mass_ratio: float
    salt_precipitation_rate_kg_s: float
    precipitated_salt_mass_kg: float
    precipitated_salt_volume_m3: float
    regime: str
    condensation_or_deliquescence_modeled: bool
    note: str


def compute_mass_budget(
    parameters: EvaporationParameters,
    *,
    equilibrium_tolerance: float = 1.0e-12,
) -> MassBudgetResult:
    """Compute an evaporation-driven steady salt precipitation budget.

    Positive evaporation is clipped at zero.  If environmental RH exceeds the
    solution water activity, the result is zero and explicitly marked as an
    unmodeled condensation/deliquescence regime.
    """

    _require_finite("equilibrium_tolerance", equilibrium_tolerance)
    if equilibrium_tolerance < 0.0:
        raise ValueError("equilibrium_tolerance must be non-negative")

    p = parameters
    c_sat = (
        p.saturation_vapor_concentration_input_kg_m3
        if p.saturation_vapor_concentration_input_kg_m3 is not None
        else saturation_vapor_concentration_kg_m3(p.temperature_k)
    )
    c_interface = p.water_activity * c_sat
    c_ambient = p.relative_humidity * c_sat
    activity_gap = p.water_activity - p.relative_humidity

    if abs(activity_gap) <= equilibrium_tolerance:
        driving = 0.0
        regime = "equilibrium"
        note = "Environmental RH equals brine water activity within tolerance."
    elif activity_gap < 0.0:
        driving = 0.0
        regime = "above_equilibrium_unmodeled"
        note = (
            "Environmental RH exceeds brine water activity; Phase 1 returns zero "
            "positive evaporation and does not model condensation, deliquescence, "
            "or dissolution."
        )
    else:
        driving = c_interface - c_ambient
        regime = "evaporation"
        note = "Positive evaporation under the imposed Phase 1 assumptions."

    flux = p.mass_transfer_coefficient_m_s * driving
    water_rate = flux * p.evaporation_area_m2
    water_mass = water_rate * p.timestep_s
    salt_to_water_ratio = p.salt_mass_fraction / (1.0 - p.salt_mass_fraction)
    salt_rate = salt_to_water_ratio * water_rate
    salt_mass = salt_to_water_ratio * water_mass
    salt_volume = salt_mass / p.solid_density_kg_m3

    return MassBudgetResult(
        saturation_vapor_concentration_kg_m3=c_sat,
        interface_vapor_concentration_kg_m3=c_interface,
        ambient_vapor_concentration_kg_m3=c_ambient,
        positive_driving_concentration_kg_m3=driving,
        evaporation_flux_kg_m2_s=flux,
        water_evaporation_rate_kg_s=water_rate,
        evaporated_water_mass_kg=water_mass,
        salt_to_water_mass_ratio=salt_to_water_ratio,
        salt_precipitation_rate_kg_s=salt_rate,
        precipitated_salt_mass_kg=salt_mass,
        precipitated_salt_volume_m3=salt_volume,
        regime=regime,
        condensation_or_deliquescence_modeled=False,
        note=note,
    )
