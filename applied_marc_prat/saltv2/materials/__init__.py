"""Material properties: the only place in ``saltv2`` where a property VALUE is written down.

    constants   molar masses, gas constant, gravity
    brine       NaCl brine: density, molar conversions, diffusivity, water activity (DuMuX port),
                viscosity, surface tension, saturation
    air         moist air: vapour diffusivity, gas density, saturated vapour mass fraction
    salt        the solid and the porous structure it builds: crystal density, birth porosity

Everything else in the package imports the value from here; a default argument elsewhere reads
``= materials.<module>.<NAME>`` and never repeats the number.  Derived quantities are computed
here once and not stored a second time.  This layer imports nothing from the rest of the package.

Provenance tags on every value: [M] measured, [L] literature / library, [A] assumption,
[B] back-calibrated, [V] validated against a reference, [N] numerical.
"""
from . import air, brine, constants, salt

__all__ = ["air", "brine", "constants", "salt"]
