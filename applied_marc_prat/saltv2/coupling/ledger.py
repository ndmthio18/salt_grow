"""Water and salt mass ledger (plan/17 P4 item 5).

The point of this module is that every kilogram of salt and water that the
closed loop moves is recorded in exactly one named channel, so that the run
can state a conservation residual instead of asserting conservation.  It
replaces the single global ``m_front`` scalar of the E1 pre-experiment
(plan/17 section 2.2 item 7 and the P4 "forbidden practices" list).

Channels (all cumulative, all in kg, all non-negative unless stated):

salt
    ``salt_grown``     solid salt put into NEW voxels by the growth law,
                       ``sum_f rho_cr (1 - eps_s) V_n A_f dt`` (section 4.2)
    ``salt_front``     salt precipitated into the pores of EXISTING voxels at
                       internal wet/dry fronts (section 4.3); it comes out of
                       brine already inside the structure
    ``salt_front_external``  salt precipitated into the pores of EXISTING wet
                       voxels by their FROZEN external faces (plan/16 M2, section
                       4.3): the face keeps evaporating, cannot add a new voxel,
                       so its salt fills the voxel behind it.  Unlike
                       ``salt_front`` this salt is drawn from the substrate pores
                       in the same step, so it belongs on the network balance
    ``salt_released_by_evaporation``  ``beta`` times the water that left through
                       WET SALT interfaces (external faces and internal wet/dry
                       fronts), ``beta = C/(1 - C)`` of the brine that feeds the face:
                       the feeding pores' in the 3-D driver (since 2026-09-19),
                       ``C_sat`` in the 2-D one.  This is what evaporation hands
                       over and it is measured from the vapour field, not from
                       the growth law, so the balance
                       ``salt_grown + salt_front + salt_front_external (+ the
                       crust channels) = salt_released_by_evaporation (+ base)``
                       is an INDEPENDENT check on the growth-law bookkeeping.
                       Before 2026-09-06 the driver placed only the ``F`` share
                       and this balance failed by ``1 - F`` (5 % at 1 mm beads,
                       24 % at 34 um); since task B the advancing face's
                       ``(1 - F)`` share is carried by the birth porosity
    ``salt_nucleated`` solid salt of the localized seed voxel a nucleation event
                       places.  The seed geometry is PRESCRIBED, so this mass is not taken out of the pore
                       that nucleated; it is a named gap, not a hidden one
    ``salt_pnm_precipitate``  salt that precipitated INSIDE the pores when a driver runs with
                       ``saltv2.pnm.precipitation`` (subflorescence); it leaves the dissolved phase,
                       takes liquid volume and narrows the throats.  Zero in the production path,
                       where the pores are only capped
    ``salt_pnm_cap``   salt removed from the network by the saturation cap
                       ``C > C_sat -> C = C_sat`` in the substrate pores
    ``salt_supplied``  salt that leaves the substrate pore network for the
                       structure, as prescribed by the growth law
    ``salt_brine_in_pores``  the part of it that is still dissolved, sitting in
                       the pore space of the new voxels (it becomes
                       ``salt_front`` when an internal front dries that voxel)
water
    ``water_evaporated``       all water leaving through the gas phase
    ``water_evaporated_base``  the share of it that left through BARE substrate
                               openings.  Its salt stays in the pore layer and
                               raises the local concentration (it is what drives
                               nucleation), so it must NOT appear in the
                               precipitation balance above
    ``water_evaporated_on_salt``  the share that left through WET SALT faces,
                               measured on exactly the faces and with exactly the
                               ``max(j, 0)`` the salt side uses.  It is NOT
                               ``water_evaporated - water_evaporated_base``: the
                               two above are signed sums, and a salt face next to
                               a still-open pore can see a negative flux (the
                               opening's fresher brine has the higher ``a_w``, so
                               its plume condenses there).  The growth closure has
                               no dissolution branch and drops those faces, so
                               only this channel can be compared with
                               ``salt_released_by_evaporation`` as an identity
geometry
    ``volume_unplaced``   solid volume the geometry engine could not place
                          (``GeometryEngine.bank``), in m^ndim, not kg

The residuals are formed by the driver through ``balance``, because which
channels must cancel is a property of the coupled model, not of the
bookkeeping.  ``check`` raises when a registered balance drifts beyond its
tolerance, so a silent conservation failure cannot survive a long run.
"""
from __future__ import annotations

from dataclasses import dataclass, field

SALT_CHANNELS = ("salt_grown", "salt_front", "salt_front_external",
                 "salt_released_by_evaporation", "salt_pnm_cap", "salt_pnm_precipitate",
                 "salt_supplied",
                 "salt_brine_in_pores", "salt_nucleated",
                 # the bare substrate evaporates too, and the brine it draws carries salt
                 # (2026-09-07).  ``salt_released_by_base`` is all of it; ``salt_base_uncovered``
                 # is the part this run did NOT place -- either because the column carries no salt
                 # (a model with prescribed seeds and no nucleation has nowhere to put it) or
                 # because the caller switched the base-salt path off as an ablation.  It is booked,
                 # never dropped silently.  The remainder precipitates into the salt voxel sitting
                 # on that column and appears in ``salt_front``.
                 "salt_released_by_base", "salt_base_uncovered",
                 # 3-D driver with a conducting crust (2026-09-11).  ``salt_crust_front`` is the
                 # salt a wet face evaporating INTO the crust puts into the sealed voxel in front
                 # of it (front tracking: the front finishes a voxel before the band behind it
                 # starts); ``salt_base`` is the base-flux salt placed into a sealed z = 0 voxel
                 # (in 3-D a bare opening's salt stays in its pore, so only the sealed case is
                 # "base salt"; the 2-D driver's bare columns book ``salt_base_uncovered``).
                 "salt_crust_front", "salt_base",
                 # salt placed by a wet cluster that has no substrate pore under any of its
                 # columns (it floats over inert substrate after the sealing front cut it off):
                 # nothing in the network supplied it, so it is booked here by name.
                 "salt_unsourced")
# 3-D driver only (2026-09-19), kept out of DEFAULT_CHANNELS so that the 2-D records keep their
# channel list.  ``salt_crystallised``: salt that left a covered pore's brine because it exceeded
# C_sat (and, at nucleation, the supersaturation the bare pore had built up) -- the structure's
# ONLY salt supply.  ``salt_interface_buffer``: the part of it not placed yet at the end of the
# run; SIGNED (a seed may hold more salt than its pore gave).  The balance is
# crystallised = supplied + nucleated + buffer.
INTERFACE_CHANNELS = ("salt_crystallised", "salt_interface_buffer")
# 3-D driver only (2026-09-21): the LIQUID mass of the pore network [kg].  ``liquid_in`` crosses
# the reservoir throats, ``liquid_out`` is what the pores lose (evaporated water, brine drawn by
# the structure, liquid displaced by in-pore precipitate), ``liquid_stored`` is what they gain
# through their own density, V d rho (pnm.network_step); in = out + stored is the network's
# balance.  ``liquid_nucleation_jump``: V d rho of a pore that is reset to C_sat when its first
# crystal appears -- instantaneous, at fixed volume, with no flux behind it (a named gap).
LIQUID_CHANNELS = ("liquid_in", "liquid_out", "liquid_stored", "liquid_nucleation_jump")
# ... and its dissolved SALT [kg]: ``salt_reservoir_in`` crosses the reservoir throats (net),
# ``salt_pore_inventory`` is the change of what the pores hold, written once at the end of the
# run from the two states.  in = crystallised + in-pore precipitate + cap + inventory change.
PORE_SALT_CHANNELS = ("salt_reservoir_in", "salt_pore_inventory")
WATER_CHANNELS = ("water_evaporated", "water_evaporated_base", "water_evaporated_on_salt")
GEOMETRY_CHANNELS = ("volume_unplaced",)
DEFAULT_CHANNELS = SALT_CHANNELS + WATER_CHANNELS + GEOMETRY_CHANNELS


@dataclass
class Balance:
    """``sum(plus) - sum(minus)`` must stay within ``atol + rtol * scale``."""

    name: str
    plus: tuple[str, ...]
    minus: tuple[str, ...]
    rtol: float = 1e-10
    atol: float = 0.0


@dataclass
class MassLedger:
    """Named cumulative totals plus registered balances."""

    totals: dict[str, float] = field(default_factory=lambda: {c: 0.0 for c in DEFAULT_CHANNELS})
    balances: dict[str, Balance] = field(default_factory=dict)

    # ------------------------------------------------------------------
    def add(self, channel: str, value: float) -> None:
        """Accumulate into ``channel``.  Unknown names are rejected: a typo
        must not create a silent second account."""
        if channel not in self.totals:
            raise KeyError(f"unknown ledger channel {channel!r}; "
                           f"known: {sorted(self.totals)}")
        v = float(value)
        if v != v:
            raise ValueError(f"non-finite contribution to {channel!r}")
        self.totals[channel] += v

    def __getitem__(self, channel: str) -> float:
        return self.totals[channel]

    # ------------------------------------------------------------------
    def register(self, name: str, plus, minus, rtol: float = 1e-10,
                 atol: float = 0.0) -> None:
        plus, minus = tuple(plus), tuple(minus)
        for c in plus + minus:
            if c not in self.totals:
                raise KeyError(f"unknown ledger channel {c!r}")
        self.balances[name] = Balance(name, plus, minus, rtol, atol)

    def residual(self, name: str) -> float:
        """Absolute imbalance [kg] of a registered balance."""
        b = self.balances[name]
        return (sum(self.totals[c] for c in b.plus)
                - sum(self.totals[c] for c in b.minus))

    def relative_residual(self, name: str) -> float:
        """Imbalance divided by the largest side of the balance (0 if empty)."""
        b = self.balances[name]
        scale = max(abs(sum(self.totals[c] for c in b.plus)),
                    abs(sum(self.totals[c] for c in b.minus)))
        if scale == 0.0:
            return 0.0
        return self.residual(name) / scale

    def check(self, name: str | None = None) -> None:
        """Raise if a registered balance has drifted beyond its tolerance."""
        names = list(self.balances) if name is None else [name]
        for n in names:
            b = self.balances[n]
            scale = max(abs(sum(self.totals[c] for c in b.plus)),
                        abs(sum(self.totals[c] for c in b.minus)))
            if abs(self.residual(n)) > b.atol + b.rtol * scale:
                raise RuntimeError(
                    f"mass balance {n!r} violated: residual {self.residual(n):.6e} kg, "
                    f"relative {self.relative_residual(n):.3e}, "
                    f"tolerance {b.atol + b.rtol * scale:.3e}")

    # ------------------------------------------------------------------
    def report(self) -> dict:
        """Flat dict for the run record: every channel and every residual."""
        out = dict(self.totals)
        for n in self.balances:
            out[f"residual_{n}"] = self.residual(n)
            out[f"relative_residual_{n}"] = self.relative_residual(n)
        return out
