"""Salt precipitating INSIDE the substrate pores, ported from the DuMuX mineralisation model.

Regime note: the production 3-D loop crystallises only OUTSIDE the substrate (the voxel domain),
and the pores are merely capped at ``C_sat`` (channel ``salt_pnm_cap``).  This module is the
opt-in alternative: the salt above saturation becomes a solid phase in the pore, which takes
away liquid volume and narrows the throats, so the network feels its own clogging
(subflorescence).  Nothing here is used unless a driver asks for it.

Ported from /home/zhixin/Dumux/dumux_theresa/dumux/pnm2pncmin2, application
2pncmin_isothermal (``spatialparams.hh``, ``problem.hh``), equilibrium branch
(``UseKineticPrecipitationReaction = false``), throat update "average source terms"
(``ThroatApproach = 2``):

    equilibrium      everything above C_sat precipitates at once, no kinetics [L]
    pore volume      V = (1 - phi_s) V_0,  phi_s = solid volume / initial pore volume
    throat area      dA = -mean over the two end pores of dphi_s, times the CURRENT area
                     -> A = A_0 exp(-(phi_a + phi_b)/2)
    throat radius    square: r = 0.5 sqrt(A);  circle: r = sqrt(A / pi)

The throat law is multiplicative, not a subtraction: their
``throatPrecipitatedAreaInTimeStep = sourceterm * A * dt / rho_solid`` multiplies the CURRENT
area (problem.hh::calculateThroatPrecipitatedAreaAverageSource), so the throat decays
exponentially and never closes on its own.  Verified against their own output at both
temperatures of the 1-D saturated case: phi_s = 0.3476 -> r 110 -> 100.845 um (theirs 100.850)
and phi_s = 0.8882 -> 88.096 um (theirs 88.095).  A linear subtraction agrees to first order and
was what this module did until 2026-09-17; it drove the throat to zero at phi_s > 0.7 and made
the network singular.

The solid is NaCl at ``rho_cr`` (2170 kg/m^3 in this project); DuMuX uses its molar solid
density of the same component.  A kinetic branch is NOT ported: their kinetic law needs a
reactive surface area per pore, which this project has no value for.
"""
from __future__ import annotations

import numpy as np

from ..materials import brine


class PoreClogged(RuntimeError):
    """A pore filled with precipitate.  The network has no closed-pore state: its volume, its
    throat radii and the step all go to zero, so the run has to end here rather than crawl."""

    def __init__(self, pores, fraction):
        self.pores = np.asarray(pores)
        self.fraction = float(fraction)
        super().__init__(f"{self.pores.size} pore(s) reached a solid volume fraction of "
                         f"{self.fraction:.4f} (pores {self.pores[:8].tolist()}...): the pore "
                         f"network has no closed-pore state")


class PorePrecipitation:
    """Equilibrium precipitation state of a pore network, mutating the network in place.

    ``layer`` is a ``PoreLayer``-like object with ``V`` (liquid volume per pore) and ``lat``
    (list of ``(a, b, r_throat, L_throat)``).  Construct it once, then call ``step`` after every
    transport step; it returns the dissolved salt mass that is left and records the solid.

    CONTINUITY.  The solid displaces liquid, and that liquid has to go somewhere: the pore
    releases ``last_volume_change`` [m^3] per step, i.e. ``last_volume_rate`` [m^3/s], which the
    caller must feed into the pressure solve as a SOURCE (``q -= last_volume_rate``).  Without it the network draws the whole
    evaporation from the reservoir and imports too much salt -- in the DuMuX saturated reference
    case the shrinking top pore supplies about 10 % of the throughput.

    The network is mutated because everything downstream (``solve_pressure``, ``salt_step``,
    the drivers' dt limiters) reads ``layer.V`` and ``layer.lat``: with precipitation on, those
    are no longer constants.  ``initial_V`` and ``initial_r`` keep the untouched geometry.
    """

    def __init__(self, layer, c_sat, rho_cr, pore_volume_change=True, throat_radius_change=True,
                 max_solid_fraction=0.99, T=brine.T_REF):
        self.layer = layer
        self.c_sat = float(c_sat)
        self.rho_cr = float(rho_cr)
        self.pore_volume_change = bool(pore_volume_change)
        self.max_solid_fraction = float(max_solid_fraction)
        self.T = float(T)
        self.throat_radius_change = bool(throat_radius_change)
        self.initial_V = np.array(layer.V, float)
        self.initial_r = np.array([t[2] for t in layer.lat], float)
        self.shape = getattr(layer, "throat_shape", "square")
        self.m_solid = np.zeros(len(self.initial_V))        # precipitated NaCl per pore [kg]
        self.ends = (np.array([t[0] for t in layer.lat], int),
                     np.array([t[1] for t in layer.lat], int))
        self.last_volume_change = np.zeros(len(self.initial_V))   # liquid volume LOST last step
        self.last_volume_rate = np.zeros(len(self.initial_V))     # and that over the step [m^3/s]

    # ------------------------------------------------------------------ state
    @property
    def solid_volume_fraction(self):
        """Solid volume over the INITIAL pore volume [-] (their ``precipitateVolumeFraction``)."""
        return self.m_solid / self.rho_cr / self.initial_V

    @property
    def throat_area(self):
        """Current throat cross-sectional area [m^2]: A_0 exp(-(phi_a + phi_b)/2)."""
        phi = self.solid_volume_fraction
        a, b = self.ends
        return self._area(self.initial_r) * np.exp(-0.5 * (phi[a] + phi[b]))

    def _area(self, r):
        r = np.asarray(r, float)
        return 4.0 * r ** 2 if self.shape == "square" else np.pi * r ** 2

    def _radius(self, area):
        area = np.maximum(np.asarray(area, float), 0.0)
        return 0.5 * np.sqrt(area) if self.shape == "square" else np.sqrt(area / np.pi)

    # ------------------------------------------------------------------- step
    def step(self, m_salt, dt=None):
        """Precipitate what is above saturation; update the pore volumes and throat radii.

        ``m_salt`` [kg] dissolved NaCl per pore -> returns the dissolved mass that is left.
        The volume the solid takes is removed from the pore's LIQUID volume, and the throats
        narrow by the same volumetric fraction, as in DuMuX's average-source-term approach.
        ``dt`` [s] is the step the salt arrived over; it turns the displaced volume into the
        RATE ``last_volume_rate`` the caller feeds to the pressure solve.  Pass it: dividing the
        volume by a later, shorter step instead makes the source blow up and the step collapse.
        """
        m_salt = np.asarray(m_salt, float)
        V = np.asarray(self.layer.V, float)
        w = brine.mass_fraction_from_salt_mass(m_salt, V, T=self.T)
        m_sat = brine.salt_mass(self.c_sat, V, T=self.T)
        dm = np.where(w > self.c_sat, m_salt - m_sat, 0.0)     # precipitates now [kg]
        self.last_volume_change[:] = 0.0
        self.last_volume_rate[:] = 0.0
        if not np.any(dm > 0.0):
            return m_salt
        self.m_solid += dm
        out = m_salt - dm
        if self.pore_volume_change:
            phi = self.solid_volume_fraction
            full = np.flatnonzero(phi >= self.max_solid_fraction)
            if full.size:
                raise PoreClogged(full, phi[full].max())
            V_new = (1.0 - phi) * self.initial_V
            self.last_volume_change = V - V_new                # liquid the solid displaced [m^3]
            if dt:
                self.last_volume_rate = self.last_volume_change / float(dt)
            self.layer.V = V_new
            # the solid took liquid volume away, so a pore that just precipitated is saturated in
            # a SMALLER volume: the difference precipitates as well (one shot, mass exact)
            sel = dm > 0.0
            m_left = brine.salt_mass(self.c_sat, V_new, T=self.T)
            self.m_solid[sel] += out[sel] - m_left[sel]
            out[sel] = m_left[sel]
        if self.throat_radius_change:
            # their approach 2: the throat loses the MEAN solid fraction of its two end pores,
            # applied to the area it has LEFT -- exponential decay, closed form in throat_area.
            # Last, so that the re-saturation above is already in m_solid.
            r_new = self._radius(self.throat_area)
            self.layer.lat = [(t[0], t[1], float(r), t[3]) for t, r in zip(self.layer.lat, r_new)]
        return out

    # ----------------------------------------------------------------- report
    def summary(self):
        return dict(solid_mass=float(self.m_solid.sum()),
                    max_solid_volume_fraction=float(self.solid_volume_fraction.max()),
                    min_porosity=float((1.0 - self.solid_volume_fraction).min()),
                    min_throat_radius=float(self._radius(self.throat_area).min()),
                    initial_throat_radius=float(self.initial_r.min()))
