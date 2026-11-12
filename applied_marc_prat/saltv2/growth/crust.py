"""Per-voxel salt-structure properties and the WET -> SEALED transition.

Sealing can trap liquid in a wet voxel; prior drying is not required.

Scope
-----
This module owns everything that plan/17 section 3.3 asks a salt voxel to
carry, and the transition that changes it.  It contains no vapour
solve, no growth law and no supply solve: the driver hands it fluxes and
masks, it hands back masks and property fields.

Physical model / numerical method / parameterisation (kept separate on
purpose, per CLAUDE.md):

* physical model -- the salt structure is a porous medium whose pores are
  sub-voxel (ideas/v2_盐域参数在体素模型中的实现.md section 1).  Its
  porosity, pore size, permeability and capillary threshold are voxel
  ATTRIBUTES initialised at birth. Front precipitation lowers porosity and
  updates relative permeability, radius, capillary threshold and diffusivity
  (explicit assumptions in update_properties). When the porosity of the
  reaction band or of the bulk reaches ``eps_min`` the voxel is SALT_SEALED:
  no liquid supply, no growth, vapour conducted through the closed band
  (``sealed_conductance``, plan/17 section 4.4 and 2-D README sections 8-12).
* numerical method -- explicit in time; the driver limits dt by ``dt_seal``
  (see ``dt_to_seal``) so that no voxel steps past ``eps_min``.
* parameterisation -- one birth closure, the Wijnhorst 2025 hierarchy of
  patchy porous salt (plan/13, the same stack as ``supply.supply_need``).
  The nanoporous-crust closure of Licsandru 2022 (an inversion, no hierarchy;
  knowledge base section 8 shows it is NOT the patchy stack scaled down in
  ``d_sub``) was carried as a second branch until 2026-09-14 and removed as
  unused; it is in the checkpoint commit 95c1acb if a crust branch is needed.

Provenance tags on the numbers below follow plan/14: [M] measured,
[L] literature, [A] assumed, [B] back-calibrated, [V] varied in a scan.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import SALT_SEALED, SALT_WET, Grid
from ..materials.air import D_V
from ..materials.brine import GAMMA
from ..materials.constants import M_H2O, R_GAS
from ..materials.salt import EPS_S, EPS_S_HIERARCHICAL, RHO_CR

PATCHY_POROUS_SALT = "PATCHY_POROUS_SALT"      # the closure's name in run records

# Knudsen correction of the dry-salt diffusivity: R_GAS and M_H2O are materials.constants'.
# (Until 2026-09-21 this module had its own rounded M_W = 0.018; D_K moved by -0.04 %.)

# Hierarchy constants of the patchy closure that no run has ever varied (2026-09-14: they were
# BirthClosure fields; the values are unchanged).  ``supply.supply_need`` imports them
# so the supply gate and the voxel properties come from one stack.
C_H = 1.0              # h_1 = C_H d_sub, thickness of the first layer [L]
ALPHA = 2.0            # layer-thickness ratio h_n = ALPHA^(n-1) h_1 [L]
A_MIN = 2e-9           # smallest pore the halving may reach [m] [A]
R_P_FACTOR = 0.5       # r_p(n) = R_P_FACTOR a_1(n): the Knudsen radius is half the pore
                       # half-width (implementation note section 3.1) [A]; D_eff is
                       # proportional to r_p in the Knudsen-limited regime
COS_THETA = 1.0        # perfectly wetting [A]
T_K = 294.15           # temperature [K] [M, 21 C]: the DEFAULT of BirthClosure.T


def knudsen_diffusivity(r_p, T: float) -> np.ndarray:
    """D_K = (r_p/3) sqrt(8 R T / (pi M_w))  [m^2/s]  (plan/17 section 4.1)."""
    return np.asarray(r_p, float) / 3.0 * np.sqrt(8.0 * R_GAS * T / (np.pi * M_H2O))


def effective_vapour_diffusivity(eps_s, r_p, D_v: float, T: float,
                                 exponent: float = 1.5) -> np.ndarray:
    """D_eff = eps_s^m (1/D_v + 1/D_K)^-1  [m^2/s]  (plan/17 section 4.1).

    ``exponent`` m defaults to 1.5, which is eps/tau with the Bruggeman tortuosity
    tau = eps^-0.5 [L]; the Bosanquet term makes the nanoporous branch Knudsen-limited.

    m is exposed because the sealed crust is evaluated at eps = eps_min ~ 0.02 [A], close
    enough to a percolation threshold that Bruggeman is an extrapolation rather than a
    measurement: near a threshold the conductance follows (eps - eps_c)^t with t ~ 2 and a
    non-zero eps_c, which is far steeper.  The plateau the crust predicts therefore has to be
    reported as a band over m, not as a single number (2026-09-07).  Changing it changes the
    physical closure and is never done implicitly.
    """
    D_K = knudsen_diffusivity(r_p, T)
    return np.asarray(eps_s, float) ** float(exponent) / (1.0 / D_v + 1.0 / D_K)


@dataclass(frozen=True)
class SaltStructure:
    """The microstructure of the porous salt: ONE record, because its two numbers are two readings
    of the same images and may not be chosen separately.

        crystal size   d_c(n) = 0.14 d_sub * pore_ratio^(n-1)   (layer n, counted from the substrate)
        pore           a_1(n) = d_c(n) / 4
        slits          2 pore_ratio^-(n-1) of length d_c(1): k_s(n) = c_k 2 a_1(n)^3 / (12 d_c(n))
        porosity       eps_s, the same on every layer

    ``pore_ratio`` = 1 is a single-scale structure: crystal, pore, permeability and capillary
    threshold do not depend on height.  A case names the structure it runs with
    (``run(salt_structure=)``); the drivers read both numbers from it and from nowhere else.
    """

    name: str
    pore_ratio: float      # crystal- and pore-size ratio of successive layers [-], in (0, 1]
    eps_s: float           # porosity the structure replicates itself with [-]

    def __post_init__(self):
        if not self.name:
            raise ValueError("a salt structure is recorded by name and needs one")
        if not (0.0 < self.pore_ratio <= 1.0):
            raise ValueError("pore_ratio must be in (0, 1]")
        if not (0.0 < self.eps_s < 1.0):
            raise ValueError("eps_s must be in (0, 1)")


# Single-scale packing: Prat (2026) baseline, eps_s = 0.43 [M, Dai], no hierarchy.
SINGLE_SCALE = SaltStructure("single_scale", 1.0, EPS_S)
# Hierarchical cauliflower: Wijnhorst, Prat & Shahidzadeh (2025), sizes halve per layer [L].
HIERARCHICAL = SaltStructure("hierarchical", 0.5, EPS_S_HIERARCHICAL)
# NOT self-consistent, kept ONLY to reproduce the run records made before 2026-09-21: the pore
# sizes of the hierarchical structure with the porosity of the single-scale one.
LEGACY_MIXED = SaltStructure("legacy_mixed", 0.5, EPS_S)


@dataclass(frozen=True)
class BirthClosure:
    """How a newly created salt voxel gets its properties.

    Patchy porous salt (implementation note section 3.1, = plan/13 stack):
        d_c(n) = d_c0 * r^(n-1),  d_c0 = 0.14 d_sub          [L, Wijnhorst 2025]
        a_1(n) = d_c(n)/4                                     [L]
        h_1    = C_H d_sub,  h_n = ALPHA^(n-1) h_1            [L] C_H = 1
        k_s(n) = c_k 2 a_1(n)^3 / (12 d_c(n))                 [L] x c_k [B]
                 (= c_k 2^n a_1(n)^3 / (12 d_c0) at r = 0.5)
        P_cth(n) = 2 gamma COS_THETA / a_1(n)                 [L]
        r_p(n) = R_P_FACTOR a_1(n)                            [A]
        r and eps_s: ``structure``, one record (SaltStructure)
    The other hierarchy constants are module constants (C_H, ALPHA, A_MIN, R_P_FACTOR,
    COS_THETA); the fields below are the ones runs and tests vary.
    """

    d_sub: float = 68e-6              # substrate grain diameter [m], the DEFAULT for voxels
                                      # whose own d_sub the caller does not supply.  Every
                                      # method below takes an optional per-voxel d_sub: with a
                                      # pore-size distribution, deriving all voxel properties
                                      # from one global mean is a listed anti-pattern
                                      # (plan/17 P9)
    structure: SaltStructure = LEGACY_MIXED   # pore-size ratio and porosity, inseparable.  The
                                      # default is the pairing of every record to 2026-09-20
    c_k: float = 3e-3                 # effective-conductance reduction [B]
    # --- the kinetic share of the growth law (task B, 2026-09-06)
    kinetic_efficiency: float | None = None   # F(Da*) of Prat (2026) eq (29), or None to keep
                                      # the paper's self-replication assumption.  When it is
                                      # given, a voxel is born at
                                      #     eps_birth = 1 - (1 - eps_s0) / F
                                      # because the (1 - F) share of the salt that evaporation
                                      # delivers cannot advance the interface and precipitates
                                      # inside the material as it is laid down (the reaction
                                      # zone delta_C sweeps through every element in
                                      # delta_C / V_n).  The ADVANCE is unchanged -- the same
                                      # V_n now carries denser material, and rho_cr (1-eps_birth)
                                      # V_n is exactly the salt evaporation handed over.  The
                                      # porosity reduction (1-eps_s0)(1-F)/F contains no grid
                                      # length, which is the point of task B.  ASSUMPTION [A]:
                                      # Prat's own open question 2 is where this salt goes
    reaction_band: float | None = None        # delta_C [m], Prat (2026) after eq (19),
                                      # sqrt(eps_s D_s,eff / ((1-C_eq) a_v k_r)), height
                                      # independent and ~0.6-3.4 um for 34-1080 um beads [L].
                                      # A face that CANNOT advance fills this band instead of
                                      # the whole voxel, so its closing time
                                      #     tau = (eps_s - eps_min) rho_cr delta_C / (beta j)
                                      # is a physical time and not a grid time
    # --- shared
    T: float = T_K                    # temperature of the Knudsen diffusivity [K]; a case that
                                      # declares its temperature hands it over (driver3d.run(T=))
    gamma: float = GAMMA              # brine surface tension [N/m] [L]
    D_v: float = D_V                  # free-air vapour diffusivity [m^2/s] [L]
    tortuosity_exponent: float = 1.5  # m in D_eff = eps^m (1/D_v + 1/D_K)^-1 [L, Bruggeman]

    @property
    def eps_s0(self) -> float:
        """Porosity of the self-replicating structure [-]: the structure's, never set here."""
        return self.structure.eps_s

    # -------------------------------------------------- patchy hierarchy
    def _d_sub(self, d_sub):
        return self.d_sub if d_sub is None else np.asarray(d_sub, float)

    def generation(self, z, d_sub=None) -> np.ndarray:
        """Layer index n >= 1 of a voxel born at height z [m] over a substrate of
        grain size ``d_sub`` (default: the closure's own).

        h_n = ALPHA^(n-1) h_1 stacked from the substrate gives
        n(z) = floor(log_ALPHA(1 + z (ALPHA-1) / h_1)) + 1, i.e. n = 1 for
        z < h_1.  Same convention as ``supply.supply_need``.
        """
        z = np.maximum(np.asarray(z, float), 0.0)
        h1 = C_H * self._d_sub(d_sub)
        n = np.floor(np.log(1.0 + z * (ALPHA - 1.0) / h1) / np.log(ALPHA)) + 1.0
        return np.clip(n, 1, 40).astype(np.int16)

    def pore_radius(self, n, d_sub=None) -> np.ndarray:
        """Internal pore half-width a_1(n) [m]."""
        n = np.asarray(n, float)
        a1 = 0.14 * self._d_sub(d_sub) / 4.0
        return np.maximum(a1 * self.structure.pore_ratio ** (n - 1.0), A_MIN)

    def permeability(self, n, d_sub=None) -> np.ndarray:
        """k_s(n) [m^2].  Set by the COARSEST pores along the flow path; kept
        independent of P_cth (CLAUDE.md parameter rule)."""
        n = np.asarray(n, float)
        d_c0 = 0.14 * self._d_sub(d_sub)
        a_n = self.pore_radius(n, d_sub)
        # W25 (before their eq 3): parallel slits of length d_c0 whose number doubles per layer,
        # 2^n, BECAUSE the crystals halve.  Written as 2 r^-(n-1) the count follows the crystal
        # size: 2^n exactly at r = 0.5, the first layer's 2 on every layer at r = 1 [A for other r]
        slits = 2.0 * self.structure.pore_ratio ** (1.0 - n)
        return self.c_k * slits * a_n ** 3 / (12.0 * d_c0)

    def capillary_threshold(self, n, d_sub=None) -> np.ndarray:
        """P_cth(n) [Pa].  Set by the FINEST outer pores, independently of k_s."""
        a_n = self.pore_radius(n, d_sub)
        return 2.0 * self.gamma * COS_THETA / a_n

    def birth_porosity(self, n, F=None) -> np.ndarray:
        """Porosity a voxel is born with.

        With ``kinetic_efficiency`` = F set, this is 1 - (1 - eps_s0)/F rather than eps_s0:
        the material passes through the reaction zone as it is laid down and keeps the
        (1 - F) share there (see the field's comment).  F <= 1 - eps_s0 would mean the salt
        cannot form a porous structure at all -- the crust limit -- and is rejected here rather
        than clipped, because a clipped porosity would silently change the mass balance.
        """
        n = np.asarray(n, float)
        eps = self.eps_s0
        F = self.kinetic_efficiency if F is None else np.asarray(F, float)
        if F is None:
            return np.full(n.shape, eps)
        F = np.broadcast_to(np.asarray(F, float), n.shape)
        if np.any(~np.isfinite(F)) or np.any(F <= 0.0) or np.any(F > 1.0):
            raise ValueError("kinetic_efficiency must be in (0, 1]")
        if np.any(F <= 1.0 - eps):
            raise ValueError(
                f"F <= 1 - eps_s = {1 - eps:.4f} for at least one voxel (min F = {F.min():.4f}): "
                "the (1 - F) share alone fills the structure as it forms, so no porous salt can "
                "be laid down.  This is the crust limit of the closure and needs its own branch, "
                "not a clipped porosity")
        return 1.0 - (1.0 - eps) / F

    def describe(self) -> dict:
        """Everything that was actually in force, for the run record."""
        return dict(structure_type=PATCHY_POROUS_SALT, d_sub=self.d_sub, gamma=self.gamma,
                    cos_theta=COS_THETA, D_v=self.D_v, T=self.T, eps_s0=self.eps_s0, c_k=self.c_k,
                    salt_structure=self.structure.name,
                    c_h=C_H, alpha=ALPHA, r_halving=self.structure.pore_ratio, a_min=A_MIN,
                    r_p_factor=R_P_FACTOR, kinetic_efficiency=self.kinetic_efficiency,
                    reaction_band=self.reaction_band,
                    tortuosity_exponent=self.tortuosity_exponent)


class CrustField:
    """The per-voxel property arrays of plan/17 section 3.3, plus the state
    transitions that act on them.

    The state array is the geometry engine's own array (shared, not copied),
    so a flip performed by the engine is immediately visible here and a
    transition performed here is immediately visible to the engine.

    Properties are written ONLY by ``birth`` (from an explicit flip mask) and
    ``deposit`` (front densification).  A voxel is never identified as "new"
    by testing ``eps_s == 0`` (plan/17 section 3.3, last sentence).
    """

    def __init__(self, grid: Grid, engine, closure: BirthClosure,
                 rho_cr: float = RHO_CR, d_sub_field: np.ndarray | None = None,
                 F_field: np.ndarray | None = None, band_field: np.ndarray | None = None):
        """``d_sub_field`` gives each lateral column the grain size of the
        substrate pore beneath it (shape ``grid.shape[:-1]``, or the full grid).
        Without it every voxel is born from ``closure.d_sub``, which is only
        right when the pore-size distribution is narrow -- plan/17 P9 lists
        "one global mean d_sub for every voxel property" as an anti-pattern."""
        self.grid = grid
        self.engine = engine
        self.closure = closure
        if d_sub_field is None:
            self.d_sub_field = None
        else:
            d = np.asarray(d_sub_field, float)
            if d.shape == grid.shape[:-1]:
                d = np.broadcast_to(d[..., None], grid.shape)
            if d.shape != grid.shape:
                raise ValueError("d_sub_field must have the lateral or the full grid shape")
            if not np.all(d > 0.0):
                raise ValueError("d_sub_field must be strictly positive")
            self.d_sub_field = d
        self.rho_cr = float(rho_cr)     # crystal density [kg/m^3]
        shape = grid.shape
        self.eps_s = np.zeros(shape)
        # Sub-grid reaction band (task B).  m_skin is the salt a STALLED face has put into the
        # outermost ``reaction_band`` metres of its own voxel, per unit face area [kg/m^2].  The
        # band closes at (eps_s - eps_min) rho_cr delta_C, a physical mass per area, so the
        # closing TIME contains no grid length.  The same salt is also taken off the voxel's bulk
        # porosity (spread over the whole voxel) so that salt_solid_mass stays exact; that bulk
        # share is delta_C / a of the band and never reaches eps_min on its own.
        self.m_skin = np.zeros(shape)
        # Optional per-voxel kinetic efficiency F(Da*); None falls back to closure.kinetic_efficiency.
        # It has to be settable at CONSTRUCTION, because __init__ births the voxels that already
        # exist (the seeds): assigning it afterwards would leave those with the box-mean F while
        # every later voxel got its own, and the salt_structure balance sees the difference
        # (3.1e-7 relative, found 2026-09-07 porting the pore-size spread to 2-D).
        self.F_field = self._as_voxel_field(F_field, "F_field")
        # Optional per-voxel reaction band delta_C [m]; None falls back to closure.reaction_band.
        # Needed as soon as the pores differ in size: delta_C ~ sqrt(d_sub), so a pore-size
        # distribution makes the band -- and with it the sealing time and the sealed conductance
        # -- a per-column quantity, exactly as F is (plan/17 P9: no voxel property from a box mean).
        self.band_field = self._as_voxel_field(band_field, "band_field")
        self.eps_birth = np.zeros(shape)
        self.k_birth = np.zeros(shape)
        self.r_birth = np.zeros(shape)
        self.P_birth = np.zeros(shape)
        self.generation = np.zeros(shape, np.int16)
        self.r_p = np.zeros(shape)
        self.k_s = np.zeros(shape)
        self.P_cth = np.zeros(shape)
        self.D_eff = np.zeros(shape)
        self.m_precipitated = np.zeros(shape)      # front-filled salt mass [kg]
        self.t_birth = np.full(shape, np.nan)
        self.t_sealed = np.full(shape, np.nan)
        # give the voxels that already exist (seeds) their birth properties
        existing = engine.solid_mask_full()
        if existing.any():
            self.birth(existing, 0.0)

    def _as_voxel_field(self, value, name: str):
        """Broadcast a lateral or full-grid field to the full grid; ``None`` stays ``None``."""
        if value is None:
            return None
        v = np.asarray(value, float)
        if v.shape == self.grid.shape[:-1]:
            v = np.broadcast_to(v[..., None], self.grid.shape)
        if v.shape != self.grid.shape:
            raise ValueError(f"{name} must have the lateral or the full grid shape")
        if not np.all(np.isfinite(v)) or np.any(v <= 0.0):
            raise ValueError(f"{name} must be finite and strictly positive")
        return np.array(v, float)

    # ------------------------------------------------------------------
    @property
    def state(self) -> np.ndarray:
        return self.engine.state

    def z_centres(self) -> np.ndarray:
        g = self.grid
        shape = [1] * g.ndim
        shape[-1] = g.n_z
        return np.broadcast_to(g.z_centres().reshape(shape), g.shape)

    def birth(self, flipped: np.ndarray, t: float):
        """Initialise the properties of the voxels in ``flipped`` from their
        birth height.  ``flipped`` must come from ``StepReport.flipped`` (or
        an explicit seed mask), never from a porosity test."""
        flipped = np.asarray(flipped, bool)
        if not flipped.any():
            return
        if np.any(flipped & np.isfinite(self.t_birth)):
            raise ValueError("birth() called twice on the same voxel")
        c = self.closure
        z = self.z_centres()[flipped]
        d_sub = None if self.d_sub_field is None else self.d_sub_field[flipped]
        n = c.generation(z, d_sub)
        self.generation[flipped] = n
        # F is a per-voxel field when the caller supplies one (plan/17 P9: a voxel's properties
        # must come from the pore that feeds it, not from the box mean).  Using the mean instead
        # broke the precipitation balance by the spread of F across the pores -- 0.11 % at
        # lambda = 1, where d_sub runs 1008-1133 um (found 2026-09-07).
        F_here = None if self.F_field is None else self.F_field[flipped]
        self.eps_s[flipped] = c.birth_porosity(n, F_here)
        a_n = c.pore_radius(n, d_sub)
        self.r_p[flipped] = a_n * R_P_FACTOR
        self.k_s[flipped] = c.permeability(n, d_sub)
        self.P_cth[flipped] = c.capillary_threshold(n, d_sub)
        self.eps_birth[flipped] = self.eps_s[flipped]
        self.k_birth[flipped] = self.k_s[flipped]
        self.r_birth[flipped] = self.r_p[flipped]
        self.P_birth[flipped] = self.P_cth[flipped]
        self.D_eff[flipped] = effective_vapour_diffusivity(
            self.eps_s[flipped], self.r_p[flipped], c.D_v, c.T, c.tortuosity_exponent)
        self.t_birth[flipped] = t

    # ------------------------------------------------------------------
    def deposit(self, dm: np.ndarray):
        """Add the salt mass ``dm`` [kg] precipitated inside each voxel during
        this step (plan/17 section 4.3): ``m_precipitated`` rises, ``eps_s``
        falls by ``dm / (rho_cr V_cell)``, and ``D_eff`` follows ``eps_s``.

        Only salt voxels may receive; the porosity is clipped at 0 and the
        caller is expected to have limited dt with ``dt_to_seal``.
        """
        dm = np.asarray(dm, float)
        if dm.shape != self.grid.shape:
            raise ValueError("dm must have the shape of the grid")
        if np.any((dm != 0.0) & ~self.engine.solid_mask_full()):
            raise ValueError("salt may only be deposited into salt voxels")
        if np.any(dm < 0.0):
            raise ValueError("front deposition must be non-negative (no dissolution in V1)")
        d_eps = dm / (self.rho_cr * self.grid.cell_volume)
        new_eps = self.eps_s - d_eps
        if np.any(new_eps < -1e-14):
            worst = float(new_eps.min())
            raise ValueError(
                f"front deposition overfills a voxel (eps_s would reach {worst:.3e}); "
                "limit dt with dt_to_seal (plan/17 P1a item 4)")
        self.m_precipitated += dm
        self.eps_s = np.maximum(new_eps, 0.0)
        changed = dm > 0.0
        if changed.any():
            self.update_properties(changed)

    def update_properties(self, mask=None):
        """Update transport after pore filling, anchored to independent birth k/Pc.

        [A] Relative Kozeny-Carman conductance, and a fixed-pore-count radius
        r/r0=sqrt(eps/eps0). These are explicit *evolution* assumptions, not
        a derivation of the measured absolute permeability from one pore size.
        No positive pore-radius floor is applied to a completely filled pore.
        """
        sel = self.engine.solid_mask_full() if mask is None else np.asarray(mask, bool)
        if not sel.any():
            return
        e, e0 = self.eps_s[sel], self.eps_birth[sel]
        ratio = np.maximum(e / e0, 0.0)
        self.r_p[sel] = self.r_birth[sel] * np.sqrt(ratio)
        self.k_s[sel] = self.k_birth[sel] * ratio ** 3 * ((1-e0)/(1-e)) ** 2
        self.P_cth[sel] = self.P_birth[sel] / np.sqrt(np.maximum(ratio, 1e-30))
        rp = self.r_p[sel]
        deff = np.zeros(e.shape)
        positive = rp > 0.0
        deff[positive] = effective_vapour_diffusivity(
            e[positive], rp[positive], self.closure.D_v, self.closure.T,
            self.closure.tortuosity_exponent)
        self.D_eff[sel] = deff

    # Numerical tolerance on the sealing threshold [-].  ``dt_to_seal`` makes the
    # limiting voxel land ON eps_min, but deposit() recomputes the porosity drop as
    # dm / (rho_cr V) with dm = rate * dt, which differs from the room it was cut
    # from by an ulp; without a tolerance ``seal`` misses the voxel and the next
    # dt_to_seal is that ulp, i.e. a zero-length step repeated forever (found by
    # the plan/17 verification cases, 2026-09-05).  1e-12 is far below any
    # physical porosity difference and far above double-precision round-off.
    EPS_TOL = 1e-12

    def dt_to_seal(self, dm_rate: np.ndarray, eps_min: float) -> float:
        """Largest dt that no voxel uses to step past ``eps_min``
        (plan/17 P1a item 4, the ``dt_seal`` limiter).  ``inf`` if nothing is
        densifying.  The analytic sealing time of section 4.4 is the same
        expression evaluated once with a constant rate.

        Voxels already within ``EPS_TOL`` of ``eps_min`` are not counted: they
        belong to ``seal`` (which the driver calls after every deposit), and
        counting them would return a zero-length step."""
        dm_rate = np.asarray(dm_rate, float)
        active = (dm_rate > 0.0) & (self.eps_s > eps_min + self.EPS_TOL)
        if not active.any():
            return np.inf
        room = (self.eps_s[active] - eps_min) * self.rho_cr * self.grid.cell_volume
        return float(np.min(room / dm_rate[active]))

    # -------------------------------------------------------------- reaction band
    def reaction_band(self) -> np.ndarray:
        """Per-voxel band thickness delta_C [m], on the full grid.

        ``band_field`` when the caller set one (pores of different size have different bands),
        otherwise ``closure.reaction_band`` broadcast.  Raises when neither exists rather than
        inventing a thickness.
        """
        if self.band_field is not None:
            return self.band_field
        if self.closure.reaction_band is None:
            raise ValueError("needs a reaction band (delta_C): set closure.reaction_band or "
                             "CrustField.band_field, or use the bulk deposit path explicitly")
        return np.broadcast_to(float(self.closure.reaction_band), self.grid.shape)

    def skin_capacity(self, eps_min: float) -> np.ndarray:
        """Salt per unit face area that closes the reaction band [kg/m^2].

        ``(eps_birth - eps_min) rho_cr delta_C``: the band is delta_C thick and closes when its
        own porosity has fallen to ``eps_min``.  Uses the voxel's BIRTH porosity, so a voxel
        whose bulk has already been densified does not get a smaller band.
        """
        band = self.reaction_band()
        cap = np.maximum(self.eps_birth - eps_min, 0.0) * self.rho_cr * band
        # Validity of the sub-grid band (2026-09-07).  Closing the band costs the voxel's BULK
        # pore space cap / (rho_cr a) = (eps_birth - eps_min) delta_C / a, which must not exceed
        # the pore space it has:
        #     delta_C / a  <=  eps_birth / (eps_birth - eps_min)
        # Above that the band is so much thicker than the voxel that closing it would need more
        # salt than the voxel can hold, and the porosity goes negative inside deposit_skin.  The
        # condition is refused here, where it can be named, instead of surfacing as an overfill.
        # NOTE this is the OPPOSITE end of the mesh from the conductance requirement of README
        # section 9.2 (a <= delta_C, so that a sealed voxel is closed over its whole height rather
        # than over delta_C / a of it).  The two meet only at a ~ delta_C, which is therefore the
        # one mesh at which the sealing model is both mass-consistent and grid-converged.
        born = self.eps_birth > eps_min
        if born.any():
            cost = cap[born] / (self.rho_cr * self.grid.a)
            if np.any(cost > self.eps_birth[born]):
                i = int(np.argmax(cost - self.eps_birth[born]))
                ratio = float(band[born][i] / self.grid.a)
                raise ValueError(
                    f"reaction band too thick for the mesh: delta_C / a = {ratio:.3f} exceeds "
                    f"eps_birth / (eps_birth - eps_min) = "
                    f"{self.eps_birth[born][i] / (self.eps_birth[born][i] - eps_min):.3f}, so "
                    "closing the band would need more salt than the voxel's pore space holds.  "
                    "Refine to a >= delta_C * (1 - eps_min/eps_birth), or coarsen; a ~ delta_C is "
                    "the mesh at which the band model and the sealed conductance are both valid")
        return cap

    def deposit_skin(self, ds: np.ndarray, eps_min: float) -> None:
        """Add ``ds`` [kg/m^2 of face] to the reaction band of each voxel.

        The same mass is removed from the voxel's bulk pore space, so the structure's total
        solid mass is right whichever way it is read; ``dt_to_close_skin`` keeps the band from
        stepping past its capacity, and ``seal_skin`` turns a closed band into SALT_SEALED.
        """
        ds = np.asarray(ds, float)
        if ds.shape != self.grid.shape:
            raise ValueError("ds must have the shape of the grid")
        if np.any(ds < 0.0):
            raise ValueError("band deposition must be non-negative (no dissolution in V1)")
        if np.any((ds != 0.0) & ~self.engine.solid_mask_full()):
            raise ValueError("salt may only be deposited into salt voxels")
        self.m_skin = self.m_skin + ds
        # bulk bookkeeping: the same salt in kg, so salt_solid_mass and the ledger agree
        self.deposit(ds * self.grid.face_area)

    def dt_to_close_skin(self, ds_rate: np.ndarray, eps_min: float) -> float:
        """Largest dt that no reaction band uses to step past its capacity ([kg/m^2/s] in)."""
        ds_rate = np.asarray(ds_rate, float)
        cap = self.skin_capacity(eps_min)
        active = (ds_rate > 0.0) & (self.m_skin < cap - self.EPS_TOL * self.rho_cr)
        if not active.any():
            return np.inf
        return float(np.min((cap[active] - self.m_skin[active]) / ds_rate[active]))

    def seal_skin(self, eps_min: float, t: float) -> np.ndarray:
        """Any salt voxel whose reaction band has closed becomes SALT_SEALED."""
        cap = self.skin_capacity(eps_min)
        salt = self.engine.solid_mask_full()
        changed = salt & (self.state != SALT_SEALED) & (self.m_skin >= cap - self.EPS_TOL * self.rho_cr)
        if changed.any():
            self.state[changed] = SALT_SEALED
            self.engine.f[changed] = 0.0
            self.t_sealed[changed] = t
        return changed

    # ------------------------------------------------------------------ states
    def seal(self, eps_min: float, t: float) -> np.ndarray:
        """Any salt voxel whose porosity has reached ``eps_min`` becomes
        SALT_SEALED (plan/17 section 4.4).

        eps_min > 0 is an assumed connectivity threshold, NOT zero pore
        volume.
        """
        salt = self.engine.solid_mask_full()
        changed = salt & (self.state != SALT_SEALED) & (self.eps_s <= eps_min + self.EPS_TOL)
        if changed.any():
            self.state[changed] = SALT_SEALED
            self.engine.f[changed] = 0.0
            self.t_sealed[changed] = t
        return changed

    # ------------------------------------------------------------------ vapour
    def sealed_conductance(self, eps_min: float, D_v: float | None = None) -> np.ndarray:
        """``D_eff / D_v`` [-] of a voxel sealed by ``seal_skin`` (added 2026-09-07).

        Physical model.  ``seal_skin`` closes the REACTION BAND -- a sub-voxel layer of
        thickness ``delta_C`` at the evaporating face (Prat 2026, after eq 19) -- not the
        whole voxel: the band reaches ``eps_min`` while the bulk pore space behind it is
        still at ``eps_s``.  The voxel's vapour conductance is therefore the two in series
        over the voxel height ``a``::

            a / D_cell = delta_C / D(eps_min, r_p(eps_min)) + (a - delta_C) / D(eps_s, r_p)

        with ``D`` the Bruggeman/Bosanquet closure ``effective_vapour_diffusivity`` and the
        pore radius carried down by the same fixed-pore-count law ``update_properties`` uses,
        ``r_p / r_birth = sqrt(eps / eps_birth)``.  A voxel sealed by ``seal`` instead (bulk
        porosity itself at ``eps_min``) has ``eps_s = eps_min`` and the two terms coincide,
        so one expression covers both sealing routes.

        Numerics.  ``delta_C >= a`` means the band fills the voxel; the bulk term is then
        absent (fraction clamped at 1).  A voxel whose bulk pore radius has closed completely
        gets conductance 0 -- a true wall -- and ``porous_argument`` leaves it out of the mask
        rather than handing ``solve_vapor`` a zero factor.

        Returned on the whole grid; only the SALT_SEALED entries are meaningful.
        """
        D_v = self.closure.D_v if D_v is None else float(D_v)
        T = self.closure.T
        f_band = np.minimum(self.reaction_band() / self.grid.a, 1.0)   # [-], per voxel
        eps_min = float(eps_min)
        # band: at eps_min, with the pore radius the shrink law gives there
        ratio = np.divide(eps_min, self.eps_birth, out=np.zeros(self.grid.shape),
                          where=self.eps_birth > 0.0)
        r_band = self.r_birth * np.sqrt(np.maximum(ratio, 0.0))
        D_band = np.zeros(self.grid.shape)
        pos = r_band > 0.0                       # unborn voxels have r = 0: D_K = 0, no division
        D_band[pos] = effective_vapour_diffusivity(eps_min, r_band[pos], D_v, T,
                                                   self.closure.tortuosity_exponent)
        # A voxel whose band already fills it needs no bulk term; one that does not, but whose
        # bulk pore radius has closed completely, is a true wall and stays out of the result.
        bulk_term = np.zeros(self.grid.shape)
        has_bulk = (f_band < 1.0) & (self.D_eff > 0.0)
        bulk_term[has_bulk] = (1.0 - f_band[has_bulk]) / self.D_eff[has_bulk]
        ok = (D_band > 0.0) & ((f_band >= 1.0) | (self.D_eff > 0.0))
        out = np.zeros(self.grid.shape)
        out[ok] = 1.0 / ((f_band[ok] / D_band[ok] + bulk_term[ok]) * D_v)
        return out

    def porous_argument(self, eps_min: float, D_v: float | None = None):
        """``(mask, factor)`` for ``vapor.solve_vapor(porous=...)``: the SEALED voxels with
        their ``sealed_conductance(eps_min)`` as ``D_eff / D_v``, so the crust conducts vapour
        through its own closed reaction band and the leak is resolved instead of
        parameterised.  Voxels whose conductance is zero stay walls (left out of the mask);
        ``None`` when nothing is sealed, which is what ``solve_vapor`` wants in that case.
        ``eps_min`` is required: there is no wall-only mode any more (removed 2026-09-14; it
        was the pre-2026-09-07 physics and needed the ``phi_vent`` leak on top).
        """
        D_v = self.closure.D_v if D_v is None else float(D_v)
        cond = self.sealed_conductance(eps_min, D_v)
        mask = (self.state == SALT_SEALED) & (cond > 0.0)
        if not mask.any():
            return None
        factor = np.ones(self.grid.shape)
        factor[mask] = cond[mask]
        return mask, factor

    # ------------------------------------------------------------------ report
    def volumes(self) -> dict:
        """Solid volume [m^ndim] in each state, for the run record."""
        V = self.grid.cell_volume
        return {"wet": float((self.state == SALT_WET).sum()) * V,
                "sealed": float((self.state == SALT_SEALED).sum()) * V}

    def salt_solid_mass(self) -> float:
        """Salt mass held by the structure [kg] = sum over salt voxels of
        (1 - eps_s) V_cell rho_cr.  Because ``deposit`` lowers eps_s by exactly
        dm / (rho_cr V_cell), this equals the birth solid of every voxel plus
        everything precipitated into their pores since."""
        salt = self.engine.solid_mask_full()
        return float(np.sum(1.0 - self.eps_s[salt])) * self.grid.cell_volume * self.rho_cr
