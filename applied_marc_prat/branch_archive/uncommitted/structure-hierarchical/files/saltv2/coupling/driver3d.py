"""The 3-D closed loop on the substrate pore network.  Regime: evaporation-wicking.

evaporation (quasi-steady vapour field, growth.vapor) -> capillary supply through the pore
network (pnm.PoreLayer, solve_pressure) -> salt transport between pores (pnm.salt_step) ->
growth of the external porous salt (Prat 2026 law, growth.geometry) with the hierarchical
supply gate (growth.supply.first_failure_height) -> densification of the reaction band and
sealing (growth.crust) -> the sealed crust conducts vapour -> evaporation.

Entry point: ``run``.  The production cases are ``saltv2.cases.scan_crust``; the named history
schema is ``HIST_COLUMNS`` / ``col``.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from ..growth import morphology as mo
from ..growth.crust import LEGACY_MIXED, T_K, BirthClosure, CrustField
from ..growth.geometry import GeometryEngine
from ..growth.grid import GAS, SALT_SEALED, SALT_WET, Grid
from ..growth.prat import prat
from ..growth.supply import first_failure_height
from ..growth.vapor import VaporProperties, solve_vapor_adaptive
from ..materials import brine
from ..materials.air import LEGACY as LEGACY_AIR, moist_air
from ..materials.brine import C_SAT, D_NACL, GAMMA, RHO_L, a_w
from ..materials.salt import RHO_CR
from ..pnm import PoreLayer, network_step, solve_pressure
from .interface import PoreOpenings
from .ledger import (DEFAULT_CHANNELS, INTERFACE_CHANNELS, LIQUID_CHANNELS, PORE_SALT_CHANNELS,
                     MassLedger)


# Values no run has ever varied (2026-09-14: they were run() keywords with these defaults).
S_NUC = 1.0                 # nucleation threshold C >= S_NUC * C_sat in a bare pore.  1.0 since
                            # 2026-09-19 (user's decision): NO supersaturation, the first crystal
                            # appears when the pore reaches C_sat, as in DuMuX.  It was 1.2 [L]
                            # (Rodriguez-Rodriguez 2026, 0.320 / 0.267; Veran-Tissoires & Prat 2014
                            # table 3, 1.05-1.26), and the pore was then reset from 1.2 C_sat to
                            # C_sat in one step -- a discontinuity no source describes: neither
                            # paper models the pore AFTER its first crystal, and Prat (2026) keeps
                            # the brine at C_eq by first-order kinetics.  Doing the relaxation
                            # properly needs a pore-structure exchange law; without one, no
                            # overshoot is the consistent choice.
F_MAX = 0.5                 # largest occupancy a candidate CELL may gain in one step [N]
DT_CAP = 900.0              # upper bound on dt [s] [N]
VAPOR_RTOL = 1e-8           # relative CG tolerance of the vapour solve [N].  E is a flux
                            # difference against a Dirichlet-dominated right-hand side, so it is
                            # sensitive to it: 1e-7 left ~1e-3 in E, 1e-6 ~1e-2, 1e-5 changed the
                            # morphology; 1e-8 gives ~2e-4 (README section 9.8).  Passed
                            # explicitly: solve_vapor_adaptive's own default is 1e-7.
SUPPLY_ITER_MAX = 50        # outer sweeps of the supply fixed point before it refuses [N]
RELEASE_TIMES = 2.0         # crystallised salt is handed to the faces over this many pore-exchange
                            # times (the network limiter's V / exchange of the step before) [N].  Not
                            # tied to dt itself: an event step of ~0 s would stall the release, and
                            # the growth limiter shrinks when the release grows (a Zeno feedback)
SALT_QUANTUM = 9.7e-26      # mass of one NaCl formula unit, M_NACL / N_A [kg]
MIN_BAND_CELLS = 6.0        # minimum voxels across the physical reaction band [N diagnostic]
DENSE_BIRTH_TOL = 1e-3      # porosity tolerance for "born at eps_min" under layer_av [N]: the regulariser mixes
                            # F across neighbouring columns, so a cell at the clamp is born within ~1e-3 of eps_min


# ---- named output schema (plan/17 P6 "output upgrade").  hist stays a 2-D array so that it
# compresses into the npz, but every column is addressed through HIST / col(), never by a bare
# index.  Records carry ``hist_columns`` by name, so a reader never depends on the position.
# SCHEMA_VERSION "plan17-P6-2" (2026-09-14) dropped the columns of the removed legacy paths:
# V_dry_open, h_max_wet, n_rewet, h_front (drying state), m_front (a ledger duplicate),
# E_vent, f_vent_open (the parameterised leak); "plan17-P6-1" records still read by name.
SCHEMA_VERSION = "plan17-P6-3"
HIST_COLUMNS = [
    "t_h",               # time [h]
    "theta",             # substrate coverage (all salt states)
    "n_clusters",        # clusters of salt (all states)
    "spanning",
    "h_max",             # [m], all salt states
    "h_mean",            # [m]
    "E_total",           # M_total / M_ref, the resolved gas-phase evaporation
    "j_tip_over_plane",
    "n_frozen",          # supply-limited clusters in THIS step
    "j_bare_over_plane",
    "A_r",               # exposed wet-face area / salt footprint area
    "E_base",            # M_base / M_ref   (pore openings, bare or under the crust)
    "E_salt",            # M_salt / M_ref   (wet salt faces)
    "E_base_bare",       # bare liquid pore-opening contribution to E_base / M_ref
    "E_base_crust",      # pore -> conducting sealed-crust contribution to E_base / M_ref
    "E_wet_gas",         # wet-salt -> free-gas contribution to E_salt / M_ref
    "E_wet_crust",       # wet-salt -> conducting sealed-crust contribution to E_salt / M_ref
    "V_wet",             # [m^3]
    "V_sealed",          # [m^3]
    "f_blocked",         # fraction of columns whose evaporating surfaces have no gas path out
    "dt",                # [s] step actually taken
    "limiter",           # index into LIMITERS: which limit set dt
    "n_sealed_cells",
    "res_water",         # relative residual of the network water balance
    "res_salt",          # relative residual of the structure salt balance
    "n_supply_iter",     # iterations the supply fixed point needed
]
HIST = {name: k for k, name in enumerate(HIST_COLUMNS)}
LIMITERS = ("grow", "network", "seal", "salt", "cap", "end")


def col(hist, name):
    """Column ``name`` of a hist array (plan/17 P6: no bare ``hist[:, i]``)."""
    return np.asarray(hist)[:, HIST[name]]


def run(lam=1, C0=0.25, c_k=3e-3, t_end=2.0 * 86400, z_top_mm=5.0, Lt_res=0.7e-3,
        verbose=True, n_side=None, vpp=7, delta_eff_mm=None, tag="", rh=None, snap_hours=None,
        cv=0.0, pitch_factor=1.0, seed=0, dt_factor=1.0, ell_over_dsub=0.0, eps_min=0.02,
        min_band_cells=MIN_BAND_CELLS, pore_radius=None, n_layers=1, dz_pore=None,
        pore_precipitation=False, substrate_dgf=None, substrate_L_box=None, d_sub=None,
        record_hours=None, salt_structure=None, T=None, layer_av=False):
    """3-D evaporation-wicking loop with the closed feedback of plan/17 (P5 + P6):

        gas-phase diffusion -> local flux j -> growth V_n = G F(Da*) j
          -> liquid demand -> supply fixed point (q -> p -> h_lim -> active faces)
          -> reaction band / bulk densification -> WET -> SEALED
          -> the sealed crust conducts vapour -> next vapour field

    Substrate and box
      lam           substrate scaling: pore radius 270 um / lam, pitch 0.7 mm / lam
      d_sub         grain diameter of the substrate [m], which the salt closures are written in
                    (crystal size 0.14 d_sub, Da*, reaction band, first hierarchy layer, ell).
                    A case value.  None: 4 r of each surface pore, the rule the generated
                    lattice was designed with (270 um -> 1080 um) [A]; a pore network has no
                    grains, so a case that knows its bead size states it here
      C0            reservoir and initial pore concentration [mass fraction]
      Lt_res        reservoir throat length [m] (22 mm = the 2014 substrate height).  With
                    n_layers > 1 it is the throat under the BOTTOM layer only; the column
                    itself resolves (n_layers - 1) * dz_pore of that depth
      n_layers      pore layers stacked downwards (1: the surface layer only, as before).  Only
                    the bottom layer touches the reservoir, which holds it at C0, at the
                    pressure datum and saturated (wicking; saturation is not a variable)
      dz_pore       vertical pore spacing [m] (None: the lateral pitch)
      pore_precipitation  let the salt above C_sat crystallise INSIDE the pores
                    (pnm.precipitation, ported from the DuMuX mineralisation model): the solid
                    takes liquid volume, narrows the throats and releases the displaced liquid
                    into the network.  Default False = the production path, where the pores are
                    only capped at C_sat (channel salt_pnm_cap) and all crystallisation happens
                    in the voxel domain.  With it on the pore concentration can never exceed
                    C_sat, so the S_NUC supersaturation trigger can never fire: nucleation then
                    starts where a pore FIRST precipitates.  The two nucleation rules are not
                    the same physics; runs with this option are not comparable to the SCAN_CRUST
                    records
      substrate_dgf  path to a DuMuX dgf to use AS the substrate network instead of the generated
                    lattice (PoreLayer.from_dgf): the file's pores, volumes, throat radii and
                    throat lengths verbatim, the bottom pore its reservoir.  Intended for
                    pnm_1d_homogeneous.dgf, the column the model was validated against.  It is a
                    single column, so the box carries ONE surface pore of that file's radius
                    (its grain size is the case's to state, see d_sub), and lam, Lt_res, n_side, cv,
                    pitch_factor, pore_radius, n_layers and dz_pore are then ignored
      n_side        pores per side of the periodic box (None -> 5*lam, 3.5 mm box)
      substrate_L_box  with substrate_dgf: lateral size of the periodic box the surface pore opens
                    into [m], i.e. the spacing of neighbouring surface pores, which a 1-D file
                    cannot state (None: the column's own pore spacing).  Equal to the opening,
                    2 r, the gas problem is one-dimensional and IS the DuMuX boundary layer
      cv, pitch_factor, seed   pore-size distribution and pitch widening of PoreLayer
      z_top_mm      height of the voxel box (salt may grow up to it)
      delta_eff_mm  physical far-field distance where X = X_top; the excess over z_top_mm
                    becomes a 1-D external segment (Robin top).  None -> delta_eff = z_top
      rh            far-field relative humidity; required (a case value, no kernel default)
      T             temperature of the experiment [K], a case value.  Given, ONE temperature
                    serves the gas side (materials.air.moist_air: D_v, W_sat, rho_g from the
                    DuMuX expressions), the brine density of the network, the precipitation
                    coefficient k_r (Arrhenius, prat.kr_at_temperature) and the Knudsen
                    diffusivity of the salt.  NOT carried by it [A]: the water activity (Pitzer
                    parameters at 25 C, as in DuMuX), D_NACL, C_SAT, the viscosity, the surface
                    tension and the pure-water density RHO_W (21 C values; each moves by well
                    under 1 % between 20 and 22 C, the surface tension and viscosity by ~2 %).
                    None: the constants of every record to 2026-09-20 -- gas side 20 C
                    (materials.air.LEGACY), brine and k_r 21 C, Knudsen 21 C
      vpp           voxels per pore pitch

    Closures
      salt_structure  the microstructure of the salt, a ``growth.crust.SaltStructure`` (pore-size
                    ratio of successive layers AND porosity, one record): SINGLE_SCALE or
                    HIERARCHICAL.  A case value.  None: LEGACY_MIXED, the pairing of every record
                    to 2026-09-20 (hierarchical pore sizes, single-scale porosity), which is not
                    self-consistent and is kept to reproduce those records
      c_k           effective-conductance reduction of the salt hierarchy [B, calibrated]
      eps_min       porosity at which a densifying voxel becomes SALT_SEALED [A]
      layer_av      B' (2026-09-22): F follows the crystals of each layer, Da*_n = Da*_1 r^(n-1),
                    and is clamped at F_min = (1 - eps_s)/(1 - eps_min); a voxel born at the clamp
                    is DENSE (eps_min) and sealed at once, its pore radius shrunk by the
                    fixed-pore-count law (crust.BirthClosure.layer_av).  Default off: unchanged
      ell_over_dsub regularisation length ell = ell_over_dsub * mean(d_sub) (P1b).  It gives the
                    lattice Laplacian growth a fastest wavelength ~ ell; without it the only
                    cut-off is the voxel size and the morphology does not converge in a.
                    0 keeps the unregularised law -- mechanism runs only, no quantitative claim.
                    ell must never be calibrated on the 2014 data.

    Numerics (plan/17 P1; no physics)
      dt_factor r   scales the accuracy limiters at once,
                    dt = min(r dt_grow, r dt_net, dt_seal, r dt_salt, r DT_CAP, dt_end)  (P1a item 4);
                    the event limiters dt_seal and dt_end are not scaled (they land ON the event).
                    A dt-convergence scan must vary r: in the baseline cases dt is set by
                    dt_net, so tightening the growth limit alone changes nothing (README
                    section 8, finding 1).  The `limiter` history column records which bound.
      snap_hours    keep a copy of the state every snap_hours [h] (None: no snapshots)
      record_hours  write a row of hist / pore_C / M_liquid_in only when t has advanced by this
                    much [h] since the last row (the first and the last state are always
                    written).  None: one row per step.  A physical interval, not a step count;
                    it changes what is recorded and nothing that is computed.  A column of a
                    thousand pores stepped a million times cannot keep every row

    The conducting crust, base salt and front tracking (ported from the 2-D driver, 2026-09-11):
      * SEALED voxels conduct vapour with ``CrustField.sealed_conductance``, so the crust leak
        is RESOLVED by the vapour solve;
      * wet faces that evaporate INTO a sealed voxel (no gas neighbour, so ``exposed_faces``
        does not see them) are picked up from the field, draw liquid like any face, and hand
        their salt to the sealed voxel in front until its bulk is down to eps_min (FRONT
        TRACKING: the front finishes a voxel before moving on); a finished voxel passes the
        salt on to the reaction band of the wet voxel behind it;
      * the pore opening beneath a sealed z = 0 voxel joins the bottom boundary, so vapour
        leaves the substrate THROUGH the crust (``fld.base_flux`` there), that water is a sink
        of the pore network and its salt precipitates into the z = 0 voxel (BASE SALT; it keeps
        lowering that voxel's porosity and hence its conductance);
      * the vapour field is refreshed whenever a sealed voxel's D_eff changes, not only when a
        state count changes.
    Where the structure's salt comes from (2026-09-19).  A pore under salt touches crystal, so
    its brine is held at C_sat and what the network carries into it beyond saturation
    crystallises (pnm.salt_step, implicit, with a complementarity: when salt diffuses away
    faster than it arrives the pore is released and falls below C_sat, and nothing grows).  That
    salt -- the supersaturation a bare pore built up before it nucleated included -- is ALL the
    structure gets; the growth law no longer charges the pores.  It waits in a per-pore
    interface buffer and is handed to the faces in proportion to the water they evaporate, at
    most beta(C_sat) per kg: Prat (2026) eq (29) is the UPPER BOUND of the closure, reached when
    the network delivers enough (any reservoir at or above saturation), and undercut when it
    does not (an undersaturated reservoir behind a short column).  What the network delivers
    BEYOND the bound stays in the buffer -- crystal at the pore-structure contact with no volume
    of its own [A]; this is what the saturation cap used to discard.  Balance:
    salt_crystallised = salt_supplied + salt_nucleated + salt_interface_buffer.
    Still open: a wet cluster no pore feeds keeps beta(C_sat) and books salt_unsourced (neither
    its water nor its salt has a source; it would really dry out); the 2-D driver has no pore
    network and keeps beta(C_sat).
    """
    if substrate_L_box is not None and substrate_dgf is None:
        raise ValueError("substrate_L_box is the box of a dgf substrate; without substrate_dgf "
                         "the box follows from lam, n_side and pitch_factor")
    if substrate_dgf is None:
        layer = PoreLayer(lam, Lt_res=Lt_res, n_side=n_side, cv=cv, pitch_factor=pitch_factor,
                          seed=seed, pore_radius=pore_radius, n_layers=n_layers, dz=dz_pore)
    else:
        # the DuMuX network itself; lam, Lt_res, n_side, cv, pitch_factor, pore_radius, n_layers
        # and dz_pore all describe the generated lattice and say nothing here
        layer = PoreLayer.from_dgf(substrate_dgf, L_box=substrate_L_box)
    surf = layer.surface                    # the pores that open into the gas domain
    a = (layer.L_box / layer.n_side) / vpp
    n_lat = int(round(layer.L_box / a)); nz = int(round(z_top_mm * 1e-3 / a))
    delta_eff = (z_top_mm if delta_eff_mm is None else delta_eff_mm) * 1e-3
    L_ext = delta_eff - nz * a
    if L_ext < -a:                                  # rounding of nz may overshoot by < 1 voxel
        raise ValueError("delta_eff_mm must be >= z_top_mm")
    L_ext = max(L_ext, 0.0)
    if rh is None:
        raise ValueError("rh (far-field relative humidity) is a case value and must be given")
    if T is not None and pore_precipitation:
        raise ValueError("T is not carried through pnm.precipitation; run that option without T")
    air = LEGACY_AIR if T is None else moist_air(T, rh)
    T_liquid = brine.T_REF if T is None else float(T)    # brine density of the network [K]
    kr = prat.KR_T21 if T is None else float(prat.kr_at_temperature(T))
    X_top = rh * air.W_sat                               # far-field vapour mass fraction
    g = Grid((n_lat, n_lat, nz), a)
    eng = GeometryEngine(g)
    # where each surface pore opens into the grid (coupling/interface.py owns the footprint)
    openings = PoreOpenings(g, layer.x[surf], layer.y[surf], layer.r[surf], layer.L_box)
    foot = openings.foot
    pore_id_map = openings.pore_id
    zidx = np.broadcast_to(np.arange(nz)[None, None, :], g.shape)
    zc_m = (zidx + 0.5) * a
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    lat_a = np.array([t_[0] for t_ in layer.lat], int)      # throat endpoints, for the dt limiter
    lat_b = np.array([t_[1] for t_ in layer.lat], int)
    if d_sub is not None and not d_sub > 0.0:
        raise ValueError("d_sub must be a positive grain diameter [m]")
    grain_size_given = None if d_sub is None else float(d_sub)     # goes into the run record
    # grain size seen by the gas domain, per surface pore: the case's, or 4 r [A]
    d_sub = 4 * layer.r[surf] if d_sub is None else np.full(surf.size, float(d_sub))
    structure = LEGACY_MIXED if salt_structure is None else salt_structure
    eps_s0 = structure.eps_s             # porosity of the self-replicating structure, read ONCE
    # Ds: the ONE salt diffusivity of the model, the pore network's (materials.brine.D_NACL =
    # 1.54e-9, DuMuX).  Until 2026-09-21 the Prat closures ran on prat2026's own default 1.3e-9 --
    # right for reproducing the paper, but a second diffusivity for the same brine in one step.
    F = np.array([prat.F(prat.Da_star(eps_s=eps_s0, d_sub=d, Ds=D_NACL, kr=kr)) for d in d_sub])
    BETA_SAT = C_SAT / (1 - C_SAT)       # salt per kg of water in SATURATED brine [-]; what a
                                         # face hands over is beta of the brine that feeds it
    props = VaporProperties(rho_D=air.rho_D, X_wet=air.X_salt, X_top=X_top)
    A_cell = g.face_area
    j_plane = air.rho_D * (air.X_salt - X_top) / delta_eff
    d_sub_mean = float(np.mean(d_sub))
    # Per-voxel salt-structure properties (plan/17 P3).  The birth hierarchy is the same stack
    # that supply_need uses, so a voxel's properties and the supply test that decides whether it
    # can be fed cannot disagree (tests/test_P3_state_kernel.py).
    # Task B (2026-09-06): the (1 - F) share the growth law cannot advance is carried by the
    # BIRTH POROSITY of the material as it is laid down, eps_birth = 1 - (1 - eps_s0)/F, and a
    # face that cannot advance fills the physical reaction band delta_C instead of a whole voxel.
    # Neither contains a grid length.  F and delta_C are per pore in principle; with cv = 0 every
    # pore is identical and the mean is exact.  A pore-size distribution needs them per voxel
    # (plan/17 P9 forbids deriving voxel properties from one global mean), which is not done yet.
    F_mean = float(np.mean(F))
    delta_C = float(np.mean([prat.boundary_layer_thickness(eps_s=eps_s0, d_sub=d, Ds=D_NACL, kr=kr) for d in d_sub]))
    # The wet-to-sealed transition uses the physical reaction-band thickness
    # delta_C.  Fewer than this many cells makes its sealing time and vapour
    # resistance lattice-dependent.  This is diagnostic only: it records a
    # resolution warning without changing the physical model or stopping a run.
    band_cells = delta_C / a
    band_resolved = bool(band_cells >= min_band_cells)
    closure = BirthClosure(d_sub=d_sub_mean, structure=structure, c_k=c_k, gamma=GAMMA,
                           kinetic_efficiency=F_mean, reaction_band=delta_C, D_v=air.D_v,
                           T=T_K if T is None else float(T), layer_av=bool(layer_av))
    # F below this cannot lay down porous material: (1 - F)(1 - eps_s0)/F >= eps_s0 - eps_min.
    # With layer_av the driver clamps F(z) here and the voxel is born dense; without it the
    # closure raises at 1 - eps_s0 as before (never reached by the recorded cases).
    F_min = (1.0 - eps_s0) / (1.0 - eps_min)
    # each column inherits the grain size of the pore it sits on, so that a run with a pore-size
    # distribution (cv > 0) does not give every voxel the properties of the mean grain
    # (plan/17 P9 anti-pattern list).  Columns between pores fall back to the mean.
    d_sub_map = np.where(pore_id_map >= 0, d_sub[np.maximum(pore_id_map, 0)], d_sub_mean)
    # F is per column for the same reason.  Even at cv = 0 the dgf pores are 270 +- 10 um, so F
    # spans about 0.3 %; using the box mean for the birth porosity while the growth used the
    # local F broke the precipitation balance by exactly that (found 2026-09-07).
    Da_map = np.asarray(prat.Da_star(eps_s=eps_s0, d_sub=d_sub_map, Ds=D_NACL, kr=kr), float)
    F_map = np.asarray(prat.F(Da_map), float)
    # The reaction band is per column for the same reason (plan/17 P9; delta_C ~ sqrt(d_sub)):
    # it sets each column's band capacity, sealing time and sealed conductance.  Both fields go
    # to the constructor so that any voxel born there gets its own value, as in the 2-D driver.
    band_map = np.asarray(prat.boundary_layer_thickness(eps_s=eps_s0, d_sub=d_sub_map, Ds=D_NACL, kr=kr), float)
    crust = CrustField(g, eng, closure, rho_cr=RHO_CR, d_sub_field=d_sub_map,
                       F_field=F_map, band_field=band_map)
    ell = ell_over_dsub * d_sub_mean
    precip = None
    if pore_precipitation:
        from ..pnm.precipitation import PoreClogged, PorePrecipitation
        precip = PorePrecipitation(layer, C_SAT, RHO_CR)
    led = MassLedger(totals={c: 0.0 for c in DEFAULT_CHANNELS + INTERFACE_CHANNELS + LIQUID_CHANNELS
                             + PORE_SALT_CHANNELS})
    # the pores' dissolved salt: what the reservoir hands over either leaves the brine (as
    # crystal at a covered pore, as in-pore precipitate, through the cap) or is still in the
    # pores.  The left side is metered at the reservoir throats, the right side comes from the
    # crystallisation and from the two states: this is the GLOBAL salt balance of the network.
    led.register("salt_pores", plus=["salt_reservoir_in"],
                 minus=["salt_crystallised", "salt_pnm_precipitate", "salt_pnm_cap",
                        "salt_pore_inventory"], rtol=1e-10)
    # the network's liquid: what crosses the reservoir throats leaves the pores or stays in them
    # as denser brine (pnm.network_step).  The inventory side -- liquid_stored + the nucleation
    # jumps = V (rho_end - rho_0) summed over the pores -- is checked from the states at the end.
    led.register("liquid_network", plus=["liquid_in"], minus=["liquid_out", "liquid_stored"],
                 rtol=1e-10)
    # what the growth law took out of the pores must equal the crystal it laid down plus the
    # salt still dissolved in the brine that fills the new pore space
    # salt_crust_front: what wet faces evaporating into the crust put into the sealed voxel in
    # front of them (front tracking); salt_base: what base evaporation through the crust puts
    # into the z = 0 voxel.  Both are drawn from the pores in the step they are placed.
    # salt_unsourced on the PLUS side: salt a cluster placed while no pore lay under any of
    # its columns (floating over inert substrate).  It was supplied by nothing the network
    # knows, so it is booked here by name rather than dropped (see pores_feeding).
    led.register("salt_network", plus=["salt_supplied", "salt_unsourced"],
                 minus=["salt_grown", "salt_brine_in_pores", "salt_front_external",
                        "salt_crust_front", "salt_base"], rtol=1e-10)
    # INDEPENDENT check (added 2026-09-06): the salt evaporation hands over, measured from the
    # vapour field, must equal the salt actually placed.  The left side comes from the growth law
    # and the right side from the flux, so this catches a share that does not sum to one -- which
    # is exactly how the missing (1 - F) survived until now.  It excludes salt_nucleated (seed
    # geometry is prescribed, a named gap) and salt_brine_in_pores (dissolved, carried in by the
    # liquid, not produced by evaporation).
    # The salt that crystallises out of a covered pore's brine is ALL the structure gets, the seed
    # included; what has crystallised and is not placed yet sits in the interface buffer.
    led.register("salt_interface", plus=["salt_crystallised"],
                 minus=["salt_supplied", "salt_nucleated", "salt_interface_buffer"], rtol=1e-10)
    led.register("salt_precipitation",
                 plus=["salt_grown", "salt_front_external", "salt_crust_front", "salt_base"],
                 minus=["salt_released_by_evaporation", "salt_released_by_base"], rtol=1e-9)

    m = brine.salt_mass(C0, layer.V, T=T_liquid)   # dissolved NaCl per pore [kg], brine density of C0
    m_pores_start = float(m.sum())
    buffer = np.zeros(layer.N)              # crystallised out of the pore's brine, not placed yet [kg]
    tau_release = None
    covered = np.zeros(layer.n_surface, bool); t_nuc = np.full(layer.n_surface, np.nan)
    frozen = {}; h_freeze = {}
    t = 0.0; hist = []; pore_C = []; M_liquid_in = []; M_ref = None; tic = time.time()
    fld = None; last_sig = None; n_steps = 0
    snaps = []; next_snap = 0.0                       # (t_h, state) every snap_hours
    next_record = 0.0                                 # [h] next hist row when record_hours is set
    n_supply_iter = 0
    n_tiny = 0; last_dt = np.nan; last_limiter = -1       # Zeno guard + step-line diagnostics
    state_at_freeze = None; t_freeze = np.nan          # state and time of the first supply failure
    n_dense_born = 0                                   # voxels born dense under layer_av
    m_cell = np.zeros(g.shape)
    stop_reason = "t_end"

    # Salt accumulated in each candidate cell before it flips [kg].  A cell fills over several
    # steps and the local F changes between them, so its birth porosity is read off the salt it
    # actually holds, eps_birth = 1 - m_cell / (rho_cr V), rather than from one step's F.  That
    # makes the structure's solid mass identically the salt the fluxes handed over, whatever the
    # engine does with the overflow volume (2026-09-07).

    per_pore = openings.per_pore

    def bottom_condition():
        """Dirichlet patches of the pore openings at a_w(C_i).

        A bare opening evaporates, and so does an opening under a SEALED z = 0 voxel: the brine
        there (capped at C_sat) evaporates through the crust, and solve_vapor gives it a base
        flux only where that voxel is active, i.e. conducting -- a voxel whose conductance has
        closed is a wall and contributes nothing.  An opening under a WET z = 0 voxel is
        liquid-connected to the structure and evaporates through the structure's own faces,
        never through the bottom."""
        Ci = brine.mass_fraction_from_salt_mass(m[surf], layer.V[surf], T=T_liquid)
        X_open = np.array([a_w(min(c, C_SAT)) * air.W_sat for c in Ci])
        return openings.bottom_condition(X_open, blocked=eng.state[:, :, 0] == SALT_WET)

    def signature(bmask, bval):
        S = eng.state_w
        # Everything the vapour solve reads, compared element by element (2026-09-15; it used to
        # compare counts and sums, which is complete only while every change is monotone --
        # GAS -> WET -> SEALED, D_eff falling -- and a coincidence of two opening values could
        # still hide).  The window height is part of it: a field solved on a shorter window must
        # not be reused after the window grew.  The sealed voxels' D_eff is part of it: base salt
        # and front tracking keep lowering a SEALED voxel's conductance without changing its
        # state (2-D driver, found 2026-09-07).  No tolerance: any change re-solves.
        d_eff = crust.D_eff[..., :eng.grid_w.n_z][S == SALT_SEALED]
        return (eng.grid_w.n_z, S.tobytes(), bmask.tobytes(), bval.tobytes(), d_eff.tobytes())

    def solve_field(bmask, bval):
        gw = eng.grid_w; nzw = gw.n_z
        porous = None
        pa = crust.porous_argument(eps_min, air.D_v)
        if pa is not None:
            porous = (pa[0][..., :nzw], pa[1][..., :nzw])
        x0 = None
        if fld is not None:                       # previous field, window may have grown
            x0 = np.full(gw.shape, X_top); k = min(fld.X.shape[-1], nzw); x0[..., :k] = fld.X[..., :k]
        return solve_vapor_adaptive(gw, eng.state_w, props, bottom=(bmask, bval), solver="cg",
                                    x0=x0, rtol=VAPOR_RTOL,
                                    top_external_length=L_ext + (nz - nzw) * a, porous=porous)

    # Active window (README E1 section 9.11): all per-step lattice work runs on the lowest
    # k_top + 1 + margin layers.  The margin holds the adaptive vapour cut (one lateral period
    # + 4) plus room for the flips of this step; k_top can only rise, by at most one layer per
    # step, so the window is recomputed at the top of every step.
    margin_w = n_lat + 8
    solid_full = eng.solid_mask_full()
    k_top = int(np.max(np.nonzero(solid_full.any(axis=(0, 1)))[0])) if solid_full.any() else -1

    while t < t_end:
        gw = eng.set_window(min(nz, k_top + 1 + margin_w))
        S = eng.state_w; zidx_w = zidx[..., :gw.n_z]; zc_w = zc_m[..., :gw.n_z]
        salt = eng.salt_mask(); solid = eng.solid_mask()
        covered[np.unique(pore_id_map[(pore_id_map >= 0) & solid[:, :, 0]])] = True
        # ---- 1-2 gas-phase domain (GAS + conducting crust) and its solve.  Refresh whenever the
        # geometry or a state changed; a "every N numerical steps" fallback is forbidden,
        # it would make the physical refresh interval depend on dt (plan/17 P1a item 5).
        bmask, bval = bottom_condition()
        sig = signature(bmask, bval)
        if fld is None or sig != last_sig:
            fld = solve_field(bmask, bval); last_sig = sig
        n_steps += 1
        if n_steps % 100 == 0:
            Cp = brine.mass_fraction_from_salt_mass(m, layer.V, T=T_liquid)
            print(f"    step {n_steps} t={t/3600:.2f} h theta={float(solid[:, :, 0].mean()):.2f} "
                  f"h_max={mo.height_stats(solid, gw)['h_max']*1e3:.2f} mm "
                  f"dt={last_dt:.1f}s lim={LIMITERS[last_limiter] if last_limiter >= 0 else '-'} "
                  f"C_pore={Cp.min():.3f}..{Cp.max():.3f} ({time.time()-tic:.0f} s)",
                  flush=True)
        if M_ref is None:
            M_ref = fld.M_total
        base = fld.base_flux
        # ---- clusters of wet salt and the pores that feed them
        lab, n_c = mo.label_clusters(salt, gw)
        lab_bottom = lab[:, :, 0]
        sel = (pore_id_map >= 0) & (lab_bottom > 0)
        pairs = np.unique(np.stack([pore_id_map[sel], lab_bottom[sel]], 1), axis=0) if sel.any() \
            else np.zeros((0, 2), int)
        cl_pores = {c: [] for c in range(1, n_c + 1)}
        for i, c in pairs:
            cl_pores[int(c)].append(int(i))
        # Pores under the cluster's COLUMN FOOTPRINT, whatever the state of the z = 0 voxel.  With
        # a conducting crust the sealing front sweeps DOWN from the top and base salt seals the
        # z = 0 voxels, so a wet cluster can end up standing on SEALED voxels: it has no wet cell
        # at z = 0, hence no entry in cl_pores, and is frozen (h_new = 0, as in 2-D: no growth
        # through a sealed base).  It still evaporates (plan/16 M2) and that brine must come from
        # somewhere -- through the sealed base, from the pores beneath its columns.  cl_pores
        # keeps deciding growth and supply; cl_pores_any decides where a frozen cluster's WATER
        # and SALT are drawn from.  Before 2026-09-12 such a cluster's salt was placed but charged
        # to nobody, and the salt_network balance raised (-3e-3 at lambda 32, 2014 conditions).
        cl_pores_any = {c: [] for c in range(1, n_c + 1)}
        if n_c:
            ii, jj, kk = np.nonzero(lab > 0)                    # every wet voxel: (column, label)
            pid_v = pore_id_map[ii, jj]; lab_v = lab[ii, jj, kk]
            on_pore = pid_v >= 0
            for i, c in set(zip(pid_v[on_pore].tolist(), lab_v[on_pore].tolist())):
                cl_pores_any[int(c)].append(int(i))
            for c in cl_pores_any:
                cl_pores_any[c].sort()
        def pores_feeding(c):
            """Pores charged for cluster c's water and salt: its wet-base pores if it has any,
            otherwise the pores under its columns; [] only for a cluster floating over inert
            substrate, whose draw is then booked as salt_unsourced (a named gap, never dropped)."""
            return cl_pores[c] if cl_pores[c] else cl_pores_any[c]

        # How a cluster's draw is shared between the pores that feed it: in proportion to the
        # OPENING AREA A_open, so that every opening of one cluster passes the same Darcy
        # velocity U = draw / A_open -- the quantity the supply check itself is written in
        # (U_of).  Before 2026-09-16 the draw was split equally, which is only the same thing
        # when all openings have the same size; the dgf radii are 270 +- 10 um even at cv = 0,
        # and a pore-size distribution (cv > 0) spreads A_open ~ r^2 much further.  Water and
        # salt use the same weights: the salt leaves the pore with the liquid drawn from it.
        # cl_pores holds OPENING numbers (0 .. n_surface - 1, the order of surf); the network is
        # indexed by pore, and surf[k] is opening k's pore.  On the generated lattice the two
        # coincide, in a dgf column they do not (its surface pore is the LAST one): until
        # 2026-09-19 the draw of a cluster on a dgf substrate was charged to the bottom pore.
        feed = {}
        for c in range(1, n_c + 1):
            pf = surf[np.asarray(pores_feeding(c), int)]
            if pf.size:
                w = layer.A_open[pf]
                feed[c] = (pf, w / w.sum())
        F_cluster = np.zeros(n_c + 1)
        for c in range(1, n_c + 1):
            if cl_pores[c]:
                F_cluster[c] = float(np.mean(F[cl_pores[c]]))
        # WHERE THE STRUCTURE'S SALT COMES FROM (2026-09-19).  A pore under salt is in contact with
        # crystal, so its brine cannot stay above C_sat: what the network carries into it beyond
        # saturation crystallises (block 7) and THAT is the structure's salt -- all of it, and
        # nothing else.  The pore is not charged by the growth law any more; it holds C_sat for as
        # long as salt arrives faster than it diffuses back down, and falls below (no growth)
        # when it does not, which is where a real crystal would start to dissolve -- a branch
        # this model does not have.  The faces share a pore's salt in proportion to the water
        # they evaporate (Prat 2026: V_n ~ j), i.e. they hand over beta_i j with
        #     beta_i = release_i / water drawn through pore i,
        # in place of eq (29)'s C_eq / (1 - C_eq); in steady state beta -> C0 / (1 - C0), which
        # is what salt conservation demands of any closure.  The buffer is released over
        # tau_release = RELEASE_TIMES pore-exchange times (numerics [N]: the salt of step n is
        # placed in the steps after it, so that the rates are known before dt is).
        exposed = []
        W_c = np.zeros(n_c + 1)                      # water evaporating on cluster c's salt [kg/s]
        for ax, s_, mask in eng.exposed_faces():
            nb = np.where(mask, gw.shift(lab, ax, s_, 0), 0)
            j = fld.face_flux(ax, s_)
            exposed.append((ax, s_, mask, nb, j))
            W_c += np.bincount(nb[mask], weights=np.maximum(j, 0.0)[mask] * A_cell, minlength=n_c + 1)
        into_crust = []
        sealed_w = S == SALT_SEALED
        for ax, s_ in gw.face_directions():
            mc = sealed_w & gw.shift(salt, ax, s_, False) & fld.active_gas
            if not mc.any():
                continue
            j = fld.face_flux(ax, s_)
            nbc = np.where(mc, gw.shift(lab, ax, s_, 0), 0)
            into_crust.append((ax, s_, mc, nbc, j))
            W_c += np.bincount(nbc[mc], weights=np.maximum(j, 0.0)[mc] * A_cell, minlength=n_c + 1)
        base_sealed = np.where(eng.state[:, :, 0] == SALT_SEALED, np.maximum(base, 0.0), 0.0)
        W_pore = np.zeros(layer.N)                   # water drawn through pore i by depositing surfaces
        W_pore[surf] = per_pore(base_sealed) * A_cell
        for c, (pf, w) in feed.items():
            W_pore[pf] += w * W_c[c]
        # below one NaCl formula unit there is nothing to release (and a geometrically emptying
        # buffer would otherwise run the rates into denormal numbers)
        release = (np.where(buffer > SALT_QUANTUM, buffer, 0.0) / tau_release
                   if tau_release is not None else np.zeros(layer.N))                 # [kg/s]
        # ... but never more than BETA_SAT: a face precipitates what the saturated brine it
        # evaporates leaves behind, which is Prat (2026) eq (29) itself, so eq (29) is the UPPER
        # BOUND of this closure and is reached whenever the network delivers enough.  What it
        # delivers beyond that (a reservoir above saturation does, by 16 % at C0 = 0.30) stays in
        # the buffer: crystal at the pore-structure contact, conserved and named, with no volume
        # of its own [A] -- the role the saturation cap used to play by discarding it.  Without
        # the bound a pore whose consumers have almost stopped evaporating would hand its whole
        # release to them (beta -> infinity) and overfill a sealed base voxel.
        # beta is the salt that becomes CRYSTAL per kg of water, as in eq (29).  A growing face
        # draws gamma times as much again, dissolved in the saturated brine that fills its new
        # pores, so the release is budgeted for (1 + gamma); a frozen face needs less than its
        # budget and the difference simply stays in the buffer.
        gamma_pore = np.zeros(layer.N)
        eps_b_pore = 1.0 - (1.0 - eps_s0) / F
        gamma_pore[surf] = eps_b_pore * RHO_L * C_SAT * F / (RHO_CR * (1 - eps_s0))
        beta_pore = np.minimum(np.divide(release, W_pore * (1.0 + gamma_pore), out=np.zeros(layer.N),
                                         where=W_pore > 0.0), BETA_SAT)
        # A cluster no pore feeds (wet salt the sealing front has cut off, standing over inert
        # substrate) keeps BETA_SAT, as it always has: neither its water nor its salt has a source
        # in this model -- it would really dry out, a state that does not exist here -- and the
        # salt is booked as salt_unsourced.  With no salt at all it would never seal and would
        # evaporate for ever, which is the worse of the two errors.
        beta_c = np.full(n_c + 1, BETA_SAT)
        salt_share = {}
        for c, (pf, w) in feed.items():
            beta_c[c] = float(np.dot(w, beta_pore[pf]))
            salt_share[c] = w * beta_pore[pf] / beta_c[c] if beta_c[c] > 0.0 else w
        zmax = ndimage.maximum(zidx_w, lab, index=np.arange(1, n_c + 1)) if n_c else np.zeros(0)

        if layer_av:
            # B': F of the LAYER a voxel at height z belongs to, Da*_n = Da*_1 r^(n-1) (the
            # specific surface follows the crystal size), clamped at F_min: dense birth
            z_centres = (np.arange(gw.n_z) + 0.5) * a
            gen_grid = closure.generation(z_centres[None, None, :], d_sub_map[..., None])
            F_grid = np.asarray(prat.F(Da_map[..., None] * closure.layer_da_factor(gen_grid)), float)
            F_grid = np.maximum(F_grid, F_min)
        else:
            F_grid = np.broadcast_to(F_map[..., None], gw.shape)

        # ---- 3 external wet faces
        # Salt bookkeeping of an evaporating face (corrected 2026-09-06,
        # ideas/计算结果系统评估_2026-09-06.md section 3).  A face evaporating j hands over
        # beta * j * A of salt (beta of the feeding brine).  The Prat (2026) growth law converts the fraction F(Da*)
        # of it into interface advance (eq 29, the kinetic efficiency); the rest is delivered
        # by the same evaporation and must go somewhere.  One rule covers the three cases:
        #     into new structure       = F * [face may advance] * beta * j * A
        #     into the pores behind it = (1 - F * [face may advance]) * beta * j * A
        # so an advancing face leaves the (1 - F) share and a frozen face leaves all of it
        # (plan/16 M2).  Putting the leftover into the pores is an ASSUMPTION [A]; it is the
        # only destination that conserves mass and it is the pore clogging that Prat's own
        # discussion says should occur (knowledge/prat2026.md, open question 2).  Until this date
        # the driver used F in all three places and dropped the remainder, so the precipitation
        # balance failed by 1 - F: 5.3 % for 1 mm beads, 24.1 % at 34 um.
        Vn_face = {}; nb_faces = []
        # ONE liquid demand per cluster, the brine its faces lose: water plus the salt that water
        # carried, of every wet face (Q_draw), plus, where a face grows, the brine that fills the
        # new pores (Q_fill_face).  The pore network loses exactly this mass, and the supply gate
        # is asked about the Darcy velocity of exactly this brine.  For a growing face it IS Prat
        # (2026) eq (33), U = U_evap + [eps + rho_cr / rho (1 - eps)] dh/dt, because the crystal of
        # the new material is the salt the face hands over (rho_cr (1 - eps_b) V_n = beta j).
        # Eq (33) is written for an advancing front and has no term for salt deposited into
        # EXISTING pore space -- a frozen face, a face evaporating into the crust -- yet that salt
        # arrives dissolved in brine that flows through the same structure.  Until 2026-09-21 the
        # gate kept eq (33) as written and saw beta j / RHO_L less brine there than the pores lost
        # (user's decision of that date: the gate counts it; the measured effect is in
        # tests/reference_cases.py).  Base evaporation under a sealed z = 0 voxel is in neither:
        # that water and its salt leave the pore through the crust, not through wet structure.
        eps_b_grid = 1.0 - (1.0 - eps_s0) / F_grid
        Q_draw = np.zeros(n_c + 1)          # water and salt of EVERY wet face, frozen or not
        Q_fill_face = {}                    # pore filling, only where a face grows
        dm_ext_rate = np.zeros(gw.shape)     # frozen EXTERNAL faces (block 4b) [kg/s]
        W_salt = 0.0                         # water evaporating at WET SALT interfaces [kg/s]
        salt_face = {}                       # salt a face hands over as crystal [kg/m^2/s]
        n_exposed = 0
        for ax, s_, mask, nb, j in exposed:
            n_exposed += int(mask.sum())
            # This growth closure has no dissolution branch. Negative vapour
            # flux (condensation or solver round-off) must not subtract solid
            # in the ledger while the geometry's positive-rate mask drops it.
            # F of the COLUMN the new voxel will occupy, so that the salt the face hands over,
            # beta j A, is exactly rho_cr (1 - eps_birth) times the volume it creates.
            salt_face[(ax, s_)] = salt_f = beta_c[nb] * np.maximum(j, 0.0)
            Vn_face[(ax, s_)] = salt_f * F_grid / (RHO_CR * (1 - eps_s0))
            nb_faces.append((ax, s_, mask, nb))
            Q_draw += np.bincount(nb[mask], weights=((j + salt_f) / RHO_L * A_cell)[mask],
                                  minlength=n_c + 1)
            Q_fill_face[(ax, s_)] = eps_b_grid * Vn_face[(ax, s_)] * A_cell
            # Water measured on the same faces and with the same max(j, 0) as the salt, so the
            # precipitation balance tests the SHARES (do they sum to one?) rather than the flux.
            W_salt += float(np.maximum(j, 0.0)[mask].sum()) * A_cell
        # Faces of wet salt that evaporate INTO the conducting crust instead of into gas
        # (2-D driver 2026-09-07).  solve_vapor already computes their flux (its
        # interface test is `active & wet-salt neighbour`, and a conducting sealed voxel is
        # active), but exposed_faces() only knows gas neighbours, so they are picked up here or
        # the evaporation the field reports in M_total would appear in no salt balance.  They
        # draw liquid and release salt like any other evaporating face; they cannot ADVANCE (the
        # voxel in front is solid), so the whole beta * j goes into pore space -- which pore
        # space is decided in block 4b (front tracking).
        crust_faces = []
        for ax, s_, mc, nbc, j in into_crust:
            salt_f = beta_c[nbc] * np.maximum(j, 0.0)
            Q_draw += np.bincount(nbc[mc], weights=((j + salt_f) / RHO_L * A_cell)[mc],
                                  minlength=n_c + 1)
            W_salt += float(np.maximum(j, 0.0)[mc].sum()) * A_cell
            crust_faces.append((ax, s_, mc, nbc, salt_f))

        # ---- 4 supply fixed point (plan/17 P5):
        #        active growth faces -> cluster demand Q_k -> pore sinks q -> network pressure p
        #        -> footprint suction s_foot -> h_lim -> active growth faces
        # The evaporation demand is unconditional (a frozen face keeps evaporating, plan/16 M2);
        # only the pore-filling demand of GROWING faces is relaxed when the top freezes, which is
        # what closes the loop.  Two nested levels, because a plain Picard sweep over the active
        # SET oscillates (the set is a step function of h_lim):
        #   inner, per cluster and at frozen p:  bisect Phi(h) = h_lim(U(h)) - h.  U(h) increases
        #     with h (more growing faces) and h_lim(U) decreases with U, so Phi is monotone and
        #     has exactly one root -- the same construction as the 2-D prototype (E2).
        #   outer, over the network:  q -> p -> s_foot, under-relaxed, until the active set stops
        #     changing and |dh_lim| < 1e-3 a.
        A_open_c = {c: float(layer.A_open[surf[cl_pores[c]]].sum()) for c in range(1, n_c + 1) if cl_pores[c]}
        dsub_c = {c: float(np.mean(d_sub[cl_pores[c]])) for c in range(1, n_c + 1) if cl_pores[c]}
        h_top_c = {c: (zmax[c - 1] + 1) * a for c in range(1, n_c + 1)}
        # per-cluster face table (face height, cumulative pore-filling demand up to that height)
        if n_c:
            nb_all = np.concatenate([nb[mask] for _, _, mask, nb in nb_faces])
            zc_all = np.concatenate([zc_w[mask] for _, _, mask, _ in nb_faces])
            qf_all = np.concatenate([Q_fill_face[(ax, s_)][mask] for ax, s_, mask, _ in nb_faces])
        z_faces_c = {}; fill_cum = {}
        for c in range(1, n_c + 1):
            idx = np.where(nb_all == c)[0] if n_c else np.zeros(0, int)
            o = np.argsort(zc_all[idx]) if len(idx) else np.zeros(0, int)
            z_faces_c[c] = zc_all[idx][o]
            fill_cum[c] = np.concatenate([[0.0], np.cumsum(qf_all[idx][o])])

        def Q_of(c, h):
            """Brine cluster c draws [m^3/s] when its faces below h grow and those above are
            frozen but still evaporating."""
            return Q_draw[c] + fill_cum[c][int(np.searchsorted(z_faces_c[c], h))]

        def U_of(c, h):
            """The same demand per unit footprint opening: the Darcy velocity the gate tests."""
            return Q_of(c, h) / A_open_c[c]

        def reach(c, h, s_foot):
            return first_failure_height(U_of(c, h), dsub_c[c], s_foot, c_k, h_top_c[c], structure)

        h_lim = np.full(n_c + 1, np.inf)
        # The network balance is on liquid MASS (pnm.solve_pressure, 2026-09-17), and the MASS a
        # pore loses does not depend on its density: a bare opening loses the water that
        # evaporates, m_dot, and a feeding pore loses RHO_L Q_k, because every cluster demand
        # is a volume of the structure's saturated brine (m_dot / RHO_L, and the pore filling).
        # q is that mass as a volume of the PORE's brine, so that q * rho_pore hands
        # solve_pressure the mass itself.  Until 2026-09-19 q was built with RHO_L throughout,
        # which left the sink wrong by rho_pore / RHO_L: 0.08 % at C_sat, but -4.6 % at C = 0.20
        # (the coupling check against DuMuX, tests/test_case_dumux1d.py) and +4 % at S_NUC C_sat.
        rho_pore = brine.density(brine.mass_fraction_from_salt_mass(m, layer.V, T=T_liquid),
                                 T=T_liquid)
        rho_res = float(brine.density(C0, T=T_liquid))
        p = np.zeros(layer.N); q = np.zeros(layer.N); Q_lat = Q_res = None
        q_base_pore = np.zeros(layer.N)                                     # [m^3/s] per pore
        # Under a sealed z = 0 voxel the water leaves through the crust and its salt precipitates
        # into that voxel (block 4d), so there the pore loses brine, water AND salt; a bare
        # opening loses the water alone and keeps the salt, which is its enrichment.
        beta_col = np.where(pore_id_map >= 0, beta_pore[surf][np.maximum(pore_id_map, 0)], 0.0)
        base_salt_flux = beta_col * base_sealed                                      # [kg/m^2/s]
        q_base_pore[surf] = per_pore(np.maximum(base, 0.0) + base_salt_flux) * A_cell / rho_pore[surf]
        if precip is not None:                             # continuity: the solid displaced liquid
            q_base_pore -= precip.last_volume_rate
        omega = 0.5                       # under-relaxation of the network pressure [numerics]
        active = {(ax, s_): mask.copy() for ax, s_, mask, _ in nb_faces}
        dh = 0.0
        for n_supply_iter in range(1, SUPPLY_ITER_MAX + 1):
            h_new = np.full(n_c + 1, np.inf)
            for c in range(1, n_c + 1):
                if not cl_pores[c]:
                    h_new[c] = 0.0                 # cut from the substrate: nothing can grow
                    continue
                # suction at the cluster's base, averaged over its wet-base openings WEIGHTED BY
                # OPENING AREA (2026-09-16): with U uniform over the cluster's openings, the
                # liquid entering the structure is the area-weighted mixture of those pores, so
                # the same weights as the draw (feed) apply here.  It was an unweighted mean.
                a_c = layer.A_open[surf[cl_pores[c]]]
                s_foot = -float(np.dot(a_c, p[surf[cl_pores[c]]]) / a_c.sum())
                if not np.isfinite(reach(c, h_top_c[c], s_foot)):
                    continue                       # even the full growth demand is supplied
                lo, hi = 0.0, h_top_c[c]
                hl_lo = reach(c, lo, s_foot)
                if np.isfinite(hl_lo) and hl_lo <= lo:
                    h_new[c] = 0.0; continue
                for _ in range(40):                # Phi(lo) >= 0 > Phi(hi), monotone
                    mid = 0.5 * (lo + hi)
                    hl = reach(c, mid, s_foot)
                    if np.isfinite(hl) and hl < mid:
                        hi = mid
                    else:
                        lo = mid
                h_new[c] = 0.5 * (lo + hi)
            # demand implied by h_new -> pore sinks -> network pressure
            Q_k = Q_draw.copy()
            for c in range(1, n_c + 1):
                if cl_pores[c]:
                    Q_k[c] = Q_of(c, h_new[c])
            # base evaporation is a sink of the pore beneath it: on a bare opening and on an
            # opening under a sealed z = 0 voxel.  fld.base_flux is zero everywhere else.
            q = q_base_pore.copy()
            for c, (pf, w) in feed.items():
                q[pf] += Q_k[c] * w * RHO_L / rho_pore[pf]
            p_new, Q_lat, Q_res = solve_pressure(layer, q * rho_pore, g_lat, g_res,
                                                 rho_pore, rho_res)
            changed = False
            for ax, s_, mask, nb in nb_faces:
                na = mask & (zc_w < h_new[nb])
                changed |= not np.array_equal(na, active[(ax, s_)])
                active[(ax, s_)] = na
            capped = lambda arr: np.minimum(np.nan_to_num(arr, posinf=nz * a), nz * a)
            dh = float(np.max(np.abs(capped(h_new) - capped(h_lim)))) if n_c else 0.0
            h_lim = h_new
            if not changed and dh < 1e-3 * a:
                p = p_new
                break
            p = (1.0 - omega) * p + omega * p_new
        else:
            np.savez_compressed(Path.cwd() / f"supply_fixedpoint_failure_{tag or 'run'}.npz",
                                state=eng.state, h_lim=h_lim, q=q, p=p, t=t)
            raise RuntimeError(
                f"supply fixed point did not converge in {SUPPLY_ITER_MAX} iterations at "
                f"t = {t/3600:.3f} h (n_clusters = {n_c}, max |dh| = {dh:.3e} m = {dh/a:.1f} voxels). "
                "plan/17 P5 forbids continuing silently; diagnostics saved")
        # final network solve at the converged demand, so that Q_lat / Q_res used by the transport
        # step and by dt_net belong to the same q as p
        p, Q_lat, Q_res = solve_pressure(layer, q * rho_pore, g_lat, g_res, rho_pore, rho_res)
        # "Supply-limited" means at least one exposed face of the cluster can no longer be fed,
        # i.e. the reach falls below the highest face the cluster offers.  Testing the reach
        # against the cluster TOP instead would miss the common case of a cluster that has
        # stopped growing while every voxel it already owns is still within reach.
        # Rebuilt every step: a stale entry from a cluster that no longer exists must not be
        # counted in n_frozen.
        frozen = {}
        for c in range(1, n_c + 1):
            z_faces = z_faces_c[c]
            frozen[c] = bool(len(z_faces) and h_lim[c] < z_faces[-1])
        h_freeze = {c: float(h_lim[c]) for c in range(1, n_c + 1)}
        if state_at_freeze is None and any(frozen.values()):
            state_at_freeze = eng.state.copy(); t_freeze = t
        S_released = 0.0                     # crystal the evaporation on wet salt hands over [kg/s]
        for ax, s_, mask, nb in nb_faces:
            Vn_face[(ax, s_)] = np.where(active[(ax, s_)], Vn_face[(ax, s_)], 0.0)
            S_released += float(salt_face[(ax, s_)][mask].sum()) * A_cell
        S_released += sum(float(salt_f[mc].sum()) for _, _, mc, _, salt_f in crust_faces) * A_cell
        # ---- 4b the salt a frozen EXTERNAL face leaves in the pores of its own voxel (section 3
        # rule).  A frozen face keeps evaporating (its demand is in Q_draw) but cannot add a
        # voxel, so all of beta * j goes into the wet voxel behind it (plan/16 M2, channel
        # salt_front_external).  It is drawn from the substrate pores in THIS step, so it enters
        # the network salt balance.  An advancing face's (1 - F) share is carried by the birth
        # porosity of the material it lays down (task B) and never passes through here.
        S_ext_k = np.zeros(n_c + 1); S_ext = 0.0
        for ax, s_, mask, nb in nb_faces:
            frozen_f = mask & ~active[(ax, s_)]
            if not frozen_f.any():
                continue
            # A face that cannot advance hands ALL of its salt to the band behind it.  An
            # advancing face is not here: its (1 - F) share is carried by EPS_B, the porosity the
            # new material is born with (task B), so it never passes through this channel.
            rate_gas = np.where(frozen_f, salt_face[(ax, s_)] * A_cell, 0.0)
            if not rate_gas.any():
                continue
            dm_ext_rate += gw.shift(rate_gas, ax, -s_, 0.0)
            S_ext_k += np.bincount(nb[frozen_f], weights=rate_gas[frozen_f], minlength=n_c + 1)
            S_ext += float(rate_gas[frozen_f].sum())
        # ---- 4c front tracking (2-D driver 2026-09-08).  A face evaporating into the crust
        # feeds the FRONT voxel -- the sealed voxel on the other side of it -- until that voxel's
        # bulk is down to eps_min.  Only a finished voxel passes the salt on to the reaction band
        # of the wet voxel behind, which then starts its own band.  The old rule sent everything
        # to the band at once, so each band closed after delta_C worth of salt while the front
        # advanced a whole voxel: front speed ~ a / delta_C, resistance per height ~ delta_C / a,
        # i.e. the result depended on the mesh.  Bulk deposits go through crust.deposit, band
        # deposits through crust.deposit_skin; the two are never mixed in one array.
        dm_crust_rate = np.zeros(gw.shape)   # into the BULK of the sealed front voxel [kg/s]
        S_crust = 0.0
        for ax, s_, mc, nbc, salt_f in crust_faces:
            rate_crust = np.where(mc, salt_f * A_cell, 0.0)
            if not rate_crust.any():
                continue
            unfinished = mc & (crust.eps_s[..., :gw.n_z] > eps_min + crust.EPS_TOL)
            dm_crust_rate += np.where(unfinished, rate_crust, 0.0)
            finished = np.where(mc & ~unfinished, rate_crust, 0.0)
            dm_ext_rate += gw.shift(finished, ax, -s_, 0.0)
            # both branches draw brine through the wet voxel behind the face in THIS step
            S_ext_k += np.bincount(nbc[mc], weights=rate_crust[mc], minlength=n_c + 1)
            S_crust += float(np.where(unfinished, rate_crust, 0.0).sum())
            S_ext += float(finished.sum())
        dm_ext_rate = np.where(solid, dm_ext_rate, 0.0)
        # ---- 4d base salt.  fld.base_flux is non-zero under a sealed z = 0 voxel only with a
        # conducting crust (bottom_condition); that water is already a pore sink (q_base_pore),
        # and its salt precipitates into the voxel standing on the opening.  A bare opening's
        # base flux leaves its salt IN the pore -- that is the enrichment the network transport
        # carries -- so nothing is placed for it here and no "uncovered" channel is needed.
        dm_base_full = np.zeros(g.shape)
        base_salt_col = base_salt_flux * A_cell
        S_base = float(base_salt_col.sum())
        dm_base_full[..., 0] = base_salt_col
        dm_deposit_rate = dm_ext_rate
        dm_bulk_full = eng.to_full(dm_deposit_rate + dm_crust_rate, 0.0) + dm_base_full

        # The regularisation redistributes the LIQUID along the front, so what it conserves is
        # the SALT that liquid carries, not the volume the salt happens to make.  Smoothing the
        # salt and converting it to a volume at the destination's own birth porosity is what
        # makes the precipitation balance exact when F varies from column to column; smoothing
        # the volume instead moves material between columns of different eps_birth and loses
        # ~0.1 % of the salt (2026-09-07).  Because eps_birth = 1 - (1 - eps_s0)/F, the volume is
        # simply the F-weighted smoothed salt over rho_cr (1 - eps_s0) -- no division needed.
        salt_rate = eng.volume_rate(                                   # [kg/s] per candidate cell
            lambda ax, s_, mask: salt_face[(ax, s_)]
            * (Vn_face[(ax, s_)] > 0.0), smooth_length=ell)
        saltF_rate = eng.volume_rate(                                  # [kg/s], weighted by F
            lambda ax, s_, mask: salt_face[(ax, s_)]
            * (Vn_face[(ax, s_)] > 0.0) * F_grid, smooth_length=ell)
        rate = saltF_rate / (RHO_CR * (1 - eps_s0))
        F_cell = np.divide(saltF_rate, salt_rate, out=np.array(F_grid, dtype=float),
                           where=salt_rate > 0.0)
        eps_b_cell = 1.0 - (1.0 - eps_s0) / F_cell
        crust.F_field[..., :gw.n_z] = F_cell

        # ---- 5 salt the structure takes out of each pore's interface buffer (frozen faces, the
        # crust front and the base included).  NOT a sink of the pore's brine: that salt left the
        # brine when it crystallised (block 7).
        q_salt = np.zeros(layer.N)
        S_k = S_ext_k.copy()
        # What the growth law PLACES, read off the smoothed fields so that it is exactly the salt
        # the faces handed over.
        S_grow = float(salt_rate.sum())
        S_brine = float((rate * eps_b_cell).sum()) * RHO_L * C_SAT
        # Per cluster the salt is exact (the regulariser conserves each source cluster's share
        # separately); the dissolved brine that fills the new pores is allocated pro rata, an
        # ALLOCATION choice whose total is exact.
        salt_c = np.zeros(n_c + 1)
        for ax, s_, mask, nb in nb_faces:
            # EXACTLY the set that salt_rate was built from.  A cluster with no pore under it has
            # F_cluster = 0, hence Vn_face = 0: it evaporates but places nothing, so its faces must
            # not enter this sum either.  Using the raw flux here instead made tot_c exceed S_grow
            # and the pro-rata allocation charged the network too little -- the salt_network
            # balance failed by 1e-5 at lambda = 4 (found 2026-09-07).
            sel_f = Vn_face[(ax, s_)] > 0.0
            if sel_f.any():
                w = (salt_face[(ax, s_)] * A_cell)[sel_f]
                salt_c += np.bincount(nb[sel_f], weights=w, minlength=n_c + 1)
        tot_c = float(salt_c.sum())
        if tot_c > 0.0:
            S_k += salt_c * ((S_grow + S_brine) / tot_c)
        S_unsourced = 0.0
        for c in range(1, n_c + 1):
            if c in feed:
                pf, _ = feed[c]
                q_salt[pf] += S_k[c] * salt_share[c]
            else:
                S_unsourced += float(S_k[c])
        # base salt is charged to the pore whose opening it left, not to a cluster
        if S_base > 0.0:
            q_salt[surf] += per_pore(base_salt_col)

        # ---- 6 time step: every limiter scaled by the same r (plan/17 P1a item 4)
        # "salt": the buffer is released over tau_release, so a step longer than that would place
        # salt that has not crystallised yet.  A consistency bound, not scaled by r.
        dt_salt = tau_release if tau_release is not None else np.inf
        # r scales the ACCURACY limiters (grow, network, salt, cap) only.  The EVENT limiters
        # (seal, end) place the step ON the event; scaling them by r < 1 would approach the event
        # geometrically and never reach it (Zeno steps, found 2026-09-05 in the plan/17
        # verification cases; tests/test_E2_feedback.py pins it).
        # The band path drains the bulk too (deposit_skin removes the same kg from the pore space
        # it adds to the band), so a voxel fed by band, front and base at once is limited on
        # their SUM (2-D driver, found 2026-09-07: limiting them apart overfilled a voxel to
        # eps_s = -5.1e-7).  dt_to_seal(.., eps_min) lands the step on the bulk sealing event;
        # once a fed voxel is at eps_min it goes inactive there and dt_to_seal(.., 0.0) keeps its
        # pore space from going negative.  Armed only when bulk deposits exist, so runs without
        # them keep the step sequence they had.
        bulk_on = dm_crust_rate.any() or dm_base_full.any()
        dt_seal = min(crust.dt_to_close_skin(eng.to_full(dm_deposit_rate, 0.0) / A_cell, eps_min),
                      crust.dt_to_seal(dm_bulk_full, eps_min) if bulk_on else np.inf,
                      crust.dt_to_seal(dm_bulk_full, 0.0) if bulk_on else np.inf)
        # Network accuracy limit: pore volume over the exchange it sees.  With one layer that is
        # the reservoir throat alone, which is what every run up to 2026-09-17 used and what the
        # production records were made with.  A resolved column has interior pores with NO
        # reservoir throat, so there the lateral and vertical throats have to be counted or the
        # limit would be infinite; the two forms are kept apart rather than silently unified,
        # because the second one is tighter and would change the single-layer step sequence.
        exchange = np.abs(Q_res) + d_res
        if layer.n_layers > 1:
            # ADVECTIVE throat exchange only.  The diffusive exchange between two neighbouring
            # pores relaxes in V / d ~ 0.3 s at lambda = 32, and resolving that would set the step
            # even though salt_step is backward Euler and damps it stably; what the step has to
            # follow is how fast the liquid itself is replaced.
            np.add.at(exchange, lat_a, np.abs(Q_lat))
            np.add.at(exchange, lat_b, np.abs(Q_lat))
        limits = (dt_factor * eng.dt_occupancy(rate, f_max=F_MAX),                          # grow
                  dt_factor * 0.25 * float(np.min(layer.V / (exchange + 1e-30))),           # network
                  dt_seal,                                                                   # seal
                  dt_salt, dt_factor * DT_CAP, t_end - t)
        limiter = int(np.argmin(limits))
        dt = float(min(limits))
        last_dt, last_limiter = dt, limiter
        # Zeno guard (2026-09-12).  A limiter that keeps returning a vanishing step -- dt_salt on
        # a drained pore, a band or a bulk that never quite closes -- makes the run spin at a
        # fixed t for ever (crust_5 sat at 16.73 h for 2000 steps).  plan/17 P5 forbids silent
        # continuation; refuse, and save what is needed to see which limiter and which voxel.
        if dt < 1e-3 and (t_end - t) > 1.0:
            n_tiny += 1
            if n_tiny >= 200:
                np.savez_compressed(Path.cwd() / f"zeno_failure_{tag or 'run'}.npz",
                                    state=eng.state, m=m, q_salt=q_salt, limits=np.array(limits, float),
                                    t=t, eps_s=crust.eps_s, m_skin=crust.m_skin, skin_cap=crust.skin_capacity(eps_min),
                                    h_lim=h_lim, dm_bulk=dm_bulk_full, dm_deposit=eng.to_full(dm_deposit_rate, 0.0))
                raise RuntimeError(
                    f"time step collapsed: {n_tiny} consecutive steps below 1 ms at t = {t/3600:.3f} h; "
                    f"limiter = {LIMITERS[limiter]}, limits = {[f'{x:.3e}' for x in limits]}; "
                    f"pore C range {float(brine.mass_fraction_from_salt_mass(m, layer.V, T=T_liquid).min()):.4f}"
                    f"..{float(brine.mass_fraction_from_salt_mass(m, layer.V, T=T_liquid).max()):.4f}; "
                    "diagnostics saved (zeno_failure_*.npz)")
        else:
            n_tiny = 0

        # ---- 7 substrate transport, saturation cap, nucleation
        # A pore under salt touches crystal, so its brine is held at C_sat and what the network
        # carries into it beyond saturation crystallises (pnm.salt_step): that is the structure's
        # salt.  It comes BEFORE in-pore precipitation -- under a structure the crystal to grow on
        # is the structure.
        touching = np.zeros(layer.N, bool); touching[surf] = covered
        # Liquid and salt of the step together (pnm.network_step, 2026-09-21): the pores' own
        # density change is a storage term of the liquid balance, so the flow is re-solved with
        # it and p, Q_lat, Q_res below are the step's, not the supply fixed point's (they differ
        # by the storage, ~1e-3 of the throughput in a transient and zero in a steady state; the
        # supply gate and the dt limiters above saw the flow without it).
        q_mass = q * rho_pore
        step = network_step(layer, m, C0, q_mass, g_lat, g_res, d_lat, d_res, dt,
                            crystal=touching, w_eq=C_SAT, flow=(p, Q_lat, Q_res), T=T_liquid)
        m, crystallised, p, Q_lat, Q_res, storage = step[:6]
        led.add("salt_reservoir_in", step.salt_in)
        m_in = float((np.where(Q_res > 0.0, rho_res, rho_pore) * Q_res).sum())
        led.add("liquid_in", m_in * dt)
        led.add("liquid_out", float(q_mass.sum()) * dt)
        led.add("liquid_stored", float(storage.sum()) * dt)
        buffer += crystallised - q_salt * dt
        led.add("salt_crystallised", float(crystallised.sum()))
        if precip is not None:
            solid_before = precip.m_solid.sum()
            try:
                m = precip.step(m, dt)
            except PoreClogged as clog:
                stop_reason = f"pore clogged ({clog.pores.size} pores at {clog.fraction:.3f})"
                led.add("salt_pnm_precipitate", float(precip.m_solid.sum() - solid_before))
                break
            led.add("salt_pnm_precipitate", float(precip.m_solid.sum() - solid_before))
            if np.any(precip.last_volume_change > 0.0):     # the throats narrowed: re-read them
                g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
                lat_a = np.array([t_[0] for t_ in layer.lat], int)
                lat_b = np.array([t_[1] for t_ in layer.lat], int)
        led.add("salt_supplied", float(q_salt.sum()) * dt)
        led.add("salt_grown", S_grow * dt)
        led.add("salt_brine_in_pores", S_brine * dt)
        led.add("salt_front_external", S_ext * dt)
        led.add("salt_crust_front", S_crust * dt)
        led.add("salt_base", S_base * dt)
        led.add("salt_unsourced", S_unsourced * dt)
        led.add("water_evaporated", fld.M_total * dt)
        led.add("water_evaporated_base", fld.M_base * dt)
        led.add("water_evaporated_on_salt", W_salt * dt)
        led.add("salt_released_by_evaporation", S_released * dt)
        led.add("salt_released_by_base", S_base * dt)
        C = brine.mass_fraction_from_salt_mass(m[surf], layer.V[surf], T=T_liquid)   # nucleation
        # ---- 8 advance the geometry, then give the new voxels their birth properties from the
        # flip mask (never from a porosity test, plan/17 section 3.3)
        m_cell[..., :gw.n_z] += salt_rate * dt
        m_view = m_cell[..., :gw.n_z]
        rep = eng.advance(rate, dt, carry=m_view)
        m_cell[..., :gw.n_z] = m_view
        if rep.flipped.any():
            flipped_full = eng.to_full(rep.flipped, False)
            crust.birth(flipped_full, t + dt)
            # override the closure's porosity by the one the accumulated salt implies
            eps_from_salt = 1.0 - m_cell[flipped_full] / (RHO_CR * g.cell_volume)
            if np.any(eps_from_salt <= 0.0) or np.any(eps_from_salt >= 1.0):
                raise RuntimeError("a voxel was born with a porosity outside (0, 1): the salt it "
                                   "accumulated does not match the volume it was given")
            crust.eps_s[flipped_full] = eps_from_salt
            crust.eps_birth[flipped_full] = eps_from_salt
            if layer_av:
                crust.shrink_at_birth(flipped_full, eps_s0)
            crust.update_properties(flipped_full)
            m_cell[flipped_full] = 0.0
            if layer_av:
                # a voxel born at the clamp holds eps_min of pore space to round-off: it is the
                # dense, dry state from birth.  The tolerance is the accumulated round-off of
                # m_cell, far above EPS_TOL and far below any porosity difference that matters.
                dense = crust.seal(eps_min + DENSE_BIRTH_TOL, t + dt)
                n_dense_born += int(dense.sum())
        led.add("volume_unplaced", rep.volume_banked)
        # without in-pore precipitation a bare pore has to become supersaturated to nucleate;
        # with it the pore can never pass C_sat, so the first crystal IS the nucleation event
        nuc = (~covered & (precip.m_solid[surf] > 0.0) if precip is not None
               else ~covered & (C >= S_NUC * C_SAT))
        for i in np.where(nuc)[0]:
            # Nucleation is local: a first crystal occupies one deterministic
            # cell at the pore centre, not the whole coarse PNM opening.  The
            # remaining opening stays a bare liquid boundary and subsequent
            # face-by-face growth decides whether the salt spreads across it.
            # ``salt_nucleated`` remains a named prescribed-seed mass gap;
            # this geometry correction does not claim to close that separate
            # sub-grid nucleation-mass problem.
            eng.state[:, :, 0][openings.nucleation_seed(i)] = SALT_WET
            covered[i] = True; t_nuc[i] = t
            # the supersaturation the pore built up IS the first crystal's salt
            m_sat = brine.salt_mass(C_SAT, layer.V[surf[i]], T=T_liquid)
            led.add("liquid_nucleation_jump", float(
                layer.V[surf[i]] * (brine.density(C_SAT, T=T_liquid)
                                    - brine.density(C[i], T=T_liquid))))
            buffer[surf[i]] += m[surf[i]] - m_sat
            led.add("salt_crystallised", float(m[surf[i]] - m_sat))
            m[surf[i]] = m_sat
        if nuc.any():
            seeded = (eng.state == SALT_WET) & ~np.isfinite(crust.t_birth)
            crust.birth(seeded, t + dt)
            # The seed's geometry is prescribed (one cell, interface.nucleation_seed), its salt is
            # not: it comes out of the pore's buffer, which goes negative if the seed holds more
            # than the pore gave and then releases nothing until growth has paid it back.
            seed_mass = (1 - crust.eps_birth[:, :, 0]) * seeded[:, :, 0] * g.cell_volume * RHO_CR
            buffer[surf] -= per_pore(seed_mass)
            if seeded[:, :, 1:].any():
                raise RuntimeError("a nucleation seed was placed outside the bottom layer")
            led.add("salt_nucleated", float(seed_mass.sum()))
        # ---- 9 pore filling by the frozen external faces (block 4b)
        if dm_deposit_rate.any():
            # Into the sub-grid reaction band (task B), not into the whole voxel: the band is
            # delta_C thick, so its closing time is a physical time.  deposit_skin also takes the
            # same mass off the bulk porosity, which keeps salt_solid_mass and the ledger exact.
            crust.deposit_skin(eng.to_full(dm_deposit_rate * dt, 0.0) / A_cell, eps_min)
        # ---- 9b bulk pore filling, not the band: the crust front finishing its own voxel and
        # the base salt.  crust.seal turns a voxel whose bulk reached eps_min into SALT_SEALED
        # (a base voxel is sealed already; its conductance simply keeps falling).
        if dm_crust_rate.any() or dm_base_full.any():
            crust.deposit((eng.to_full(dm_crust_rate, 0.0) + dm_base_full) * dt)
            crust.seal(eps_min, t + dt)
        t += dt
        tau_release = RELEASE_TIMES * 0.25 * float(np.min(layer.V / (exchange + 1e-30)))

        # ---- 10 a voxel whose reaction band has closed becomes SALT_SEALED
        crust.seal_skin(eps_min, t)

        # ---- 11 re-solve the field that goes into the record, so that E(t) belongs to state(t)
        # and not to the state one step earlier (plan/17 P1a item 6).  The next iteration reuses
        # it through the same signature test, so this costs no extra solve.
        salt = eng.salt_mask(); solid = eng.solid_mask()
        covered[np.unique(pore_id_map[(pore_id_map >= 0) & solid[:, :, 0]])] = True
        cols_z = np.nonzero(solid.any(axis=(0, 1)))[0]
        k_top = int(cols_z.max()) if len(cols_z) else -1
        if k_top + 1 + (n_lat + 4) > gw.n_z and gw.n_z < nz:
            # salt grew into the solver margin this step: widen the window before the record solve
            gw = eng.set_window(min(nz, k_top + 1 + margin_w))
            S = eng.state_w; zidx_w = zidx[..., :gw.n_z]; zc_w = zc_m[..., :gw.n_z]
            salt = eng.salt_mask(); solid = eng.solid_mask()
        bmask, bval = bottom_condition()
        sig = signature(bmask, bval)
        if sig != last_sig:
            fld = solve_field(bmask, bval); last_sig = sig
        base = fld.base_flux
        # Source-resolved evaporation diagnostic.  This changes neither the
        # vapour field nor the growth law: it only partitions the conservative
        # M_base + M_salt field that is being recorded at this same time.
        sealed_base = eng.state[:, :, 0] == SALT_SEALED
        bare_base = bmask & (eng.state[:, :, 0] == GAS)
        M_base_bare = float((base[bare_base] * A_cell).sum())
        M_base_crust = float((base[sealed_base] * A_cell).sum())
        M_wet_gas = M_wet_crust = 0.0
        for face_j in fld.face_fluxes.values():
            M_wet_gas += float((face_j[S == GAS] * A_cell).sum())
            M_wet_crust += float((face_j[S == SALT_SEALED] * A_cell).sum())
        if snap_hours is not None and t / 3600 >= next_snap:
            snaps.append((t / 3600, eng.state.copy())); next_snap += snap_hours

        # ---- 12 record
        jt = max(v.max() for v in fld.face_fluxes.values()) if salt.any() else 0.0
        hs = mo.height_stats(solid, gw)
        vols = crust.volumes()
        # blocked column (P7a): the column's EVAPORATING surfaces (wet faces, bare openings) have no
        # gas path to the far field, or the column has no evaporating surface left (fully sealed).
        # "No far-field-connected gas anywhere in the column" is wrong: the air above a sealed roof
        # is connected, so that test never fires once the roof sits below the box top
        # (plan/17 verification case C, 2026-09-05).
        gas_now = S == GAS
        evap_gas = np.zeros(gw.shape, bool)
        for ax, s_ in gw.face_directions():
            evap_gas |= gas_now & gw.shift(salt, ax, s_, False)
        evap_gas[:, :, 0] |= gas_now[:, :, 0] & bmask
        # only columns that carry salt or a pore opening can be blocked; inert substrate between the
        # pores has no evaporating surface and must not be charged with the leak
        has_surface = solid.any(axis=gw.z_axis) | bmask
        blocked = has_surface & ~(evap_gas & fld.active_gas).any(axis=gw.z_axis)
        f_blocked = float(blocked.mean())
        # structure salt balance: what the voxels hold, plus the occupancy that has been
        # injected but not yet flipped, plus what the engine could not place, must equal what
        # the growth law, the front and the nucleation seeds put in (plan/17 P4).
        struct = crust.salt_solid_mass()
        # in-flight solid: eps_birth is per voxel since 2026-09-07, so the partly filled cells are
        # summed with their OWN porosity (from the F field the last step wrote) and only the
        # engine's unplaced bank, which has no location, uses the structure's mean.
        _sf = eng.solid_mask_full()
        mean_eps_b = float(np.mean(crust.eps_birth[_sf])) if _sf.any() else eps_s0
        # in flight = the salt sitting in cells that have not flipped yet, exactly as booked
        pending = float(m_cell.sum()) + eng.bank_carry
        res_salt = ((struct + pending - led["salt_grown"] - led["salt_front_external"]
                     - led["salt_crust_front"] - led["salt_base"]
                     - led["salt_nucleated"]) / max(struct, 1e-300))
        # the network balances liquid MASS since 2026-09-17, so the residual is on mass too: a
        # volumetric one would now show the reservoir/pore density difference (2.6 % at C0 = 0.30)
        # (liquid in) - (liquid the pores lose) - (liquid they keep as denser brine)
        res_water = ((m_in - float(q_mass.sum()) - float(storage.sum()))
                     / max(abs(float(q_mass.sum())), 1e-300))
        box_top = hs["h_max"] >= (nz * a - a)       # salt reached the box top (one voxel below)
        if record_hours is not None and t / 3600 < next_record and t < t_end and not box_top:
            continue                                # between two rows; the last state is always written
        if record_hours is not None:
            next_record = t / 3600 + record_hours
        hist.append((t / 3600, float(solid[:, :, 0].mean()), mo.n_clusters(solid, gw),
                     mo.spanning(solid, gw) if solid.any() else False, hs["h_max"], hs["h_mean"],
                     fld.M_total / M_ref, jt / j_plane, sum(frozen.values()),
                     float(base[bmask].mean()) / j_plane if bmask.any() else np.nan,
                     n_exposed / max(int(solid[:, :, 0].sum()), 1),
                     fld.M_base / M_ref, fld.M_salt / M_ref,
                     M_base_bare / M_ref, M_base_crust / M_ref,
                     M_wet_gas / M_ref, M_wet_crust / M_ref,
                     vols["wet"], vols["sealed"], f_blocked,
                     dt, limiter, float((S == SALT_SEALED).sum()),
                     res_water, res_salt, n_supply_iter))
        pore_C.append(brine.mass_fraction_from_salt_mass(m, layer.V, T=T_liquid))   # at t, like the hist row
        M_liquid_in.append(m_in)
        if box_top:
            stop_reason = "box top"
            break
    hist = np.array(hist, float)
    led.add("salt_interface_buffer", float(buffer.sum()))
    eng.set_window(None)                       # hand back full-grid masks to the caller
    led.add("salt_pore_inventory", float(m.sum()) - m_pores_start)
    led.check("salt_network")
    led.check("salt_pores")
    led.check("liquid_network")
    led.check("salt_interface")
    led.check("salt_precipitation")
    if verbose:
        share = ", ".join(f"{name} {float(np.mean(col(hist, 'limiter') == k)):.2f}"
                          for k, name in enumerate(LIMITERS)) if len(hist) else "-"
        print(f"[vapor lam={lam} C0={C0} c_k={c_k:.0e} Lt={Lt_res*1e3:.1f}mm {tag}] grid {g.shape}, "
              f"a={a*1e6:.2f} um, z_box={nz*a*1e3:.2f} mm, L_ext={L_ext*1e3:.2f} mm, "
              f"RH={X_top/air.W_sat:.2f}, j_plane={j_plane:.2e}, ell={ell*1e6:.1f} um, "
              f"{n_steps} steps, {time.time()-tic:.0f} s | dt set by: {share} | "
              f"salt-network residual {led.relative_residual('salt_network'):.2e}")
        print(f"  reaction-band resolution: delta_C/a = {band_cells:.2f} cells "
              f"(criterion >= {min_band_cells:g}: {'PASS' if band_resolved else 'WARNING'})")
    if snap_hours is not None:
        snaps.append((t / 3600, eng.state.copy()))
    # every input of the run (older records carry more keys: the legacy flags removed 2026-09-14)
    params = dict(lam=lam, C0=C0, c_k=c_k, t_end=t_end, z_top_mm=z_top_mm, Lt_res=Lt_res,
                  layer_av=bool(layer_av), F_min=F_min, n_dense_born=n_dense_born,
                  n_side=n_side, vpp=vpp, delta_eff_mm=delta_eff_mm, rh=rh, cv=cv,
                  pitch_factor=pitch_factor, seed=seed, dt_factor=dt_factor,
                  n_layers=n_layers, dz_pore=dz_pore, pore_precipitation=pore_precipitation,
                  substrate_dgf=("" if substrate_dgf is None else str(substrate_dgf)),
                  substrate_L_box=substrate_L_box, d_sub=grain_size_given,
                  salt_structure=structure.name, T=T, D_v=air.D_v, W_sat=air.W_sat,
                  rho_g=air.rho_g, kr=kr,
                  record_hours=record_hours,
                  ell_over_dsub=ell_over_dsub, eps_min=eps_min, tag=tag,
                  pore_radius=pore_radius,
                  S_nuc=S_NUC, f_max=F_MAX, dt_cap=DT_CAP, vapor_rtol=VAPOR_RTOL,
                  supply_iter_max=SUPPLY_ITER_MAX, min_band_cells=min_band_cells)
    # per-pore histories, one row per hist row: salt mass fraction of EVERY pore (layer.z gives
    # the heights) and the liquid mass the reservoir supplies [kg/s]
    return dict(layer=layer, grid=g, eng=eng, crust=crust, hist=hist,
                pore_C=np.array(pore_C, float).reshape(len(hist), layer.N),
                pore_z=layer.z.copy(), M_liquid_in=np.array(M_liquid_in, float),
                interface_buffer=buffer.copy(),
                hist_columns=HIST_COLUMNS, t_nuc=t_nuc, h_cl=dict(h_freeze), frozen=frozen,
                a=a, z_box=nz * a, L_ext=L_ext, h_freeze=dict(h_freeze), snaps=snaps,
                state_at_freeze=state_at_freeze, t_freeze=t_freeze, n_steps=n_steps,
                stop_reason=stop_reason,
                M_ref=M_ref,
                resolution=dict(reaction_band_m=delta_C, voxel_edge_m=a,
                                reaction_band_cells=band_cells,
                                min_reaction_band_cells=float(min_band_cells),
                                reaction_band_resolved=band_resolved),
                precipitation=(precip.summary() if precip is not None else None),
                ledger=led.report(), closure=closure.describe(), ell=ell, params=params)


def summarize(res, label):
    """Named-column summary (plan/17 P6: no bare hist[:, i] outside HIST)."""
    h = res["hist"]; tn = res["t_nuc"] / 3600; ok = np.isfinite(tn)
    c = lambda name: col(h, name)
    span = c("spanning")
    span_txt = f"{c('t_h')[np.argmax(span > 0)]:.1f} h" if (span > 0).any() else "never"
    print(f"  {label:34s} nucleated {int(ok.sum())}/{len(tn)} (first {np.nanmin(tn) if ok.any() else np.nan:.2f} h) | "
          f"N_c max {int(c('n_clusters').max())} final {int(c('n_clusters')[-1])} | "
          f"theta final {c('theta')[-1]:.2f} | span {span_txt} | "
          f"h_max {c('h_max')[-1]*1e3:.2f} mm, h_mean {c('h_mean')[-1]*1e3:.2f} mm | "
          f"E final {c('E_total')[-1]:.2f} (max {c('E_total').max():.2f}) | "
          f"tip j/j_plane max {c('j_tip_over_plane').max():.1f} | "
          f"bare j/j_plane final {c('j_bare_over_plane')[-1]:.2f} | frozen {int(c('n_frozen')[-1])} | "
          f"t_end {c('t_h')[-1]:.1f} h")
    print(f"  {'':34s} sealed {c('V_sealed')[-1]*1e9:.3f} mm3 "
          f"({int(c('n_sealed_cells')[-1])} cells), f_blocked {c('f_blocked')[-1]:.2f} | "
          f"dt mean {c('dt').mean():.0f} s, supply iters max {int(c('n_supply_iter').max())} | "
          f"residual water {np.abs(c('res_water')).max():.1e}, salt {np.abs(c('res_salt')).max():.1e}")


def case_label(cs):
    return (f"lam={cs['lam']} C0={cs['C0']} c_k={cs['c_k']:.0e} Lt={cs.get('Lt_res',0.7e-3)*1e3:.1f}mm {cs.get('tag','')}"
            + (f" cv={cs['cv']} seed={cs.get('seed',0)}" if cs.get('cv', 0) > 0 else "")
            + (f" pf={cs['pitch_factor']}" if cs.get('pitch_factor', 1.0) != 1.0 else "")
            + f" vpp={cs.get('vpp', 7)}"
            + (f" ell/d_sub={cs['ell_over_dsub']}" if cs.get('ell_over_dsub') else ""))
