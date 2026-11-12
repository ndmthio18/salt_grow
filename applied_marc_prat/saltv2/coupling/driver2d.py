"""E2: minimal 2-D lattice model of the morphology -> evaporation feedback (plan/16).

Regime: evaporation-wicking.  Skeleton follows Veran-Tissoires & Prat (2014) section 6.2 (2-D lattice,
gas-phase Laplace field, wet salt faces at the saturated-solution vapour pressure, dry salt as walls,
hydraulic reach deciding wet/dry), with three replacements:
  * growth speed on every wet exposed face from the Prat (2026) law V_n = G F(Da*) j  (instead of
    adding one cell at the largest flux),
  * hydraulic reach from the hierarchical stack of plan/13 (first_failure_height), instead of a
    Darcy solve in the salt with a single k_s,
  * the sealing rule of plan/16 M2: a wet face that can no longer grow (frozen top, or cluster cut
    from the base) keeps evaporating, and the salt it precipitates fills its own cell; when the cell
    porosity reaches eps_min the cell becomes SEALED (a wall for vapour and liquid).  The sealing
    time follows from the salt mass balance, no new parameter.
  * the sealed crust conducts vapour (2026-09-07): a SEALED voxel conducts through its own closed
    reaction band (CrustField.sealed_conductance) and the brine surface under it keeps evaporating
    and precipitating into its band, so the crust thickens downward instead of the run freezing.
    (Until 2026-09-14 the driver also carried the earlier wall model with the calibrated leak
    phi_vent = 0.07 of 2014 section 6.3; that path is in the checkpoint commit 95c1acb.)
  * the sealing front finishes a voxel before it moves on (2026-09-08).  Closing the band costs
    (eps_birth - eps_min) rho_cr delta_C per unit area, which is the salt a delta_C-THICK layer
    needs; declaring the whole voxel sealed on that payment advances the front by ``a`` for the
    price of ``delta_C``, so the front ran a/delta_C times too fast and left the material
    a/delta_C times too porous -- both grid quantities (README section 11).  The front now keeps
    depositing into the voxel it has sealed until that voxel's BULK reaches eps_min, and only then
    does the voxel behind it start.  delta_C still sets when a voxel first seals (task B, no grid
    length in the onset); the front speed and the swept material's porosity no longer contain one.
  * base salt (2026-09-07): the BARE SUBSTRATE evaporates as well, and the brine feeding it carries
    salt.  Until now that salt was dropped without being booked.  It is now routed into the salt
    voxel standing on the same column (Wijnhorst 2024: once the crust top loses hydraulic contact
    the crystallisation moves to the crust base), and where the column carries no salt it is booked
    as ``salt_base_uncovered`` -- a model with prescribed seeds and no nucleation has nowhere to
    put it, and that gap must be visible in the ledger rather than hidden in a rounding error.
  * the 2-D history is an 8-column array: t_h, E_vap, E_tot, coverage, h_max, f_sealed_top,
    V_sealed, n_frozen_faces.  E_tot == E_vap since the parameterised leak went; the column is
    kept so that older records and their readers keep working.

The bare substrate between structures is a wet plane at a_w(C_sat) (2014's X_ve), seeds are
prescribed (no nucleation model), the box is periodic laterally and the far field sits at
delta_eff = 18 mm through the Robin top boundary.  2-D: spanning is trivial, so only E(t), coverage,
height and sealed fraction are read off, exactly the quantities Fig. 23 of 2014 constrains.
"""
from __future__ import annotations

import time

import numpy as np

from ..growth import morphology as mo
from ..growth.crust import LEGACY_MIXED, BirthClosure, CrustField
from ..growth.geometry import GeometryEngine
from ..growth.grid import GAS, SALT_SEALED, Grid
from ..growth.prat import prat
from ..growth.supply import first_failure_height
from ..growth.vapor import VaporProperties, solve_vapor
from ..materials.air import RHO_D, W_SAT, X_SALT
from ..materials.brine import C_SAT, D_NACL, RHO_L
from ..materials.salt import RHO_CR
from ..pnm import PoreLayer
from .ledger import MassLedger

SEALED = SALT_SEALED

# Values no run has ever varied (2026-09-14: they were run() keywords with these defaults).
U_FACTOR = 2.0     # a 2-D structure exposes half the surface per unit footprint of a 3-D pillar
                   # of the same width (2h/w versus 4h/w), so the liquid demand per footprint
                   # area is scaled by 2 to keep the hydraulic reach comparable with the 3-D
                   # runs (driver3d) at the same c_k [N]
DT_CAP = 900.0     # upper bound on dt [s] [N]
F_MAX = 0.5        # largest occupancy a candidate CELL may gain in one step [N]


def run(d_sub=68e-6, n_seeds=5, cells_per_pitch=5, z_box_mm=5.0, delta_eff_mm=18.0, rh=0.07,
        c_k=3e-3, eps_min=0.02, t_end_h=150.0, snap_hours=6.0, verbose=True,
        dt_factor=1.0, ell_over_dsub=0.0, tortuosity_exponent=1.5, cv=0.0, pore_seed=0,
        seed_fraction=0.77, salt_structure=None):
    """salt_structure: the microstructure of the salt (``growth.crust.SaltStructure``: pore-size
    ratio AND porosity, one record), a case value.  None: LEGACY_MIXED, the pairing of every
    record to 2026-09-20, not self-consistent and kept to reproduce them.

    Numerical controls (plan/17 P1, all of them pure numerics -- no physics):
      dt_factor r  scales the accuracy limiters at once (dt = min(r dt_grow, dt_seal, r DT_CAP,
                   dt_end)), which is what a dt-convergence scan must vary: the growth limit is
                   GeometryEngine.dt_occupancy with F_MAX and changing that alone does not move
                   the step when another limiter binds
      ell_over_dsub  regularisation length ell = ell_over_dsub * d_sub [P1b].  0 disables it and
                   leaves the growth law without a fastest wavelength (not for quantitative use).
                   ell must not be calibrated on the 2014 data.

    tortuosity_exponent  m in D_eff = eps^m (1/D_v + 1/D_K)^-1.  1.5 [L, Bruggeman] is the default
                   and the only value used for results; exposed because the crust is evaluated at
                   eps_min ~ 0.02, where Bruggeman is an extrapolation rather than a measurement,
                   so the plateau has to be reported as a band over m rather than a single number.
    seed_fraction  width of each prescribed seed as a fraction of the pore pitch (default 0.77,
                   the opening fraction of the 5x5 dgf).  It fixes the surface coverage at t = 0,
                   and it is the SAME for every grain size, which means the model starts every run
                   with a continuous salt layer (theta = 0.83 at t = 0, 1.00 within hours) and has
                   no patchy morphology at all: crusty and patchy differ in this model only by
                   whether that layer seals, never by whether it covers.  Eloukabi's mechanism is
                   the opposite -- fine beads give short structures that spread laterally and
                   coalesce, coarse beads give tall isolated mounds with bare substrate between --
                   so the grain-size dependence that mechanism carries cannot appear here.  The
                   parameter exists to measure how much of the missing dependence lives in the
                   coverage; a physical coverage needs the nucleation model (2014 section 5,
                   project module 2), which is still deferred.
    cv, pore_seed  pore-size distribution (plan/17 P9 ported to 2-D, 2026-09-07).  The pores are
                   drawn by ``pnm.PoreLayer`` -- the same dgf generator, the same
                   lognormal, the same opening-area rescaling and clipping as the 3-D driver --
                   and each column inherits the grain size of the pore it stands on.  cv = 0 still
                   gives the dgf's own spread (270 +- 10 um before the lambda scaling), which is
                   why d_sub, F and delta_C are per column at every cv, not only above zero."""
    pitch = 0.7e-3 * d_sub / 1080e-6                    # same substrate scaling as driver3d (lam = 1080 um / d_sub)
    a = pitch / cells_per_pitch
    nx = n_seeds * cells_per_pitch; nz = int(round(z_box_mm * 1e-3 / a))
    g = Grid((nx, nz), a); eng = GeometryEngine(g)
    delta_eff = delta_eff_mm * 1e-3; L_ext = max(delta_eff - nz * a, 0.0)
    X_top = rh * W_SAT
    props = VaporProperties(rho_D=RHO_D, X_wet=X_SALT, X_top=X_top)
    j_plane = RHO_D * (X_SALT - X_top) / delta_eff
    # ---- per-pore grain size.  PoreLayer is the generator the 3-D driver uses; lam is the same
    # substrate scaling (pitch 0.7 mm / lam, r = dgf median / lam), so 2-D and 3-D draw from one
    # distribution.  n_side = n_seeds gives n_seeds^2 pores; the first n_seeds of them are this
    # row's, an unbiased subsample of the same draw.
    lam = 1080e-6 / d_sub
    layer = PoreLayer(lam=lam, n_side=n_seeds, cv=cv, seed=pore_seed)
    d_sub_pore = 4.0 * np.asarray(layer.r[:n_seeds], float)          # [m], one per pore
    pore_of_col = np.arange(n_seeds * cells_per_pitch) // cells_per_pitch
    d_sub_col = d_sub_pore[pore_of_col]                              # [m] per column
    d_sub_mean = float(d_sub_col.mean())
    structure = LEGACY_MIXED if salt_structure is None else salt_structure
    eps_s0 = structure.eps_s             # porosity of the self-replicating structure, read ONCE
    # Ds: the ONE salt diffusivity of the model (materials.brine.D_NACL), not prat2026's own
    # 1.3e-9 default, which belongs to that package's reproduction of the paper (2026-09-21)
    F_col = np.asarray(prat.F(prat.Da_star(eps_s=eps_s0, d_sub=d_sub_col, Ds=D_NACL)), float)
    band_col = np.asarray(prat.boundary_layer_thickness(eps_s=eps_s0, d_sub=d_sub_col, Ds=D_NACL), float)
    eps_b_col = 1.0 - (1.0 - eps_s0) / F_col                          # birth porosity per column
    F = float(F_col.mean())
    Gpref = C_SAT / (RHO_CR * (1 - eps_s0) * (1 - C_SAT))
    BETA_SALT = C_SAT / (1 - C_SAT)      # salt left behind per kg of water evaporated [-]
    # seeds: one cell high, 0.77 pitch wide (opening fraction of the 5x5 dgf), centred on each pitch
    w_seed = max(1, int(round(seed_fraction * cells_per_pitch)))
    for i in range(n_seeds):
        eng.seed(int((i + 0.5) * cells_per_pitch), height_cells=1, width_cells=w_seed)
    # per-voxel salt properties (growth.crust): birth porosity eps_s0, densified at frozen faces.
    # The hierarchy of k_s / P_cth is not used here -- the 2-D prototype takes its hydraulic reach
    # from first_failure_height directly -- but eps_s, the sealing rule and the ledger are shared
    # with the 3-D driver, so they live in one place.
    # Task B (2026-09-06), same closure as the 3-D driver: the (1 - F) share the growth law
    # cannot advance is carried by the birth porosity eps_birth = 1 - (1 - eps_s0)/F, and a face
    # that cannot advance fills the physical reaction band delta_C (Prat 2026, after eq 19)
    # instead of a whole voxel.  Neither contains a grid length.
    delta_C = float(band_col.mean())
    closure = BirthClosure(d_sub=d_sub_mean, structure=structure, c_k=c_k,
                           kinetic_efficiency=F, reaction_band=delta_C,
                           tortuosity_exponent=tortuosity_exponent)
    EPS_B = float(closure.birth_porosity(np.array(1.0)))   # box-mean birth porosity, REPORTING only
    # Every voxel property that depends on the grain size is per column, never from the box mean
    # (plan/17 P9).  F drives the growth rate and the birth porosity; delta_C drives the band
    # capacity, the sealing time and the sealed conductance.  Both are handed to the constructor,
    # because it births the seeds and they must not get the box mean either.
    F_grid = np.broadcast_to(F_col[:, None], g.shape).copy()
    eps_b_grid = np.broadcast_to(eps_b_col[:, None], g.shape).copy()
    crust = CrustField(g, eng, closure, rho_cr=RHO_CR, d_sub_field=d_sub_col,
                       F_field=F_col, band_field=band_col)
    ell = ell_over_dsub * d_sub_mean
    zidx = np.broadcast_to(np.arange(nz)[None, :], g.shape); zc = (zidx + 0.5) * a
    M_ref = solve_vapor(g, np.full(g.shape, GAS, np.int8), props, top_external_length=L_ext).M_total
    def state_signature(wet, sealed):
        """What the vapour field depends on.  Cell counts alone are NOT enough once base salt is
        being placed: it keeps lowering an already-SEALED voxel's porosity, and therefore its
        conductance, without changing any count.  Leaving the field un-refreshed then let a voxel
        that had filled to eps_s = 0 keep drawing base flux, and the next deposit overfilled it
        (found 2026-09-07).  The conductance itself therefore enters the signature -- no
        tolerance, so any change re-solves."""
        return (int(wet.sum()), int(sealed.sum()), float(crust.D_eff[sealed].sum()))

    def vapour_porous():
        """solve_vapor's ``porous`` argument for the CURRENT state.  Rebuilt at every solve
        because a voxel's conductance is fixed only once it is sealed (its band is then full,
        so nothing more is deposited into it) -- but which voxels are sealed keeps changing,
        and that is exactly what the (n_wet, n_sealed) signature below detects."""
        return crust.porous_argument(eps_min)

    # Salt accumulated in each candidate cell before it flips.  A voxel's birth porosity is read
    # off this, not off a box mean or off the F of the last step that touched it (plan/17 P9).
    m_cell = np.zeros(g.shape)
    t = 0.0; hist = []; snaps = []; next_snap = 0.0; fld = None; last_sig = None; tic = time.time()
    seal_events = 0
    led = MassLedger()
    led.register("salt_structure", plus=["salt_supplied"],
                 minus=["salt_grown", "salt_front", "salt_nucleated"], rtol=1e-10)
    # Independent precipitation balance (2026-09-06, same correction as the 3-D driver): the
    # salt evaporation hands over must equal the salt actually placed.  salt_base_uncovered sits
    # on the PLUS side: it is released on a bare column but cannot be placed (no salt stands
    # there), so it closes the balance exactly like a deposit would.
    led.register("salt_precipitation", plus=["salt_grown", "salt_front", "salt_base_uncovered"],
                 minus=["salt_released_by_evaporation", "salt_released_by_base"], rtol=1e-9)
    led.add("salt_nucleated",
            float((1 - eps_b_grid)[eng.solid_mask()].sum()) * g.cell_volume * RHO_CR)
    while t / 3600 < t_end_h:
        wet = eng.salt_mask(); sealed = eng.state == SEALED
        sig = state_signature(wet, sealed)
        # refresh whenever the geometry or a state changed (plan/17 P1a item 5: a "every N
        # numerical steps" fallback would make the physical refresh interval depend on dt)
        if fld is None or sig != last_sig:
            fld = solve_vapor(g, eng.state, props, solver="direct", top_external_length=L_ext,
                              porous=vapour_porous())
            last_sig = sig
        # clusters of wet salt; only those touching the substrate are supplied
        lab, n_c = mo.label_clusters(wet, g)
        touches = np.zeros(n_c + 1, bool)
        touches[np.unique(lab[:, 0])] = True; touches[0] = False
        # growth speed and the two demand terms per face: the brine every wet face loses (water
        # AND the salt it carried, whether the face grows or not), and the pore filling of new
        # salt.  For a growing face the sum is Prat (2026) eq (33): its crystal term
        # rho_cr (1 - eps_b) V_n IS the salt BETA_SALT j.  Until 2026-09-21 a face that cannot
        # advance counted for its water only (same change and reason as in the 3-D driver).
        Vn_face = {}; A_foot = np.bincount(lab[:, 0], minlength=n_c + 1) * g.face_area
        faces = []
        Q_draw = np.zeros(n_c + 1)                     # sum over faces of (j + salt)/rho_l * a (always)
        W_salt = 0.0                                   # water evaporating at wet-salt faces [kg/s]
        grow_terms = []                                # (nb, zc-on-face, growth demand per face) for the fixed point
        for ax, s_, mask in eng.exposed_faces():
            nb = np.where(mask, g.shift(lab, ax, s_, 0), 0)
            j = fld.face_flux(ax, s_)
            Vn = Gpref * F_grid * j          # F at the CANDIDATE cell: that is where the salt lands
            brine = (j + BETA_SALT * np.maximum(j, 0.0)) / RHO_L * g.face_area
            Q_draw += np.bincount(nb[mask], weights=brine[mask], minlength=n_c + 1)
            grow_terms.append((nb[mask], zc[mask], (eps_b_grid * Vn * g.face_area)[mask]))
            Vn_face[(ax, s_)] = Vn; faces.append((ax, s_, mask, nb, j))
            W_salt += float(np.maximum(j, 0.0)[mask].sum()) * g.face_area
        # Faces of wet salt that evaporate into the CONDUCTING CRUST instead of into gas
        # (2026-09-07).  solve_vapor already computes their flux -- its interface
        # test is `active & wet-salt neighbour`, and a conducting sealed voxel is active -- but
        # exposed_faces() only knows gas neighbours, so the driver has to pick them up here or
        # the evaporation the field reports in M_total appears in no salt balance.
        # They draw liquid and release salt like any other evaporating face; they cannot ADVANCE
        # (the voxel in front is solid), so their whole BETA_SALT * j goes into the reaction band
        # of the wet voxel behind, exactly as for a frozen face.
        crust_faces = []
        crust_cells = eng.state == SEALED
        for ax, s_ in g.face_directions():
            m = crust_cells & g.shift(wet, ax, s_, False) & fld.active_gas
            if not m.any():
                continue
            j = fld.face_flux(ax, s_)
            nb = np.where(m, g.shift(lab, ax, s_, 0), 0)
            brine = (j + BETA_SALT * np.maximum(j, 0.0)) / RHO_L * g.face_area
            Q_draw += np.bincount(nb[m], weights=brine[m], minlength=n_c + 1)
            W_salt += float(np.maximum(j, 0.0)[m].sum()) * g.face_area
            crust_faces.append((ax, s_, m, j))
        nb_all = np.concatenate([t[0] for t in grow_terms]); zc_all = np.concatenate([t[1] for t in grow_terms])
        qg_all = np.concatenate([t[2] for t in grow_terms])

        def U_of(c, h):
            """Demand per footprint area when faces of cluster c below h grow and faces above are stationary."""
            sel = (nb_all == c) & (zc_all < h)
            return U_FACTOR * (Q_draw[c] + qg_all[sel].sum()) / A_foot[c]

        # hydraulic reach as a fixed point: h* = h_lim(U(h*)).  U(h) increases with h (more growing
        # faces), h_lim(U) decreases with U, so Phi(h) = h_lim(U(h)) - h is monotone: one root.
        # (Codex review 2026-09-04, finding 1: the reach must be consistent with the demand it implies.)
        # Each cluster is fed by the pores in its own footprint, so its supply test uses their
        # mean grain size, not the box mean.  (Mean, not minimum: the pores under one footprint
        # feed it in parallel.  Clusters that touch no pore never grow -- touches[] is False.)
        d_sub_cluster = np.full(n_c + 1, d_sub_mean)
        for c in range(1, n_c + 1):
            sel_c = lab[:, 0] == c
            if sel_c.any():
                d_sub_cluster[c] = float(d_sub_col[sel_c].mean())
        h_lim = np.full(n_c + 1, np.inf)
        zmax = np.zeros(n_c + 1)
        if n_c:
            from scipy import ndimage
            zmax[1:] = ndimage.maximum(zidx, lab, index=np.arange(1, n_c + 1))
        for c in range(1, n_c + 1):
            if not touches[c]:
                h_lim[c] = 0.0                                  # cut from the base: nothing grows
                continue
            if A_foot[c] <= 0:
                continue
            h_top = (zmax[c] + 1) * a
            hl_grow = first_failure_height(U_of(c, h_top), d_sub_cluster[c], 0.0, c_k, h_top, structure)
            if not np.isfinite(hl_grow):
                continue                                        # full growth demand is supplied: no limit
            lo, hi = 0.0, h_top                                 # Phi(lo) >= 0 (stationary demand), Phi(hi) < 0
            hl_lo = first_failure_height(U_of(c, lo), d_sub_cluster[c], 0.0, c_k, h_top, structure)
            if np.isfinite(hl_lo) and hl_lo <= lo:
                h_lim[c] = 0.0; continue
            for _ in range(25):
                mid = 0.5 * (lo + hi)
                hl = first_failure_height(U_of(c, mid), d_sub_cluster[c], 0.0, c_k, h_top, structure)
                if np.isfinite(hl) and hl < mid:
                    hi = mid
                else:
                    lo = mid
            h_lim[c] = 0.5 * (lo + hi)
        # frozen faces: above h_lim (or in unsupplied clusters).  They keep evaporating and their
        # precipitate fills the cell (plan/16 M2) instead of adding volume.
        frozen_faces = []
        for ax, s_, mask, nb, j in faces:
            frozen = mask & (zc >= h_lim[nb])
            Vn_face[(ax, s_)] = np.where(frozen, 0.0, Vn_face[(ax, s_)])
            frozen_faces.append((ax, s_, frozen, j))
        # ---- the salt an evaporating face leaves in the pores of the cell behind it
        # (corrected 2026-09-06, ideas/计算结果系统评估_2026-09-06.md section 3; same rule as the
        # 3-D driver).  A face evaporating j hands over BETA_SALT * j * A of salt; the growth law
        # turns the fraction F(Da*) of it into advance and only where the face may advance, so
        #     into the pores = (1 - F * [not frozen]) * BETA_SALT * j * A.
        # A frozen face therefore leaves ALL of it (plan/16 M2) and a growing face leaves the
        # (1 - F) kinetic share.  Before this date both used F and the remainder was dropped.
        dm_rate = np.zeros(g.shape)       # into a voxel's reaction band  [kg/s]
        dm_front = np.zeros(g.shape)      # into the bulk of the voxel the front is inside [kg/s]
        for ax, s_, frozen, j in frozen_faces:
            # Only faces that CANNOT advance reach the band; an advancing face leaves its
            # (1 - F) share in the material it lays down, through EPS_B (task B).
            rate_gas = np.where(frozen, np.maximum(j, 0.0) * BETA_SALT * g.face_area, 0.0)
            dm_rate += g.shift(rate_gas, ax, -s_, 0.0)
        # A face evaporating through the crust feeds the FRONT voxel -- the sealed voxel on
        # the other side of it -- until that voxel's bulk is down to eps_min.  Only a finished
        # voxel passes the salt on to the wet voxel behind, which then starts its own band.
        for ax, s_, m, j in crust_faces:
            rate_crust = np.where(m, np.maximum(j, 0.0) * BETA_SALT * g.face_area, 0.0)
            unfinished = m & (crust.eps_s > eps_min + crust.EPS_TOL)
            dm_front += np.where(unfinished, rate_crust, 0.0)
            dm_rate += g.shift(np.where(m & ~unfinished, rate_crust, 0.0), ax, -s_, 0.0)
        dm_rate = np.where(eng.solid_mask(), dm_rate, 0.0)
        # ---- salt the BARE SUBSTRATE releases.  fld.base_flux is non-zero only where the z = 0
        # cell is GAS or a conducting crust cell, so `solid0` splits it exactly: covered columns
        # precipitate into the voxel standing on them, bare columns have nowhere to put it.
        dm_base = np.zeros(g.shape)
        base_salt_rate = BETA_SALT * np.maximum(fld.base_flux, 0.0) * g.face_area   # [kg/s] per column
        solid0 = eng.solid_mask_full()[..., 0]      # deposit() checks the full-grid mask
        S_base_uncovered = float(base_salt_rate[~solid0].sum())
        S_base_released = float(base_salt_rate.sum())
        dm_base[..., 0] = np.where(solid0, base_salt_rate, 0.0)
        # ---- time step (plan/17 P1a item 4).  r scales the ACCURACY limiters (growth occupancy,
        # cap); the EVENT limiters (sealing, end of run) are not scaled: they place the step ON the
        # event, and scaling them by r < 1 only approaches the event geometrically (a Zeno
        # sequence of ever smaller steps, found 2026-09-05 in the verification cases).
        # ---- regularisation on the SALT, not on the volume (plan/17 P9, ported to 2-D
        # 2026-09-07).  Smoothing the volume moves material into columns whose birth porosity
        # differs, which does not conserve salt: the salt_structure balance caught it at 1.8e-6
        # the moment the pores stopped being identical.  What the liquid actually redistributes
        # laterally is the salt it carries, so the salt is smoothed and the volume follows from
        # the birth porosity where the salt LANDS:
        #     volume = salt / (rho_cr (1 - eps_birth))  and  1 - eps_birth = (1 - eps_s0) / F
        # hence volume = smoothed(salt * F) / (rho_cr (1 - eps_s0)).
        # Faces with Vn_face == 0 are excluded: a frozen face's salt goes to the band instead,
        # and counting it here would place it twice.
        def _salt_face(ax, s_, mask, weight=None):
            v = BETA_SALT * np.maximum(fld.face_flux(ax, s_), 0.0) * (Vn_face[(ax, s_)] > 0.0)
            return v if weight is None else v * weight
        salt_rate = eng.volume_rate(lambda ax, s_, mask: _salt_face(ax, s_, mask),
                                    smooth_length=ell)                      # [kg / (m s)]
        saltF_rate = eng.volume_rate(lambda ax, s_, mask: _salt_face(ax, s_, mask, F_grid),
                                     smooth_length=ell)
        rate = saltF_rate / (RHO_CR * (1 - eps_s0))
        # F at the cell is the salt-weighted mean of the F's that fed it; where nothing arrives,
        # the column's own F stands so that a cell born later still has a valid value.
        F_cell = np.divide(saltF_rate, salt_rate, out=np.array(F_grid, float), where=salt_rate > 0.0)
        crust.F_field = F_cell
        # The band path drains the bulk too -- deposit_skin removes the same kg from the pore
        # space it adds to the band -- so a cell fed by BOTH paths must be limited on their SUM,
        # not on each separately (found 2026-09-07: limiting them apart overfilled a voxel to
        # eps_s = -5.1e-7).  dt_to_seal(.., eps_min) places the step on the bulk sealing event;
        # once every fed cell is already at or below eps_min it goes inactive and
        # dt_to_seal(.., 0.0) takes over, keeping the pore space from going negative.
        # Only armed when base salt is actually being placed: with dm_base = 0 the band alone can
        # never overfill a voxel (its capacity costs the bulk only (eps_b - eps_min) delta_C / a),
        # so runs before 2026-09-07 keep the step sequence they had.
        dm_bulk = dm_rate + dm_base + dm_front
        dt = min(dt_factor * eng.dt_occupancy(rate, f_max=F_MAX),
                 crust.dt_to_close_skin(dm_rate / g.face_area, eps_min),
                 crust.dt_to_seal(dm_bulk, eps_min) if (dm_base.any() or dm_front.any()) else np.inf,
                 crust.dt_to_seal(dm_bulk, 0.0) if (dm_base.any() or dm_front.any()) else np.inf,
                 dt_factor * DT_CAP, t_end_h * 3600 - t)
        m_cell += salt_rate * dt
        rep = eng.advance(rate, dt, carry=m_cell)                            # 7 advance geometry
        if rep.flipped.any():
            # 8 new voxels get their properties from the flip mask, never from a porosity test
            crust.birth(rep.flipped, t + dt)
            # ... and their birth porosity from the salt the cell actually accumulated, which is
            # exact even when F changed between the steps that filled it.
            eps_from_salt = 1.0 - m_cell[rep.flipped] / (RHO_CR * g.cell_volume)
            crust.eps_s[rep.flipped] = eps_from_salt
            crust.eps_birth[rep.flipped] = eps_from_salt
            crust.update_properties(rep.flipped)
            m_cell[rep.flipped] = 0.0
        # salt_grown is read off the FACE RATE, not off the flipped voxels: salt that has
        # precipitated but not yet completed a voxel is still precipitated, and the precipitation
        # balance has to see it in the same step as the evaporation that produced it.
        # the salt the advancing faces hand over, read off the FACE rate: salt that has
        # precipitated but not yet completed a voxel is still precipitated, and the balance has to
        # see it in the same step as the evaporation that produced it.
        S_grow = BETA_SALT * g.face_area * sum(
            float((np.maximum(j, 0.0) * (Vn_face[(ax, s_)] > 0.0))[mask].sum())
            for ax, s_, mask, nb, j in faces)
        led.add("salt_grown", S_grow * dt)
        led.add("water_evaporated", fld.M_total * dt)
        led.add("salt_released_by_evaporation", BETA_SALT * W_salt * dt)
        led.add("volume_unplaced", rep.volume_banked)
        if dm_rate.any():                                # 8 deposition into the reaction band
            crust.deposit_skin(dm_rate * dt / g.face_area, eps_min)
            led.add("salt_front", float(dm_rate.sum()) * dt)
            seal_events += int(crust.seal_skin(eps_min, t + dt).sum())       # 9 state transitions
        led.add("salt_released_by_base", S_base_released * dt)
        led.add("salt_base_uncovered", S_base_uncovered * dt)
        if dm_front.any() or dm_base.any():              # 8b bulk pore filling, not the band:
            crust.deposit((dm_front + dm_base) * dt)     #    the front finishing its own voxel,
            led.add("salt_front", float((dm_front + dm_base).sum()) * dt)   # and the base salt
            seal_events += int(crust.seal(eps_min, t + dt).sum())
        t += dt
        wet = eng.salt_mask(); sealed = eng.state == SEALED; salt_all = wet | sealed
        # 11 re-solve the field that is written to the record, so that E(t) belongs to state(t)
        # and not to the state one step earlier (plan/17 P1a item 6).  The next iteration reuses
        # it through the same signature test, so this costs no extra solve.
        sig = state_signature(wet, sealed)
        if sig != last_sig:
            fld = solve_vapor(g, eng.state, props, solver="direct", top_external_length=L_ext,
                              porous=vapour_porous())
            last_sig = sig
        # blocked column: no gas cell of the column is connected to the far field (its wet faces, if
        # any, evaporate into sealed pockets and contribute nothing to M_total: no double counting).
        # Columns roofed by a neighbour's overhang count as blocked as well.
        # blocked column: the column's evaporating surfaces (wet faces, bare substrate) have no gas
        # path to the far field, or none are left.  Testing "no far-field gas anywhere in the
        # column" never fires once the sealed roof sits below the box top (air above it is
        # connected) -- corrected 2026-09-05 together with the 3-D driver.
        gas_now = eng.state == GAS
        evap_gas = np.zeros(g.shape, bool)
        for ax, s_ in g.face_directions():
            evap_gas |= gas_now & g.shift(wet, ax, s_, False)
        evap_gas[:, 0] |= gas_now[:, 0]                       # bare substrate evaporates
        blocked = ~(evap_gas & fld.active_gas).any(axis=1)
        f_sealed_top = float(blocked.mean())
        E_vap = fld.M_total / M_ref
        E_tot = E_vap                     # the leak is resolved; the column stays for old readers
        hs = mo.height_stats(salt_all, g)
        hist.append((t / 3600, E_vap, E_tot, float(salt_all[:, 0].mean()), hs["h_max"], f_sealed_top,
                     float(sealed.sum()) * g.cell_volume, int(sum(fr.sum() for _, _, fr, _ in frozen_faces))))
        if snap_hours and t / 3600 >= next_snap:
            snaps.append((t / 3600, eng.state.copy())); next_snap += snap_hours
        if hs["h_max"] >= nz * a - a:
            break
    hist = np.array(hist, float)
    if verbose:
        print(f"[E2 d_sub={d_sub*1e6:.0f} um] grid {g.shape} a={a*1e6:.1f} um z_box={nz*a*1e3:.1f} mm "
              f"L_ext={L_ext*1e3:.1f} mm j_plane={j_plane:.2e} F={F:.2f} | {len(hist)} steps {time.time()-tic:.0f} s | "
              f"t_end {hist[-1,0]:.0f} h theta {hist[-1,3]:.2f} h_max {hist[-1,4]*1e3:.2f} mm E_vap {hist[-1,1]:.2f} "
              f"E_tot {hist[-1,2]:.2f} sealed-top {hist[-1,5]:.2f} sealed cells {seal_events}")
    # what the structure actually holds, plus the solid still in flight (partly filled cells and
    # the engine's bank), which salt_grown already counted
    # m_cell IS the salt still in flight in the partly filled cells, and bank_carry is the part
    # that went to the engine's bank with its volume.  No porosity conversion is needed.
    led.add("salt_supplied", crust.salt_solid_mass() + float(m_cell.sum()) + eng.bank_carry)
    led.check("salt_structure")
    led.check("salt_precipitation")
    return dict(hist=hist, snaps=snaps, a=a, grid=g, state=eng.state.copy(), d_sub=d_sub,
                j_plane=j_plane, eps_s=crust.eps_s.copy(), m_precipitated=crust.m_precipitated.copy(),
                m_skin=crust.m_skin.copy(), skin_capacity=crust.skin_capacity(eps_min),
                eps_birth=EPS_B, delta_C=delta_C, d_sub_col=d_sub_col.copy(),
                F_col=F_col.copy(), band_col=band_col.copy(), eps_b_col=eps_b_col.copy(),
                ledger=led.report(), closure=closure.describe(), ell=ell,
                params=dict(tortuosity_exponent=tortuosity_exponent, cv=cv, pore_seed=pore_seed,
                            seed_fraction=seed_fraction, w_seed=w_seed,
                            salt_structure=structure.name,
                            d_sub=d_sub, n_seeds=n_seeds, cells_per_pitch=cells_per_pitch,
                            z_box_mm=z_box_mm, delta_eff_mm=delta_eff_mm, rh=rh, c_k=c_k,
                            eps_min=eps_min, t_end_h=t_end_h, dt_factor=dt_factor,
                            ell_over_dsub=ell_over_dsub,
                            U_factor=U_FACTOR, dt_cap=DT_CAP, f_max=F_MAX))

