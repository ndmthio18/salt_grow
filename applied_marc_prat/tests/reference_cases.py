"""The small production-path cases whose outputs are pinned bit for bit (cleanup of 2026-09-14).

Regime: evaporation-wicking.  Each case is a driver call that finishes in seconds; its full
output (history, ledger, per-voxel state and salt masses) was captured from the checkpoint tree
95c1acb ("checkpoint before cleanup") into ``tests/data/ref_*.npz`` by ``make_references.py``,
and ``test_regression.py`` re-runs the case and demands exact equality.  The point is that the
cleanup removes code, never behaviour: a single ulp of difference is a defect, not noise.

Intended behaviour changes since, each re-captured on purpose and recorded here:
  * 2026-09-15  ref_3d_small: surface-pore openings are the square top face of the QSPNM cube
    instead of the inscribed disk (saltv2/interface.py, runs/E1_pre/README.md section 12).  At
    vpp 3 the square fills the whole pitch; the 2-D references are not affected.
  * 2026-09-15  all three re-captured on frodo (platform move, not a behaviour change): the 2-D
    ones from a0db925, ref_3d_small from f0ff577.  The local-machine captures differ from these
    at the ulp level only.
  * 2026-09-16  ref_3d_small: a cluster's draw is split between its feeding pores by opening
    area instead of equally, and s_foot is area-weighted (dfb5db9, README E1 section 12).  Same
    steps, states and sealed count; ledger ~1e-6 relative.
  * 2026-09-17  ref_3d_small re-captured at 26f9331, covering two intended changes that the
    reference had not followed: the localized nucleation seed with source-resolved evaporation
    and a resolved reaction band (06e6a99), and the DuMuX brine properties with square throats
    (cf5d793: D_S 1.5e-9 -> 1.54e-9, throat area 4 r^2 instead of pi r^2, density rho(w)).
    1323 -> 1100 steps, 648 -> 516 sealed voxels, salt_grown and salt_released_by_evaporation
    down 16 / 17 %, water_evaporated_base up 3.1x (the localized seed leaves more bare opening),
    and the new channel salt_pnm_precipitate appears at 0 (the option is off).
  * 2026-09-17  all three re-captured at 26f9331 for the ledger's channel LIST, not for any
    value: water_evaporated_on_salt joins it (tests/test_kinetic_share.py), and the 2-D pair
    still carried four channels the 2026-09-14 cleanup had removed.  Every array and every
    shared ledger value of ref_2d_small / ref_2d_cv is bit for bit what it was; only the names
    moved.
  * 2026-09-17  ref_3d_small: the network balances liquid MASS instead of liquid VOLUME
    (pre_experiment.solve_pressure, runs/E1_pre/README.md section 13.7).  The reservoir's brine
    is denser than the pores' (1232.65 vs 1200.97 at C0 = 0.30 and C_sat), so the advective salt
    import falls 2.6 % and salt_pnm_cap -- a difference of two larger numbers -- falls 9.2 %;
    salt_grown is unchanged to 1e-6, 1100 -> 1097 steps, the same 516 sealed voxels.  The 2-D
    references are bit for bit identical (feedback2d has no pore network).
  * 2026-09-17  all three again, for the water activity: a_w is the DuMuX reference case's
    simplified Pitzer model instead of a straight line to 0.75 (saltv2/brine.water_activity,
    runs/E1_pre/README.md section 13.5).  a_w(C_sat) 0.7500 -> 0.7530, so the 2-D cases, whose
    salt surfaces all sit at C_sat, simply evaporate 0.15 / 0.55 % more.  The 3-D case also has
    BARE openings, where the curve differs in shape rather than at the anchor, and its feedback
    absorbs the rest: same 1100 steps and 516 sealed voxels, salt_grown -0.10 %,
    water_evaporated_base +0.12 %, water_evaporated -0.09 %.  The normalised E_total is
    unchanged at t = 0, because the driver divides by a plane flux built from the same a_w.
  * 2026-09-19  ref_3d_small: the pore's liquid-mass sink no longer carries rho_pore / RHO_L
    (saltv2/coupling/driver3d.py, runs/E1_pre/README.md section 15.5; found by the coupling check
    against DuMuX, where it was 4.6 %).  At this case's C0 = 0.30 the factor was +2.7 % at the
    start and ~0 at C_sat: 1097 -> 1093 steps, the same final state and 516 sealed voxels,
    salt_grown +0.003 %, water_evaporated +0.22 %, salt_pnm_cap -0.16 %.  The 2-D references are
    bit for bit identical (no pore network).
  * 2026-09-19  ref_3d_small again: brine leaves a feeding pore at the pore's concentration and
    carries the salt of every wet face, growing or not (563004c, README section 15.7).  C0 = 0.30
    keeps the covered pores at the C_sat cap, so beta itself does not move: same final state and
    516 sealed voxels, 1093 -> 1109 steps, salt_grown +0.20 %, water_evaporated +1.1 %; the
    larger liquid draw of the non-growing faces imports more reservoir salt, salt_pnm_cap +19 %.
    The 2-D references are bit for bit identical.
  * 2026-09-19  ref_3d_small once more: a pore under salt is held at C_sat and what crystallises
    there is the structure's only salt, at most beta(C_sat) per kg of water (4017046, README
    section 15.8).  C0 = 0.30 delivers more than that, so the faces still get eq (29): 1109 ->
    1087 steps, 516 -> 507 sealed voxels, salt_grown -0.49 %.  The ledger gains salt_crystallised
    and salt_interface_buffer; salt_pnm_cap is 0 (its 27 % of the reservoir salt now sit, named,
    in the buffer: 25 % of what crystallised).  The 2-D references are bit for bit identical.
  * 2026-09-21  all three re-captured on frodo from a clean clone of aaa7058 (branch
    baseline-selfconsistent, stage A; 1027 / 624 / 1689 rows), ONCE for the five intended
    changes below.  Each was measured on its own, on one machine, against the step before it
    (runs/E1_pre/README.md section 18); the ledger of the 3-D case also gains the liquid and
    pore-salt channels (liquid_in / _out / _stored / _nucleation_jump, salt_reservoir_in,
    salt_pore_inventory):
    - A2, the supply gate counts the salt of the faces that do not advance (one brine demand
      serves the gate and the network; user's decision, closes the open item of 2026-09-19):
      ref_3d_small 1087 -> 1023 steps, 507 -> 448 sealed voxels, h_max 0.109 -> 0.095 mm (one
      voxel layer), salt_grown -13.0 %, water_evaporated -12.0 %, first seal unchanged at
      6.75 h.  ref_2d_small salt_grown -2.7 %, 205 -> 229 sealed; ref_2d_cv +2.0 %, 761 -> 786.
      The DuMuX-column record (no frozen face) is bit for bit identical.
    - A3, one salt diffusivity: the Prat closures (Da*, reaction band) take the pore network's
      D_NACL = 1.54e-9 instead of prat2026's own default 1.3e-9.  Da* -15.6 % at every grain
      size, F 0.759 -> 0.744 at 34 um, so the material is laid down ~2 % denser:
      ref_3d_small salt_grown +2.1 %, same 448 sealed voxels and h_max, first seal 6.75 -> 6.82 h;
      ref_2d_small +2.1 % (229 -> 221 sealed), ref_2d_cv +2.3 % (786 -> 818).
    - A4, the network conserves liquid and salt in a transient (pnm.network_step: pore storage
      V d rho in the liquid balance, V d(rho w) in the salt balance): ref_3d_small salt_grown
      -0.004 %, same 1027 steps and 448 sealed voxels; the ledger gains liquid_in, liquid_out,
      liquid_stored and liquid_nucleation_jump.  The 2-D references (no pore network) are bit
      for bit identical.
    - A5, run(T=) declares the case temperature (gas side from the DuMuX expressions,
      materials.air.moist_air; brine density, k_r and Knudsen at the same T).  The references
      run with T = None, the legacy constants, and move only through crust.M_W = 0.018 ->
      materials.constants.M_H2O (D_K -0.04 %, the sealed crust's conductance): ref_3d_small
      and ref_2d_small < 0.001 % in every channel, same steps and states; ref_2d_cv 1676 -> 1689
      steps, 818 -> 831 sealed, salt_grown +0.03 %.

The kwargs below are written in the OLD signature (with the legacy flags at their production
values) so that the reference can be regenerated from the checkpoint.  ``call`` strips the
keywords the current driver no longer accepts, by name and only those, so a typo cannot pass
silently.
"""
from __future__ import annotations

import inspect
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DATA = Path(__file__).resolve().parent / "data"

# keywords that the cleanup removes from the drivers; passing them at their production value
# to the checkpoint tree reproduces exactly what the cleaned tree does without them
LEGACY_3D = dict(crust_vapour=True, phi_vent=0.0, dry=False, seal=True, base_salt=True,
                 conserve_salt=True, adaptive_box=True, allow_rewet=False)
LEGACY_2D = dict(crust_vapour=True, phi_vent=0.0, seal=True, base_salt=True)

CASES = {
    # the test_crust_vapour_3d fixture box: 2 x 2 pores, vpp 3, 0.5 mm, 24 h
    "ref_3d_small": ("3d", dict(lam=32, C0=0.30, Lt_res=22e-3, n_side=2, vpp=3, delta_eff_mm=18.0,
                                rh=0.07, z_top_mm=0.5, tag="test", verbose=False, ell_over_dsub=1.0,
                                eps_min=0.05, c_k=1e-6, t_end=24 * 3600, **LEGACY_3D)),
    # the 2-D grain-size box of runs/E3_grainsize at 34 um, 22 h
    "ref_2d_small": ("2d", dict(d_sub=34e-6, n_seeds=4, cells_per_pitch=6, z_box_mm=0.6, rh=0.07,
                                c_k=1e-5, eps_min=0.02, t_end_h=22.0, verbose=False,
                                ell_over_dsub=0.25, snap_hours=0.0, **LEGACY_2D)),
    # the same box with a lognormal pore-size spread (per-column d_sub / F / band path)
    "ref_2d_cv": ("2d", dict(d_sub=34e-6, n_seeds=8, cells_per_pitch=6, z_box_mm=0.6, rh=0.07,
                             c_k=1e-5, eps_min=0.02, t_end_h=24.0, verbose=False,
                             ell_over_dsub=0.25, snap_hours=0.0, cv=0.4, pore_seed=3,
                             **LEGACY_2D)),
}


def driver(kind):
    if kind == "3d":
        from saltv2.coupling import driver3d as v
        return v.run
    from saltv2.coupling import driver2d as fb
    return fb.run


def call(kind, kwargs):
    """Run the case on whichever driver signature is current."""
    fn = driver(kind)
    accepted = inspect.signature(fn).parameters
    legacy = LEGACY_3D if kind == "3d" else LEGACY_2D
    kw = {k: val for k, val in kwargs.items() if k in accepted or k not in legacy}
    return fn(**kw)


def capture(kind, res):
    """Everything the regression compares, as plain arrays."""
    led = res["ledger"]
    out = dict(hist=np.asarray(res["hist"], float),
               ledger_names=np.array(sorted(led)), ledger_values=np.array([led[k] for k in sorted(led)]),
               eps_s=np.asarray(res["eps_s"]) if kind == "2d" else res["crust"].eps_s,
               m_precipitated=np.asarray(res["m_precipitated"]) if kind == "2d" else res["crust"].m_precipitated,
               m_skin=np.asarray(res["m_skin"]) if kind == "2d" else res["crust"].m_skin)
    if kind == "3d":
        out.update(hist_columns=np.array(list(res["hist_columns"])), state=res["eng"].state,
                   t_birth=res["crust"].t_birth, t_sealed=res["crust"].t_sealed,
                   t_nuc=np.asarray(res["t_nuc"]), n_steps=np.array(int(res["n_steps"])))
    else:
        out.update(state=np.asarray(res["state"]), d_sub_col=np.asarray(res["d_sub_col"]),
                   n_steps=np.array(len(res["hist"])))
    return out


def meta():
    import scipy
    try:
        git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                             text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git = "unknown"
    return dict(numpy_version=np.__version__, scipy_version=scipy.__version__,
                platform=platform.platform(), python=platform.python_version(), git=git)
