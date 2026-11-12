"""The self-consistent baseline: evaporation -> wicking -> transport -> crystallisation at C_sat ->
porous growth by the Prat law, with ONE microstructure behind every property of the salt and
nothing that was put in to make a crust.

Regime: evaporation-wicking -- a saturated bead column fed by a Dirichlet reservoir at its
bottom, one surface pore opening into a periodic box one pore spacing wide.  The columns are the
resolved ones of ``vtp2014_column`` (the four bead classes of Veran-Tissoires & Prat 2014, their
own network geometry) and the DuMuX reference column.  Purpose (user's decision of 2026-09-20):
NOT agreement with the experiments but a coupled model that is consistent with itself; the six
criteria are evaluated by ``runs/E1_pre/baseline.py``.

What differs from ``vtp2014_column`` (the archived crust chain), each stated here and not
inherited:
  salt_structure  the caller's, SINGLE_SCALE or HIERARCHICAL -- never the legacy pairing
  c_k = 1         the cubic law as written; no conductance reduction [L].  The keyword stays so
                  that a later scan can ask what k_s a viscous-limitation hypothesis would need
  T               the temperature of the experiment, carried by every closure that has one
  delta_eff       re-derived from the same table-1 rates with the 22 C gas constants
Expected, and the point of the exercise: no bead size forms a crust -- with the cubic law as
written the viscous drop inside the salt is negligible, and what is left of the supply need
(gravity, the column's own viscous drop) stays below the capillary threshold over any height
the box allows -- so the structure keeps growing.
That agrees with the 1 mm class and disagrees with the fine classes; it is a parameter-free
statement about "Prat law + capillary supply through one microstructure".
"""
import numpy as np

from ..materials import air, brine
from ..materials.brine import C_SAT
from . import dumux1d, vtp2014_column as vtp

# ---- the experiment (section 6.1 of the paper; geometry: vtp2014_column) -----------------------
T = 295.15                  # [K, M] 22 C
RH = vtp.RH                 # 0.07 [M]
C0 = vtp.C0                 # 0.25 [M]
BEADS = vtp.BEADS
# their table 1, second series (recess 15 mm, L_s = 22 mm): measured evaporation rate of the
# tube before efflorescence [kg/s, M], as transcribed and corrected in ideas/plan/14 section 12
TABLE1_RATE = {"1mm": 1.80e-8, "100-125um": 1.74e-8, "50-62um": 1.75e-8, "5-50um": 1.9e-8}
TUBE_AREA = np.pi * (19e-3) ** 2        # [m^2, M] inner cross-section of the tube


def calibrated_delta_eff(T_K):
    """Far-field distance [m] at which pure diffusion from a saturated-brine plane gives the mean
    of the four measured rates: delta = rho_g D_v W_sat (a_w(C_sat) - RH) / j, the procedure of
    ideas/plan/14 section 12 with the gas constants of temperature ``T_K`` [B].  It is NOT a
    geometric length: it absorbs the free convection in the enclosure (factor 1.49 on the
    paper's own pure-diffusion rate)."""
    gas = air.moist_air(T_K, RH)
    numerator = gas.rho_D * gas.W_sat * (float(brine.a_w(C_SAT)) - RH)
    return float(np.mean([numerator / (rate / TUBE_AREA) for rate in TABLE1_RATE.values()]))


DELTA_EFF_MM = round(calibrated_delta_eff(T) * 1e3, 1)     # 20.6 mm [B]; 18.1 with the 20 C
                            # constants: the 14 % of rho_g D_v W_sat between 20 and 22 C that
                            # the old value had absorbed

# ---- closure and numerics ---------------------------------------------------------------------
C_K = 1.0                   # [L] the cubic law of Wijnhorst 2025 as written
EPS_MIN = 0.02              # sealing porosity [A]; no voxel is expected to reach it
ELL_OVER_DSUB = 0.25        # regularisation length [N]
VPP = 7                     # 5 of 7 cells = the 0.7134 d_b opening (see vtp2014_column)
VPP_FINE = 14
Z_TOP_MM = dict(vtp.Z_TOP_MM)
T_END_H = dict(vtp.T_END_H)
RECORD_HOURS = 0.01


def case(name, dgf_path, salt_structure, c_k=C_K, vpp=VPP):
    """Keyword arguments of ``driver3d.run`` for bead class ``name`` (column file:
    ``vtp2014_column.write_dgf``).  ``salt_structure`` has no default on purpose."""
    d_b = BEADS[name]
    return dict(substrate_dgf=str(dgf_path), substrate_L_box=vtp.SPACING_OVER_DB * d_b, d_sub=d_b,
                C0=C0, rh=RH, T=T, delta_eff_mm=DELTA_EFF_MM, z_top_mm=Z_TOP_MM[name], vpp=vpp,
                salt_structure=salt_structure, c_k=c_k, eps_min=EPS_MIN,
                ell_over_dsub=ELL_OVER_DSUB, t_end=T_END_H[name] * 3600.0,
                record_hours=RECORD_HOURS, tag=f"baseline_{name}_{salt_structure.name}")


def dumux_column(salt_structure, c_k=C_K):
    """The DuMuX reference column with its own far field (delta = 1.6 mm, RH 0.2) and its own
    temperature, 293.15 K -- the state its gas constants were always at, so ``T`` moves the
    liquid side and k_r by one degree and the gas side by ~0.1 % (rho_g is now averaged with
    a salt surface instead of a pure-water one)."""
    kw = dict(dumux1d.RUN_3D_GROWTH)
    kw.update(salt_structure=salt_structure, c_k=c_k, T=dumux1d.T_DUMUX,
              tag=f"baseline_dumux1d_{salt_structure.name}")
    return kw
