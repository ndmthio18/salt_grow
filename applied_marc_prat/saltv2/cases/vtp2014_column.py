"""A resolved pore column under the bead-size series of Veran-Tissoires & Prat (2014, JFM 749).

Regime: evaporation-wicking -- a saturated bead column fed by a Dirichlet reservoir of nearly
saturated brine at its bottom, evaporating at its top into still air (their second series:
recess delta = 15 mm, RH 7 %, 22 C, four bead classes; section 6.1, figures 22 and 23).

The substrate is ONE vertical column of identical pores, as long as the bead column, built with
the geometry of the paper's OWN pore network (section 4.3 and appendix), and the voxel box is
one pore spacing wide with one surface pore in it.  Every surface pore is therefore the same
pore: the case can say whether and when the structure over a pore seals and what that does to
the evaporation rate (figure 23), and it CANNOT say anything about which pores nucleate, how far
apart the patches stand or how they compete -- coverage carries no patchy/crusty information
here.  User's decisions of 2026-09-20: fully resolved column, opening = the pore body's face.
"""
from pathlib import Path

import numpy as np

from ..materials import brine
from ..pnm import PoreLayer, boundary_layer_evaporation, write_column_dgf

# ---- the experiment ---------------------------------------------------------------------------
TUBE_LENGTH = 50e-3         # hollow cylinder holding the beads [m, M] (section 2)
RECESS = 15e-3              # bead surface below the rim, their delta [m, M]
L_COLUMN = TUBE_LENGTH - RECESS     # height of the bead column = the salt transport length
                            # [m, D].  Their table 2 agrees: Pe = J L / (rho A eps D*) = 1.48 at
                            # J = 1.80e-8 kg/s gives L = 35 mm.  It is NOT L_s = 22 mm, the height
                            # of the bead surface over the reservoir level (a hydrostatic head).
C0 = 0.25                   # reservoir and initial salt mass fraction [-, M]
RH = 0.07                   # enclosure, LiBr solution [-, M]
DELTA_EFF_MM = 18.0         # far-field distance that reproduces their table 1 evaporation rates
                            # [mm, B] (ideas/plan/14 section 12: mean of the four L_s = 22 mm rows
                            # 18.1 mm, +-4.3 %; absorbs free convection and a 20/22 C mismatch)
# Bead classes [m, M]: the middle of each range.  The third class is 100-120 um in the paper's
# text, 100-125 um in its tables 1-2 and the legend of figure 23, and 200-220 um on the row label
# of figure 22; the tables are followed.  "crusty"/"patchy" is the paper's own call (figure 22).
BEADS = {"5-50um": 27.5e-6, "50-62um": 56e-6, "100-125um": 112.5e-6, "1mm": 1062.5e-6}
OBSERVED = {"5-50um": "crusty", "50-62um": "crusty", "100-125um": "patchy", "1mm": "patchy"}

# ---- their pore network, per bead diameter d_b (section 4.3, appendix) [L] ----------------------
EPS_BED = 0.363                     # porosity of the packing (table 2 footnote)
SPACING_OVER_DB = 1.0               # cubic network, node distance = one bead diameter
THROAT_RADIUS_OVER_DB = 0.15        # mean of their truncated Gaussian, [0.09, 0.21] d_b
PORE_VOLUME_OVER_DB3 = EPS_BED      # mean element volume, [0.3, 0.426] d_b^3
# This project's pores are cubes that open into the gas domain with a whole face (user's choice
# over the paper's 0.3 d_b model openings), so the cube of that volume fixes the radius and what
# is left of the spacing is the throat [A]:
PORE_RADIUS_OVER_DB = 0.5 * PORE_VOLUME_OVER_DB3 ** (1.0 / 3.0)            # 0.3567
THROAT_LENGTH_OVER_DB = SPACING_OVER_DB - 2.0 * PORE_RADIUS_OVER_DB        # 0.2866
# DECLARED MISMATCH.  The paper calibrates its throat conductances on the continuum values
# (g = alpha r^4 / d_b with alpha = 1.29 for Kozeny-Carman, g_D = beta D_s r^2 / d_b with
# beta = 10.76 for eps^1.4 D_s, D_s = 1.3e-9).  Here the DuMuX rules of saltv2.pnm act on the
# geometry above (square duct, D A / L with D = 1.54e-9): 1.5 x their permeability and 1.54 x
# their eps D*, so the column Peclet number is 1.0 where theirs is 1.48.  Both deliver
# more than 96 % of the advected salt to the surface pore, so the supply is the same within 2 %.

# ---- closure and numerics: the production values (saltv2.cases.scan_crust.CRUST_COMMON) --------
C_K = 1e-5                  # [B]; the second set uses C_K_RECALIBRATED
C_K_RECALIBRATED = 3e-5     # [B] runs/E2_feedback/README.md section 12, on the conserving model
EPS_MIN = 0.02
ELL_OVER_DSUB = 0.25
VPP = 7                     # voxels per pore spacing: the opening 0.7134 d_b rasterises to 5 of 7
                            # cells = 0.7143 d_b (10 cells would give 8 = 0.8 d_b, 25 % too much
                            # wet area); VPP_FINE = 14 gives 10 of 14, the same opening
VPP_FINE = 14
Z_TOP_MM = {"5-50um": 1.0, "50-62um": 2.0, "100-125um": 3.0, "1mm": 8.0}
T_END_H = {"5-50um": 60.0, "50-62um": 150.0, "100-125um": 150.0, "1mm": 150.0}
RECORD_HOURS = 0.01         # a row every 36 s: the fine columns take millions of steps


def column_geometry(name):
    """The column of one bead class: what ``write_column_dgf`` takes, plus d_b itself."""
    d_b = BEADS[name]
    return dict(d_b=d_b, n_pores=int(round(L_COLUMN / (SPACING_OVER_DB * d_b))),
                spacing=SPACING_OVER_DB * d_b, pore_radius=PORE_RADIUS_OVER_DB * d_b,
                pore_volume=PORE_VOLUME_OVER_DB3 * d_b ** 3,
                throat_radius=THROAT_RADIUS_OVER_DB * d_b,
                throat_length=THROAT_LENGTH_OVER_DB * d_b)


def write_dgf(name, directory):
    """Write the column of bead class ``name`` into ``directory``; returns the path."""
    geo = column_geometry(name); geo.pop("d_b")
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    return write_column_dgf(directory / f"vtp2014_column_{name}.dgf", **geo)


def case(name, dgf_path, c_k=C_K, vpp=VPP):
    """Keyword arguments of ``saltv2.coupling.driver3d.run`` for bead class ``name``.  The grain
    size of the salt closures IS the bead diameter, and neighbouring surface pores stand one bead
    diameter apart."""
    d_b = BEADS[name]
    return dict(substrate_dgf=str(dgf_path), substrate_L_box=SPACING_OVER_DB * d_b, d_sub=d_b,
                C0=C0, rh=RH, delta_eff_mm=DELTA_EFF_MM, z_top_mm=Z_TOP_MM[name], vpp=vpp,
                c_k=c_k, eps_min=EPS_MIN, ell_over_dsub=ELL_OVER_DSUB,
                t_end=T_END_H[name] * 3600.0, record_hours=RECORD_HOURS,
                tag=f"vtp2014_{name}")


def bare_column_peclet(name, dgf_path):
    """Peclet number of the bare column [-], n ln(1 + Q/g) for n identical throats of diffusive
    conductance g, with Q = water / rho the liquid that replaces what one box floor evaporates as
    a fully wet plane (a porous surface much finer than the far-field distance evaporates like
    one; the same Q = J A / rho as saltv2.cases.pnm_column)."""
    layer = PoreLayer.from_dgf(dgf_path)
    d_lat, d_res = layer.diffusive()
    g = float(d_lat[0])
    if not (np.allclose(d_lat, g) and np.isclose(d_res[layer.bottom][0], g)):
        raise ValueError("bare_column_peclet assumes identical throats")
    box = (SPACING_OVER_DB * BEADS[name]) ** 2
    water = float(boundary_layer_evaporation(C0, box, DELTA_EFF_MM * 1e-3, RH))     # [kg/s]
    Q = water / float(brine.density(C0))
    return layer.N * float(np.log1p(Q / g))
