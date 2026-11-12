"""The DuMuX pnm2pncmin2 1-D reference case (the network the PNM was validated against).

Regime: evaporation-wicking -- a saturated column fed by a Dirichlet reservoir at the bottom,
water leaving the top pore through a stagnant boundary layer.  The DuMuX inputs themselves are
runs/E1_pre/dumux_reference/params_1d_{dilute,saturated}.input; the values below are those inputs
in this project's variables.
"""
from ..materials import brine
from ..materials.air import W_SAT
from ..pnm.network import DGF_DUMUX_1D as DGF      # 6 cubic pores, the bottom one the reservoir
from ..pnm.network import PoreLayer

DELTA_BL = 1.6e-3           # stagnant boundary layer over the pore opening [m, L]: QSPNM's
                            # delta_top and DuMuX's FreeFlow.BoundaryLayerThickness, the same 1.6 mm
X_REF = 0.00293             # far-field vapour mass fraction [-, L] (QSPNM; DuMuX's
                            # ReferenceWaterVaporMoleFraction = 4.7e-3 in mass form)
RH_ENV = X_REF / W_SAT      # the same far field as a relative humidity, 0.19959
T_DUMUX = 293.15            # SpatialParams.Temperature the reference runs were made with [K]

X_NACL_DILUTE = 0.05        # InitialSaltMoleFraction of the pure-transport run [-]
X_NACL_SATURATED = 0.095    # ... of the run that reaches saturation [-]
C0_DILUTE = float(brine.mass_fraction(X_NACL_DILUTE))          # 0.14583
C0_SATURATED = float(brine.mass_fraction(X_NACL_SATURATED))    # 0.25402

_column = PoreLayer.from_dgf(DGF)
OPENING = float(2.0 * _column.r[_column.surface][0])      # the surface pore's cube face, 0.4 mm
PORE_SPACING = float(_column.dz)                          # 0.6 mm
# Grain size the salt closures are written in [m, A].  The file has no grains; its cubic pores
# (body 200 um, throat 110 um, spacing 600 um) are within 11 % of the voids of a simple cubic
# packing of spheres whose diameter is the pore spacing (body 0.366 d = 220 um, throat
# 0.207 d = 124 um), so the spacing is the grain size this network is consistent with.  Until
# 2026-09-20 the driver's default 4 r = 800 um was used (user's decision to change it).
D_SUB = PORE_SPACING

# The 3-D closed loop (saltv2.coupling.driver3d.run) on that network, at the reference case's own
# far field: the DuMuX column is the substrate, and the 1.6 mm boundary layer is the whole gas
# domain (a 1.0 mm voxel box plus a 0.6 mm external segment).  The closure values are those of the
# production cases (saltv2.cases.scan_crust).  There is NO DuMuX counterpart to compare this run
# with -- their case has no external salt -- so it is a case of this model, not a validation.
RUN_3D = dict(substrate_dgf=DGF, d_sub=D_SUB, C0=C0_SATURATED, rh=RH_ENV, delta_eff_mm=DELTA_BL * 1e3,
              z_top_mm=1.0, vpp=5, c_k=1e-5, eps_min=0.02, ell_over_dsub=0.25,
              t_end=12 * 3600, tag="dumux1d")

# ---- the column coupled to the external salt structure (2026-09-18) -------------------------

# Coupling check.  With the box as wide as the opening the gas problem is one-dimensional and the
# voxel vapour solve IS the DuMuX top boundary condition (rho_g D_v A (X_pore - X_ref) / delta),
# so the whole chain vapour solve -> pore sink -> pressure -> salt transport -> a_w can be held
# against their dilute run, which never nucleates (top pore 0.1997 < C_sat = 0.264).
# The case states the temperature their run was made at (since 2026-09-21): the brine density of
# the network is then taken at 293.15 K as theirs is, not at the project's 21 C.  Measured on the
# deviation below: 3.02e-3 without it, 2.72e-3 with it (rho_g also moves by +0.1 %, because
# moist_air averages it with a salt surface; the two contributions were not separated).
RUN_3D_CHECK = dict(substrate_dgf=DGF, substrate_L_box=OPENING, d_sub=D_SUB, C0=C0_DILUTE, rh=RH_ENV,
                    T=T_DUMUX, delta_eff_mm=DELTA_BL * 1e3, z_top_mm=1.0, vpp=4, c_k=1e-5,
                    eps_min=0.02, ell_over_dsub=0.25, t_end=40000.0, tag="dumux1d_check")
# their dilute run at 20 C, pores top -> bottom at t = 40000 s (runs/E1_pre/README.md section 13)
DUMUX_DILUTE_STEADY = (0.19969, 0.18770, 0.17635, 0.16560, 0.15544)

# The column under the far field of the one-pore runs of 2026-09-16 (runs/E1_pre/
# single_pore_growth_figures.py: RH 7 %, delta_eff 18 mm, 150 h), with the reservoir at the DuMuX
# saturated case's C0 = 0.254 instead of their 0.34 -- that value sat above saturation to skip
# the wait, whereas here nucleation has to come out of the column's own enrichment.
# Neighbouring surface pores are taken as far apart as the column's pores are [A, isotropic
# packing]: the opening covers 44 % of the box floor.  vpp = 12 puts 8 x 8 voxels of 50 um on it.
# In-pore precipitation stays off, as on the production path.
# RESULT.  The column Peclet number is only 0.064 (5.8e-12 kg/s through five throats of 200 um x
# (220 um)^2), so the surface pore would settle at C0 exp(0.064) = 0.2707.  With the nucleation
# threshold at 1.2 C_sat = 0.3168 (until 2026-09-19) that never nucleated; with no supersaturation
# it reaches C_sat after 0.47 h and then grows slowly for the whole 150 h without sealing, handing
# over 0.15 kg of salt per kg of water -- 42 % of eq (29), because Pe sits just above the
# ln(C_sat / C0) = 0.039 below which a crystal on this column cannot be fed at all.
RUN_3D_FARFIELD = dict(substrate_dgf=DGF, substrate_L_box=PORE_SPACING, d_sub=D_SUB,
                       C0=C0_SATURATED, rh=0.07, delta_eff_mm=18.0, z_top_mm=3.5, vpp=12,
                       c_k=1e-5, eps_min=0.02,
                       ell_over_dsub=0.25, t_end=150 * 3600, tag="dumux1d_farfield")

# Growth on the column: the same box and closure under the reference case's OWN far field, whose
# 1.6 mm boundary layer evaporates eight times faster (column Pe 0.56): C_sat after 0.02 h, growth
# to the box top region within 5 h, sealing.  The boundary layer is the whole gas domain, which
# caps the box at 1.5 mm.  Once covered, the surface pore holds C_sat and the column below falls
# linearly to C0.  Like DuMuX the pore crystallises AT C_sat; unlike DuMuX the crystal is the
# external structure and not a solid fraction inside the pore, so there is still no counterpart.
RUN_3D_GROWTH = dict(RUN_3D_FARFIELD, rh=RH_ENV, delta_eff_mm=DELTA_BL * 1e3, z_top_mm=1.5,
                     t_end=48 * 3600, tag="dumux1d_growth")
# ... and with the crust sealing two decades earlier, as in the second one-pore run (the shell)
RUN_3D_GROWTH_SHELL = dict(RUN_3D_GROWTH, c_k=1e-6, tag="dumux1d_growth_shell")
