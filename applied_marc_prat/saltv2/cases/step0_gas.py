"""Step 0 of the surface-first verification ladder: gas-phase evaporation from ONE liquid-filled
square surface pore -- nothing else.

Physical model.  A square opening of side PORE_SIDE, filled with NaCl brine at C_BRINE, sits in an
impermeable plane; above it a 3-D gas domain, periodic laterally with period L_BOX, in which water
vapour diffuses at steady state,  div(rho_g D_v grad X) = 0,  j = -rho_g D_v grad X.  Boundaries:
X = a_w(C_BRINE) W_sat(T) on the opening, zero normal flux on the solid, X = RH W_sat(T) at the far
field, which lies DELTA_EFF above the plane.  No salt, no nucleation, no liquid flow, no
temperature gradient, no convection.

Numerical model.  ``growth.vapor.solve_vapor`` (uniform voxels, finite volumes, direct solve) on a
box of height H_BOX; the gas between the box top and the far field, DELTA_EFF - H_BOX, is the
laterally uniform 1-D segment in series (the Robin top of ``fv.assemble_laplace``), exact once the
field at the box top is laterally uniform.  The lateral Fourier modes of the field decay as
exp(-2 pi z / L_BOX), so the box height is a NUMERICAL choice; the physical transfer distance is
DELTA_EFF in every configuration.  The central question of Step 0 is whether the result depends on
that choice (configurations A: 2 mm, B: 4 mm).

Pore-opening rasterisation.  The opening is the production footprint of
``coupling.interface.PoreOpenings`` (cells whose centre lies within half a side of the pore
centre).  With A_VOXEL = 50 um and a 540 um opening centred in the 700 um box that is 10 x 10
cells = 500 um: the voxelised area is 14.3 % smaller than the physical one.  This is reported,
not corrected (``rasterisation``); the half-voxel-shifted pore (11 x 11 cells, +3.8 %) is the
supplementary diagnostic of how much the evaporation rate depends on that rule.

Provenance tags: [M] measured / stated by the experiment, [L] literature or lattice geometry,
[B] calibrated, [N] numerical choice.

    from saltv2.cases.step0_gas import solve_plane, solve_pore, rasterisation
"""
from __future__ import annotations

import numpy as np

from ..coupling.interface import PoreOpenings
from ..growth import fv
from ..growth.grid import GAS, Grid
from ..growth.vapor import VaporField, VaporProperties, solve_vapor
from ..materials.air import moist_air
from ..materials.brine import C_SAT, a_w

# ---- case values ---------------------------------------------------------------------------
T = 295.15                  # [K, M] 22 C
RH = 0.07                   # far-field relative humidity [-, M]
C_BRINE = C_SAT             # NaCl mass fraction in the opening, 0.264 [-, M]
PORE_SIDE = 540e-6          # side of the square opening [m, L]: 2 x the 270 um QSPNM pore radius
L_BOX = 700e-6              # lateral period [m, L]: the QSPNM pitch
DELTA_EFF = 20.6e-3         # effective vapour-transfer distance plane -> far field [m, B]
A_VOXEL = 50e-6             # gas voxel size [m, N]
H_BOX = {"A": 2.0e-3, "B": 4.0e-3}      # numerical gas-domain heights [m, N]

assert C_BRINE == 0.264


def gas_properties():
    """Gas-side numbers at (T, RH): (MoistAir, VaporProperties)."""
    air = moist_air(T, RH, a_w_surface=float(a_w(C_BRINE)))
    props = VaporProperties(rho_D=air.rho_D, X_wet=float(a_w(C_BRINE)) * air.W_sat,
                            X_top=RH * air.W_sat)
    return air, props


def build(h_box, a=A_VOXEL, shift=0.0):
    """Grid of height ``h_box`` [m], the pore opening (centre shifted by ``shift`` [m] in both
    lateral directions, zero for the production placement) and the external segment length."""
    n_lat = L_BOX / a
    n_z = h_box / a
    if abs(n_lat - round(n_lat)) > 1e-9 or abs(n_z - round(n_z)) > 1e-9:
        raise ValueError("L_BOX and h_box must be integer multiples of the voxel size")
    grid = Grid((int(round(n_lat)), int(round(n_lat)), int(round(n_z))), a)
    openings = PoreOpenings(grid, [0.5 * L_BOX + shift], [0.5 * L_BOX + shift], [0.5 * PORE_SIDE], L_BOX)
    L_ext = DELTA_EFF - h_box
    if L_ext < 0.0:
        raise ValueError("the box is taller than the effective transfer distance")
    return grid, openings, L_ext


def rasterisation(h_box=H_BOX["A"], a=A_VOXEL, shift=0.0):
    """V0.4: physical vs voxelised opening area [m^2] and the relative error."""
    grid, openings, _ = build(h_box, a, shift)
    A_phys = PORE_SIDE ** 2
    A_raster = float(openings.voxel_area()[0])
    n_cells = int(openings.foot[0].sum())
    return dict(A_phys=A_phys, A_raster=A_raster, n_cells=n_cells,
                side_raster=float(np.sqrt(n_cells)) * a,
                area_rel_err=(A_raster - A_phys) / A_phys)


def discrete_residual(grid, props, fld: VaporField, bottom, L_ext):
    """||A x - b|| / ||b|| of the assembled system at the returned field: the solved field
    satisfies the discrete conservation equation of every gas cell (independent of the
    flux bookkeeping of ``solve_vapor``)."""
    state = np.full(grid.shape, GAS, int)
    rho_D = fv.as_field(props.rho_D, grid.shape, "rho_D")
    active = fv.top_connected(grid, state == GAS)
    mask, values = bottom
    system = fv.assemble_laplace(grid, rho_D, active, np.zeros(grid.shape, bool), props.X_wet,
                                 props.X_wet, props.X_top, bottom_mask=mask, bottom_values=values,
                                 top_external_length=L_ext)
    x = fld.X[active]
    return float(np.linalg.norm(system.matrix @ x - system.rhs) / np.linalg.norm(system.rhs))


def _outputs(grid, props, fld: VaporField, bottom, L_ext, h_box, A_wet):
    A_box = fld.projected_area
    E = fld.M_base
    dX = props.X_wet - props.X_top
    out = dict(
        h_box=h_box, L_ext=L_ext, delta_eff=DELTA_EFF, a_voxel=grid.a, shape=list(grid.shape),
        rho_D=float(props.rho_D), X_wet=props.X_wet, X_top=props.X_top,
        E=E,                                    # total evaporation rate [kg/s]
        E_top=fld.M_top,                        # vapour leaving through the far-field boundary [kg/s]
        A_box=A_box,                            # projected periodic area [m^2]
        A_wet=A_wet,                            # rasterised wet area [m^2]
        j_projected=E / A_box,                  # [kg m^-2 s^-1] per projected area
        j_wet=E / A_wet,                        # [kg m^-2 s^-1] per wet (rasterised) area
        imbalance=fld.imbalance,                # (E - E_top) / E_top
        discrete_residual=discrete_residual(grid, props, fld, bottom, L_ext),
        # series-resistance reading: E = rho_D dX A_box / (delta_eff + L_spread)
        L_spread=float(props.rho_D * dX * A_box / E - DELTA_EFF),
    )
    return out


def solve_plane(h_box):
    """Reference 1: the whole bottom plane wet.  Adds the analytical 1-D flux and its error."""
    grid, _, L_ext = build(h_box)
    _, props = gas_properties()
    state = np.full(grid.shape, GAS, int)
    mask = np.ones(grid.shape[:2], bool)
    bottom = (mask, np.full(grid.shape[:2], props.X_wet))
    fld = solve_vapor(grid, state, props, bottom=bottom, top_external_length=L_ext)
    out = _outputs(grid, props, fld, bottom, L_ext, h_box, fld.projected_area)
    j_1d = props.rho_D * (props.X_wet - props.X_top) / DELTA_EFF
    out.update(j_1d=j_1d, rel_err_1d=(out["j_projected"] - j_1d) / j_1d,
               X_profile=fld.X[0, 0, :].tolist())
    return out, fld


def solve_pore(h_box, shift=0.0):
    """Reference 2: only the rasterised opening wet, the rest of the plane a no-flux solid.
    ``shift`` = A_VOXEL / 2 is the supplementary rasterisation diagnostic (11 x 11 cells)."""
    grid, openings, L_ext = build(h_box, shift=shift)
    _, props = gas_properties()
    state = np.full(grid.shape, GAS, int)
    bottom = openings.bottom_condition(np.array([props.X_wet]), np.zeros(grid.shape[:2], bool))
    fld = solve_vapor(grid, state, props, bottom=bottom, top_external_length=L_ext)
    A_wet = float(openings.voxel_area()[0])
    out = _outputs(grid, props, fld, bottom, L_ext, h_box, A_wet)
    ras = rasterisation(h_box, grid.a, shift)
    out.update(A_phys=ras["A_phys"], n_cells=ras["n_cells"], side_raster=ras["side_raster"],
               area_rel_err=ras["area_rel_err"], shift=shift,
               j_wet_physical_area=out["E"] / ras["A_phys"],
               # lateral uniformity at the box top: the Robin segment is exact only if this is ~0
               top_nonuniformity=float(np.ptp(fld.X[..., -1]) / (props.X_wet - props.X_top)),
               # the flux of the wet cells: edge enhancement of an isolated patch
               base_flux_max_over_mean=float(fld.base_flux[openings.foot[0]].max()
                                             / fld.base_flux[openings.foot[0]].mean()))
    return out, fld
