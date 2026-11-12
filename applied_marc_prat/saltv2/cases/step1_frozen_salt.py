"""Step 1 of the surface-first verification ladder: evaporation from a PRESCRIBED wet external
geometry -- one frozen square-prism salt pillar over the single surface pore of Step 0.

Physical model.  The Step 0 gas problem (``step0_gas``: T = 22 C, RH = 0.07, saturated brine
C = 0.264, 540 um square opening in a 700 um periodic cell, far field at DELTA_EFF = 20.6 mm) with a
square prism of the opening's footprint and prescribed height h standing on the opening.  Its
exposed faces are either WET (Dirichlet X = a_w(C_sat) W_sat, unlimited supply) or INERT (zero
normal flux); the opening under the pillar never evaporates into the gas (it is covered).  Three
prescribed conditions:
    A0  h = 0, the bare opening of Step 0;
    A1  wet top face, inert side faces;
    A2  wet top and wet side faces.
No liquid or salt transport in the pillar, no capillarity, no growth, no heat: this is NOT a model
of the salt material, it measures how a given wet geometry changes the gas-side evaporation.

Numerical model.  The production kernel ``growth.vapor.solve_vapor`` makes EVERY face of a
SALT_WET voxel wet, so it can represent A0 and A2 but not A1.  ``solve_faces`` below is the
case-local face-selective form of the same discretisation: the pillar voxels are removed from the
domain (``fv.assemble_laplace`` with no internal Dirichlet cells, no bottom Dirichlet), and the
Dirichlet half-cell conductance of ``fv._half_scale`` -- coefficient x face area / (a / 2) -- is
added only on the selected (gas cell, direction) faces.  ``cross_check`` validates it: with the top
and side faces selected it must reproduce ``solve_vapor`` on the SALT_WET pillar, and with h = 0 the
bare opening through ``solve_vapor(bottom=...)``, both to round-off.  Nothing in ``saltv2.growth``
is changed.

Grid.  A_VOXEL = 20 um resolves the 540 um opening (27 cells) and the 700 um period (35 cells)
exactly, and every height of HEIGHTS is a multiple of it; A_VOXEL_COARSE = 50 um is the Step 0
grid (10 x 10 cells = 500 um) for the regression and the limited grid comparison.

    from saltv2.cases.step1_frozen_salt import run_case, cross_check
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, cg, spsolve

from ..coupling.interface import PoreOpenings
from ..growth import fv
from ..growth.grid import GAS, SALT_WET, Grid
from ..growth.vapor import solve_vapor
from .step0_gas import DELTA_EFF, L_BOX, PORE_SIDE, gas_properties

# ---- case values ---------------------------------------------------------------------------
A_VOXEL = 20e-6                     # [m, N] primary grid: 35 x 35 lateral, 27 x 27 opening
A_VOXEL_COARSE = 50e-6              # [m, N] the Step 0 grid
H_BOX = 2.0e-3                      # [m, N] explicit gas-domain height; L_ext = DELTA_EFF - H_BOX
H_BOX_CHECK = 4.0e-3                # [m, N] the one domain-truncation check (A2, h = 0.5 mm)
HEIGHTS = (0.0, 0.1e-3, 0.3e-3, 0.5e-3)     # prescribed pillar heights [m, N]
CONDITIONS = ("A0", "A1", "A2")
CG_RTOL = 1e-12


def build(h_box=H_BOX, a=A_VOXEL):
    """Grid, the opening footprint (lateral bool) and the external segment length."""
    n_lat = L_BOX / a; n_z = h_box / a
    if abs(n_lat - round(n_lat)) > 1e-9 or abs(n_z - round(n_z)) > 1e-9:
        raise ValueError("L_BOX and h_box must be integer multiples of the voxel size")
    grid = Grid((int(round(n_lat)), int(round(n_lat)), int(round(n_z))), a)
    openings = PoreOpenings(grid, [0.5 * L_BOX], [0.5 * L_BOX], [0.5 * PORE_SIDE], L_BOX)
    L_ext = DELTA_EFF - h_box
    if L_ext < 0.0:
        raise ValueError("the box is taller than the effective transfer distance")
    return grid, openings.foot[0], L_ext


def pillar_mask(grid, foot, h):
    """Voxels of the prism: the footprint, layers k < h / a (h must be a multiple of a)."""
    n_h = h / grid.a
    if abs(n_h - round(n_h)) > 1e-9:
        raise ValueError("the pillar height must be a multiple of the voxel size")
    n_h = int(round(n_h))
    if n_h >= grid.n_z:
        raise ValueError("the pillar reaches the box top")
    mask = np.zeros(grid.shape, bool)
    mask[..., :n_h] = foot[..., None]
    return mask


def wet_faces(grid, foot, h, sides):
    """The wet (gas cell, direction) faces as {(axis, step): lateral-or-3-D bool mask} in the
    convention of ``VaporField.face_fluxes``: the mask marks the GAS cell, the wet surface is its
    face towards ``step`` along ``axis``.  Top face: the gas cell right above the pillar (or, for
    h = 0, the gas cell on the opening -- its bottom face IS the opening).  Side faces: the gas
    cells laterally next to a pillar voxel.  Each physical face appears exactly once."""
    pillar = pillar_mask(grid, foot, h)
    n_h = int(round(h / grid.a))
    faces = {}
    top = np.zeros(grid.shape, bool); top[..., n_h] = foot
    faces[(grid.z_axis, -1)] = top
    if sides and n_h > 0:
        for axis in grid.lateral_axes:
            for step in (+1, -1):
                m = grid.shift(pillar, axis, step, False) & ~pillar
                faces[(axis, step)] = m
    return pillar, faces


def _solve(A, b, solver, rtol):
    """Returns (x, CG iteration count; 0 for the direct solver)."""
    n_iter = 0
    if solver == "direct":
        x = spsolve(A, b)
    elif solver == "cg":
        diag = A.diagonal()
        M = LinearOperator(A.shape, lambda v: v / diag)
        count = [0]
        x, info = cg(A, b, M=M, rtol=rtol, maxiter=200000, callback=lambda xk: count.__setitem__(0, count[0] + 1))
        n_iter = count[0]
        if info != 0:
            raise RuntimeError(f"CG did not converge (info={info})")
    else:
        raise ValueError("solver must be 'direct' or 'cg'")
    if not np.all(np.isfinite(x)):
        raise RuntimeError("non-finite solution")
    return x, n_iter


def solve_faces(grid, props, pillar, faces, L_ext, solver="direct", rtol=CG_RTOL):
    """Steady vapour diffusion on the gas cells outside ``pillar`` with Dirichlet X_wet on the
    selected faces only, zero flux on every other solid face, the 1-D segment L_ext at the top.
    Returns (outputs dict, X field with NaN in the pillar)."""
    a = grid.a
    rho_D = fv.as_field(props.rho_D, grid.shape, "rho_D")
    active = fv.top_connected(grid, ~pillar)
    if not np.array_equal(active, ~pillar):
        raise RuntimeError("gas cells disconnected from the top")
    system = fv.assemble_laplace(grid, rho_D, active, np.zeros(grid.shape, bool), props.X_wet,
                                 props.X_wet, props.X_top, bottom_mask=np.zeros(grid.shape[:-1], bool),
                                 top_external_length=L_ext)
    # face-selective Dirichlet terms: conductance = rho_D * a^2 / (a / 2), as fv._half_scale
    add = np.zeros(system.matrix.shape[0])
    for (axis, step), m in faces.items():
        if np.any(m & ~active):
            raise ValueError("a wet face must belong to an active gas cell")
        np.add.at(add, system.index[m], (2.0 * rho_D * a)[m])
    A = (system.matrix + sp.diags(add)).tocsr()
    b = system.rhs + add * props.X_wet
    x, n_iter = _solve(A, b, solver, rtol)
    residual = float(np.linalg.norm(A @ x - b) / np.linalg.norm(b))
    X = np.full(grid.shape, np.nan); X[active] = x
    # fluxes on the selected faces (positive = evaporation), per face set
    E_face, j_face = {}, {}
    for key, m in faces.items():
        j = np.zeros(grid.shape); j[m] = (2.0 * rho_D / a)[m] * (props.X_wet - X[m])
        E_face[key] = float(j[m].sum() * a * a); j_face[key] = j
    top_cells = active[..., -1]
    E_far = float(((rho_D[..., -1] / (0.5 * a + L_ext)) * (X[..., -1] - props.X_top))[top_cells].sum() * a * a)
    E_top = E_face[(grid.z_axis, -1)]
    E_side = sum(v for k, v in E_face.items() if k[0] != grid.z_axis)
    E_total = E_top + E_side
    j_top = j_face[(grid.z_axis, -1)]; m_top = faces[(grid.z_axis, -1)]
    side_vals = np.concatenate([j_face[k][faces[k]] for k in faces if k[0] != grid.z_axis]) if E_side else np.array([])
    out = dict(E_total=E_total, E_top=E_top, E_side=E_side, E_farfield=E_far,
               imbalance=(E_total - E_far) / E_far, residual=residual,
               n_unknown=int(A.shape[0]), solver=solver, cg_iterations=n_iter, cg_rtol=rtol if solver == "cg" else None,
               top_max_over_mean=float(j_top[m_top].max() / j_top[m_top].mean()),
               side_max_over_mean=float(side_vals.max() / side_vals.mean()) if side_vals.size else None,
               n_side_faces=int(side_vals.size))
    return out, X, j_face


def run_case(cond, h, a=A_VOXEL, h_box=H_BOX, solver="direct"):
    """One prescribed condition at height h [m].  Returns (summary, X, face fluxes, masks)."""
    if cond not in CONDITIONS:
        raise ValueError("cond must be A0, A1 or A2")
    if cond == "A0" and h != 0.0:
        raise ValueError("A0 is the bare opening, h = 0")
    grid, foot, L_ext = build(h_box, a)
    _, props = gas_properties()
    pillar, faces = wet_faces(grid, foot, h, sides=(cond == "A2"))
    out, X, j_face = solve_faces(grid, props, pillar, faces, L_ext, solver)
    n_cells = int(foot.sum()); side_raster = np.sqrt(n_cells) * a
    A_open_phys = PORE_SIDE ** 2; A_open_raster = n_cells * a * a
    A_top_wet = A_open_raster                                       # voxel top = voxel footprint
    A_side_phys = 4.0 * PORE_SIDE * h if cond == "A2" else 0.0
    A_side_raster = 4.0 * side_raster * h if cond == "A2" else 0.0  # = n_side_faces a^2
    if cond == "A2" and abs(A_side_raster - out["n_side_faces"] * a * a) > 1e-18:
        raise RuntimeError("side-face count does not match the prism")
    out.update(cond=cond, h=h, a=a, h_box=h_box, L_ext=L_ext, delta_eff=DELTA_EFF, shape=list(grid.shape),
               A_projected=L_BOX ** 2, A_open_phys=A_open_phys, A_open_raster=A_open_raster,
               A_top_wet_phys=A_open_phys, A_top_wet=A_top_wet,
               A_side_wet_phys=A_side_phys, A_side_wet=A_side_raster,
               A_total_wet_phys=A_open_phys + A_side_phys, A_total_wet=A_top_wet + A_side_raster,
               j_projected=out["E_total"] / L_BOX ** 2, j_top_mean=out["E_top"] / A_top_wet,
               j_side_mean=(out["E_side"] / A_side_raster) if A_side_raster else None,
               j_wet_mean=out["E_total"] / (A_top_wet + A_side_raster),
               X_wet=props.X_wet, X_top=props.X_top, rho_D=float(props.rho_D),
               # 1-D planar control: a fully wet plane raised by h sees delta_eff - h
               ratio_1d_control=DELTA_EFF / (DELTA_EFF - h))
    return out, X, j_face, dict(pillar=pillar, faces=faces, foot=foot)


def cross_check(h=0.3e-3, a=A_VOXEL_COARSE, h_box=H_BOX, solver="direct"):
    """The face-selective assembly against the production kernel: (i) A2 = solve_vapor on a
    SALT_WET pillar (every face wet, base covered); (ii) h = 0 = solve_vapor(bottom=opening).
    Returns the relative differences of E and of the field."""
    grid, foot, L_ext = build(h_box, a)
    _, props = gas_properties()
    out = {}
    # (i) A2
    pillar, faces = wet_faces(grid, foot, h, sides=True)
    o2, X2, _ = solve_faces(grid, props, pillar, faces, L_ext, solver)
    state = np.where(pillar, SALT_WET, GAS)
    f = solve_vapor(grid, state, props, bottom=(np.zeros(grid.shape[:2], bool), np.zeros(grid.shape[:2])),
                    top_external_length=L_ext, solver=solver, rtol=CG_RTOL)
    sel = ~pillar
    out["A2_E_rel"] = (o2["E_total"] - f.M_salt) / f.M_salt
    out["A2_X_rel_max"] = float(np.max(np.abs(X2[sel] - f.X[sel])) / props.X_wet)
    out["A2_E_base_kernel"] = f.M_base                      # must be 0: the opening is covered
    out["A2_E_top_kernel"] = float((f.face_flux(grid.z_axis, -1) * grid.a ** 2).sum())
    out["A2_E_top_rel"] = (o2["E_top"] - out["A2_E_top_kernel"]) / out["A2_E_top_kernel"]
    # (ii) h = 0
    pillar0, faces0 = wet_faces(grid, foot, 0.0, sides=False)
    o0, X0, _ = solve_faces(grid, props, pillar0, faces0, L_ext, solver)
    f0 = solve_vapor(grid, np.full(grid.shape, GAS, int), props,
                     bottom=(foot, np.full(grid.shape[:2], props.X_wet)), top_external_length=L_ext, solver=solver, rtol=CG_RTOL)
    out["A0_E_rel"] = (o0["E_total"] - f0.M_base) / f0.M_base
    out["A0_X_rel_max"] = float(np.max(np.abs(X0 - f0.X)) / props.X_wet)
    out["A0_E_kernel"] = f0.M_base
    return out
