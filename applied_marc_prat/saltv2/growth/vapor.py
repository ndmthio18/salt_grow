"""Quasi-steady vapour diffusion on a prescribed salt geometry.

The geometry is fixed for one solve: wet salt is a Dirichlet surface at ``X_wet``, sealed
salt may join the domain as a porous conductor (``porous=``), bare pore openings are given
per column (``bottom=``) and the far field sits ``top_external_length`` beyond the box top.
The solver returns the vapour mass-fraction field and face-resolved evaporation fluxes; it
does not advance salt or invoke the Prat growth law.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg, spsolve

from . import fv
from .fv import as_field, assemble_laplace, top_connected
from .grid import GAS, SALT_WET, Grid


@dataclass(frozen=True)
class VaporProperties:
    """Gas-side properties in the project's vapour mass-fraction convention."""

    rho_D: float | np.ndarray = 3.0e-5  # rho_gas * D_v [kg m^-1 s^-1]
    X_wet: float = 0.0198               # saturated NaCl solution at 21 C
    X_top: float = 0.00185              # approximately RH = 7 %

    def __post_init__(self):
        if not np.isfinite(self.X_wet) or not np.isfinite(self.X_top):
            raise ValueError("vapour boundary values must be finite")
        if self.X_wet < self.X_top:
            raise ValueError("X_wet must be greater than or equal to X_top")


@dataclass(frozen=True)
class VaporField:
    """Vapour solution and conservative flux accounting.

    ``face_fluxes[(axis, step)]`` is an array aligned with gas cells.  At a
    gas cell adjacent to wet salt in that direction it stores the positive
    evaporation flux ``j`` [kg m^-2 s^-1] across that face; elsewhere it is
    zero.  This convention matches ``GeometryEngine.exposed_faces``.
    """

    X: np.ndarray
    active_gas: np.ndarray
    face_fluxes: dict[tuple[int, int], np.ndarray] = field(repr=False)
    base_flux: np.ndarray
    top_flux: np.ndarray
    M_base: float
    M_salt: float
    M_top: float
    imbalance: float
    projected_area: float

    @property
    def M_total(self) -> float:
        """Total evaporation rate from bare substrate plus wet salt [kg/s]."""
        return self.M_base + self.M_salt

    def face_flux(self, axis: int, step: int) -> np.ndarray:
        return self.face_fluxes[(axis, step)]


def solve_vapor(
    grid: Grid,
    state: np.ndarray,
    props: VaporProperties | None = None,
    bottom: tuple[np.ndarray, np.ndarray] | None = None,
    solver: str = "direct",
    x0: np.ndarray | None = None,
    rtol: float = 1e-10,
    top_external_length: float = 0.0,
    porous: tuple[np.ndarray, np.ndarray] | None = None,
) -> VaporField:
    """Solve the vapour problem for a fixed 2-D or 3-D geometry.

    ``bottom=(mask, values)`` (lateral arrays) makes only the marked columns of
    the substrate plane evaporate, at the given vapour mass fractions (bare
    pore openings with their own water activity); the rest of the plane is a
    no-flux solid. Default: the whole plane is wet at ``X_wet``.

    ``top_external_length`` L_ext [m] puts ``X_top`` at distance L_ext beyond
    the box top (series 1-D segment, see ``assemble_laplace``), so that the
    far-field distance delta_eff = box height + L_ext can exceed the box.

    ``porous=(mask, factor)`` adds the masked cells (the sealed crust) to the
    diffusion domain with coefficient ``rho_D * factor`` per cell (``factor``
    is a full-grid array or scalar, the sealed voxel's conductance over D_v).
    They are neither sources nor walls: a wet-salt face next to a porous cell
    evaporates into it, and the top-connectivity is evaluated on gas plus
    porous cells.  Default None reproduces the gas-only problem exactly.
    """
    props = VaporProperties() if props is None else props
    state = np.asarray(state)
    if state.shape != grid.shape:
        raise ValueError("state shape must match grid.shape")

    gas = state == GAS
    wet_salt = state == SALT_WET
    wet_values = np.broadcast_to(np.asarray(props.X_wet, float), grid.shape)
    rho_D = as_field(props.rho_D, grid.shape, "rho_D")
    if porous is None:
        active = top_connected(grid, gas)
    else:
        porous_mask, porous_factor = porous
        porous_mask = np.asarray(porous_mask, dtype=bool)
        if porous_mask.shape != grid.shape:
            raise ValueError("porous mask shape must match grid.shape")
        if np.any(porous_mask & (gas | wet_salt)):
            raise ValueError("porous cells must be neither GAS nor SALT_WET")
        factor = np.broadcast_to(np.asarray(porous_factor, float), grid.shape)
        if not np.all(np.isfinite(factor[porous_mask])) or np.any(factor[porous_mask] <= 0.0):
            raise ValueError("porous factor must be finite and strictly positive on the mask")
        rho_D = np.where(porous_mask, rho_D * factor, rho_D)
        active = top_connected(grid, gas | porous_mask)
    bottom_mask = bottom_values = None
    if bottom is not None:
        bottom_mask, bottom_values = bottom
        bottom_values = np.broadcast_to(np.asarray(bottom_values, float), grid.shape[:-1])
    system = assemble_laplace(
        grid,
        rho_D,
        active,
        wet_salt,
        wet_values,
        props.X_wet,
        props.X_top,
        bottom_mask=bottom_mask,
        bottom_values=bottom_values,
        top_external_length=top_external_length,
    )
    if solver == "direct":
        solution = spsolve(system.matrix, system.rhs)
    elif solver == "cg":
        # SPD system: conjugate gradients with Jacobi preconditioning; x0 is a
        # full-grid field (previous solution) used as warm start where active.
        A = system.matrix
        diag = A.diagonal()
        M = LinearOperator(A.shape, lambda v: v / diag)
        start = None
        if x0 is not None:
            start = np.asarray(x0, float)[active]
            start = np.where(np.isfinite(start), start, props.X_top)
        solution, info = cg(A, system.rhs, x0=start, M=M, rtol=rtol, maxiter=20000)
        if info != 0:
            raise RuntimeError(f"vapour CG did not converge (info={info})")
    else:
        raise ValueError("solver must be 'direct' or 'cg'")
    if not np.all(np.isfinite(solution)):
        raise RuntimeError("vapour solve produced a non-finite solution")

    # Salt is held at X_wet.  Disconnected gas is NaN by design: it was
    # removed from the PDE and must not silently contribute to a flux.
    X = np.full(grid.shape, np.nan)
    X[wet_salt] = wet_values[wet_salt]
    X[active] = solution

    sizes, area, _ = fv.mesh_metrics(grid)
    a_lat = area[grid.ndim - 1]              # horizontal face area of each column
    dz = sizes[grid.ndim - 1]
    base_flux_full = np.zeros(grid.shape)
    base_cells = np.zeros(grid.shape, dtype=bool)
    base_cells[..., 0] = active[..., 0]
    if bottom_mask is not None:
        base_cells[..., 0] &= bottom_mask
    base_value = props.X_wet if bottom_values is None else np.broadcast_to(bottom_values[..., None], grid.shape)[base_cells]
    base_flux_full[base_cells] = (
        (2.0 * rho_D / np.broadcast_to(dz, grid.shape))[base_cells]
        * (base_value - X[base_cells])
    )
    base_flux = base_flux_full[..., 0]

    top_flux_full = np.zeros(grid.shape)
    top_cells = np.zeros(grid.shape, dtype=bool)
    top_cells[..., -1] = active[..., -1]
    top_flux_full[top_cells] = (
        (rho_D / np.broadcast_to(0.5 * dz + top_external_length, grid.shape))[top_cells]
        * (X[top_cells] - props.X_top)
    )
    top_flux = top_flux_full[..., -1]

    face_fluxes: dict[tuple[int, int], np.ndarray] = {}
    M_salt = 0.0
    for axis, step in grid.face_directions():
        salt_neighbour = grid.shift(wet_salt, axis, step, False)
        interface = active & salt_neighbour
        flux = np.zeros(grid.shape)
        flux[interface] = (
            (2.0 * rho_D / np.broadcast_to(sizes[axis], grid.shape))[interface]
            * (grid.shift(wet_values, axis, step, 0.0)[interface] - X[interface])
        )
        face_fluxes[(axis, step)] = flux
        M_salt += float((flux * np.broadcast_to(area[axis], grid.shape))[interface].sum())

    col_area = np.broadcast_to(a_lat, grid.shape)[..., 0]
    M_base = float((base_flux * col_area).sum())
    M_top = float((top_flux * col_area).sum())
    imbalance = (
        (M_base + M_salt - M_top) / M_top if M_top != 0.0 else 0.0
    )
    return VaporField(
        X=X,
        active_gas=active,
        face_fluxes=face_fluxes,
        base_flux=base_flux,
        top_flux=top_flux,
        M_base=M_base,
        M_salt=M_salt,
        M_top=M_top,
        imbalance=float(imbalance),
        projected_area=float(col_area.sum()),
    )


def solve_vapor_adaptive(grid: Grid, state: np.ndarray, props: VaporProperties, bottom=None,
                         solver: str = "cg", x0=None, rtol: float = 1e-7,
                         top_external_length: float = 0.0, porous=None):
    """``solve_vapor`` on the lowest layers only, the rest as the 1-D external segment.

    Above the salt the field becomes laterally uniform within about one lateral period
    (exponential decay of the lateral modes), so the box can be cut one lateral period plus
    four cells above the highest non-gas voxel and the removed height added to
    ``top_external_length``.  This is the same series-resistance argument that the Robin top
    already uses (tests/test_V1c_solver.py).
    Returned arrays are embedded in the full grid: fluxes above the cut are zero (no salt
    there), ``active_gas`` above the cut is the gas mask, and ``X`` above the cut is the
    1-D linear profile (used only as the next warm start).  The solve cost scales with the
    number of active layers and, for CG, roughly with their number again through the
    iteration count, so a structure much lower than the box is solved many times faster.
    """
    state = np.asarray(state)
    non_gas = state != GAS
    if porous is not None:
        non_gas = non_gas | np.asarray(porous[0], bool)
    k_top = int(np.max(np.nonzero(non_gas.any(axis=tuple(range(grid.ndim - 1))))[0])) if non_gas.any() else -1
    margin_cells = max(grid.shape[:-1]) + 4
    nz_act = min(grid.n_z, k_top + 1 + margin_cells)
    if nz_act >= grid.n_z:
        return solve_vapor(grid, state, props, bottom=bottom, solver=solver, x0=x0, rtol=rtol,
                           top_external_length=top_external_length, porous=porous)
    sub = Grid(grid.shape[:-1] + (nz_act,), grid.a)
    cut = (nz_act - 1)
    sl = (Ellipsis, slice(0, nz_act))
    sub_porous = None if porous is None else (np.asarray(porous[0], bool)[sl], np.asarray(porous[1], float)[sl])
    sub_x0 = None if x0 is None else np.asarray(x0, float)[sl]
    L_ext = top_external_length + (grid.n_z - nz_act) * grid.a
    f = solve_vapor(sub, state[sl], props, bottom=bottom, solver=solver, x0=sub_x0, rtol=rtol,
                    top_external_length=L_ext, porous=sub_porous)
    X = np.full(grid.shape, np.nan)
    X[sl] = f.X
    # 1-D linear continuation to X_top over the external length (warm start only)
    z_above = (np.arange(nz_act, grid.n_z) - cut) * grid.a
    X_cut = np.where(np.isfinite(f.X[..., -1]), f.X[..., -1], props.X_top)
    frac = z_above / (L_ext + 0.5 * grid.a)
    X[..., nz_act:] = X_cut[..., None] + (props.X_top - X_cut[..., None]) * frac
    active = np.zeros(grid.shape, bool); active[sl] = f.active_gas
    active[..., nz_act:] = (state[..., nz_act:] == GAS)
    faces = {}
    for key, arr in f.face_fluxes.items():
        full = np.zeros(grid.shape); full[sl] = arr; faces[key] = full
    return VaporField(X=X, active_gas=active, face_fluxes=faces, base_flux=f.base_flux, top_flux=f.top_flux,
                      M_base=f.M_base, M_salt=f.M_salt, M_top=f.M_top, imbalance=f.imbalance,
                      projected_area=f.projected_area)

