"""Finite-volume building blocks shared by the transport solvers.

V1a needs only a steady scalar diffusion equation on the gas cells, but the
assembly is kept independent of vapour properties so that later Darcy and
transport modules can reuse the same lattice conventions.

All lateral axes are periodic and the last axis is vertical.  Cell-centre to
cell-centre faces use the harmonic mean of the two diffusivities.  A
Dirichlet value imposed at a cell face uses the gas-side half-cell distance,
and therefore twice the corresponding one-cell conductance.

Meshes are uniform (``Grid.a`` in every direction); ``mesh_metrics`` returns the
sizes and face areas as scalar-valued broadcast arrays so the assembly is written
once in metric form.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from .grid import Grid


@dataclass(frozen=True)
class LinearSystem:
    """Sparse diffusion system and the full-grid map of its unknowns."""

    matrix: sp.csr_matrix
    rhs: np.ndarray
    index: np.ndarray
    active: np.ndarray


def mesh_metrics(grid: Grid):
    """Cell sizes and face metrics of the uniform mesh.

    Returns ``(sizes, area, dist)`` where for each axis
      ``sizes[k]``  cell size along axis k [m] as a ``(1,)*d`` broadcast array,
      ``area[k]``   face area normal to axis k [m^(d-1)], same shape,
      ``dist[k]``   the cell sizes again (the caller forms centre-to-centre distances).
    """
    d = grid.ndim
    sizes = [np.full((1,) * d, grid.a) for _ in range(d)]
    area = [np.full((1,) * d, grid.a ** (d - 1)) for _ in range(d)]
    return sizes, area, sizes


def as_field(value, shape, name: str) -> np.ndarray:
    """Return a finite, strictly positive full-grid coefficient field."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        arr = np.full(shape, float(arr))
    else:
        arr = np.broadcast_to(arr, shape).astype(float, copy=False)
    if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
        raise ValueError(f"{name} must be finite and strictly positive")
    return arr


def top_connected(grid: Grid, open_mask: np.ndarray) -> np.ndarray:
    """Open cells face-connected to the top boundary.

    Gas pockets sealed by salt do not participate in the steady diffusion
    problem; keeping them would add unconstrained rows and make the matrix
    singular.  Periodic lateral axes, face connectivity: scipy.ndimage.label
    plus a union-find merge across the periodic seams (tests/test_V1c_solver.py
    checks it against a plain breadth-first search).
    """
    from scipy import ndimage
    open_mask = np.asarray(open_mask, dtype=bool)
    if open_mask.shape != grid.shape:
        raise ValueError("open_mask shape must match grid.shape")
    lab, n = ndimage.label(open_mask, structure=ndimage.generate_binary_structure(grid.ndim, 1))
    if n == 0:
        return np.zeros(grid.shape, dtype=bool)
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for ax in grid.lateral_axes:
        lo = np.take(lab, 0, axis=ax); hi = np.take(lab, -1, axis=ax)
        touch = (lo > 0) & (hi > 0)
        if touch.any():
            for a, b in np.unique(np.stack([lo[touch], hi[touch]], 1), axis=0):
                ra, rb = find(int(a)), find(int(b))
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(a) for a in range(n + 1)])
    top_roots = np.unique(roots[np.take(lab, -1, axis=grid.z_axis)])
    top_roots = top_roots[top_roots > 0]
    keep = np.isin(roots, top_roots)
    keep[0] = False
    return keep[lab]



def assemble_laplace(
    grid: Grid,
    coefficient,
    active: np.ndarray,
    internal_dirichlet: np.ndarray,
    internal_value: float,
    bottom_value: float,
    top_value: float,
    bottom_mask: np.ndarray | None = None,
    bottom_values: np.ndarray | None = None,
    top_external_length: float = 0.0,
) -> LinearSystem:
    """Assemble ``div(coefficient grad x) = 0`` on active cells.

    ``internal_dirichlet`` marks blocked cells whose faces carry a prescribed
    value (wet salt in V1a).  Other inactive neighbours are no-flux walls.
    The bottom and top values apply at the corresponding domain faces.

    Optional ``bottom_mask`` / ``bottom_values`` (lateral arrays) restrict the
    bottom Dirichlet condition to the marked columns with per-column values
    (bare pore openings); unmarked bottom faces are no-flux (solid substrate).
    Both None reproduces the fully wet bottom plane at ``bottom_value``.

    ``top_external_length`` L_ext >= 0 places ``top_value`` at distance L_ext
    beyond the top face instead of on it (Robin form: the 1-D diffusion segment
    between the box top and the far field, of the same coefficient, is added in
    series). Exact when the field is laterally uniform at the box top, i.e.
    when the box top lies more than about one lateral period above the salt.
    """
    active = np.asarray(active, dtype=bool)
    internal_dirichlet = np.asarray(internal_dirichlet, dtype=bool)
    if active.shape != grid.shape or internal_dirichlet.shape != grid.shape:
        raise ValueError("mask shapes must match grid.shape")
    if np.any(active & internal_dirichlet):
        raise ValueError("active and internal_dirichlet masks must be disjoint")

    coeff = as_field(coefficient, grid.shape, "coefficient")
    internal_values = np.broadcast_to(np.asarray(internal_value, float), grid.shape)
    if not np.all(np.isfinite(internal_values[internal_dirichlet])):
        raise ValueError("internal Dirichlet values must be finite")
    n_unknown = int(active.sum())
    if n_unknown == 0:
        raise ValueError("no open cell is connected to the top boundary")

    index = np.full(grid.shape, -1, dtype=np.int64)
    index[active] = np.arange(n_unknown)
    diagonal = np.zeros(n_unknown)
    rhs = np.zeros(n_unknown)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    values: list[np.ndarray] = []
    sizes, area, _ = mesh_metrics(grid)

    def _half_scale(axis):
        """area / (half cell size) along ``axis`` -- a Dirichlet value on a cell face."""
        return area[axis] / (0.5 * sizes[axis])

    def _face_scale(axis, step):
        """area / (centre-to-centre distance) for the face towards ``step`` along ``axis``."""
        return area[axis] / sizes[axis]

    for axis, step in grid.face_directions():
        neighbour_active = grid.shift(active, axis, step, False)
        pair = active & neighbour_active
        if pair.any():
            # The fill value is irrelevant outside ``pair``; a scalar keeps
            # Grid.shift well-defined at the non-periodic vertical boundary.
            neighbour_coeff = grid.shift(coeff, axis, step, 1.0)
            conductance = (
                2.0 * coeff * neighbour_coeff / (coeff + neighbour_coeff)
            ) * _face_scale(axis, step)
            row = index[pair]
            col = grid.shift(index, axis, step, -1)[pair]
            val = conductance[pair]
            np.add.at(diagonal, row, val)
            rows.append(row)
            cols.append(col)
            values.append(-val)

        neighbour_dirichlet = grid.shift(
            internal_dirichlet, axis, step, False
        )
        boundary = active & neighbour_dirichlet
        if boundary.any():
            row = index[boundary]
            conductance = (coeff * _half_scale(axis))[boundary]
            np.add.at(diagonal, row, conductance)
            neighbour_value = grid.shift(internal_values, axis, step, 0.0)
            np.add.at(rhs, row, conductance * neighbour_value[boundary])

    bottom = np.zeros(grid.shape, dtype=bool)
    bottom[..., 0] = active[..., 0]
    if bottom_mask is not None:
        bottom_mask = np.asarray(bottom_mask, dtype=bool)
        if bottom_mask.shape != grid.shape[:-1]:
            raise ValueError("bottom_mask must be a lateral array")
        bottom[..., 0] &= bottom_mask
    if bottom.any():
        row = index[bottom]
        conductance = (coeff * _half_scale(grid.ndim - 1))[bottom]
        if bottom_values is None:
            values_here = bottom_value
        else:
            values_here = np.broadcast_to(bottom_values, grid.shape[:-1])[bottom[..., 0]]
        np.add.at(diagonal, row, conductance)
        np.add.at(rhs, row, conductance * values_here)

    if not np.isfinite(top_external_length) or top_external_length < 0.0:
        raise ValueError("top_external_length must be finite and >= 0")
    top = np.zeros(grid.shape, dtype=bool)
    top[..., -1] = active[..., -1]
    if top.any():
        row = index[top]
        # half-cell distance to the top face plus the external 1-D segment
        conductance = (coeff * area[grid.ndim - 1]
                       / (0.5 * sizes[grid.ndim - 1] + top_external_length))[top]
        np.add.at(diagonal, row, conductance)
        np.add.at(rhs, row, conductance * top_value)

    rows.append(np.arange(n_unknown))
    cols.append(np.arange(n_unknown))
    values.append(diagonal)
    matrix = sp.csr_matrix(
        (np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n_unknown, n_unknown),
    )
    return LinearSystem(matrix, rhs, index, active)
