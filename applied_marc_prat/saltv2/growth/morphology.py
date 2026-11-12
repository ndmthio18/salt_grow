"""Morphology observables of a salt set (T section 19), dimension-agnostic.

All functions take a boolean ``salt`` mask shaped like the Grid and return
plain numbers or arrays; nothing here depends on physics.

Definitions
-----------
label_clusters      face-connected components, periodic along lateral axes
spanning            a cluster whose footprint, in the NON-periodic sense,
                    touches both ends of a lateral axis (2-D: covers the
                    whole line, i.e. theta = 1 -- see plan/11 section 2.3)
height_stats        mean and rms of column heights over covered columns, and
                    the maximum height
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .grid import Grid


def label_clusters(salt: np.ndarray, grid: Grid):
    """(labels, n): face-connected components with lateral periodicity.
    Labels 1..n, 0 = background."""
    struct = ndimage.generate_binary_structure(salt.ndim, 1)
    lab, n = ndimage.label(salt, structure=struct)
    if n <= 1:
        return lab, n
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for ax in grid.lateral_axes:
        lo = np.take(lab, 0, axis=ax)
        hi = np.take(lab, -1, axis=ax)
        touch = (lo > 0) & (hi > 0)
        for a, b in zip(lo[touch].ravel(), hi[touch].ravel()):
            ra, rb = find(int(a)), find(int(b))
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(a) for a in range(n + 1)])
    uniq = np.unique(roots[1:])
    remap = np.zeros(n + 1, int)
    remap[uniq] = np.arange(1, len(uniq) + 1)
    return remap[roots[lab]], int(len(uniq))


def n_clusters(salt: np.ndarray, grid: Grid) -> int:
    return label_clusters(salt, grid)[1]


def spanning(salt: np.ndarray, grid: Grid) -> bool:
    """True if some face-connected cluster (non-periodic labelling) touches
    both ends of at least one lateral axis."""
    struct = ndimage.generate_binary_structure(salt.ndim, 1)
    lab, n = ndimage.label(salt, structure=struct)
    for ax in grid.lateral_axes:
        lo = set(np.unique(np.take(lab, 0, axis=ax))) - {0}
        hi = set(np.unique(np.take(lab, -1, axis=ax))) - {0}
        if lo & hi:
            return True
    return False


def height_stats(salt: np.ndarray, grid: Grid):
    """dict(h_mean, h_rms, h_max) [m]; means over covered columns."""
    h = grid.column_max_height(salt)
    covered = h > 0
    if not covered.any():
        return {"h_mean": 0.0, "h_rms": 0.0, "h_max": 0.0}
    hc = h[covered]
    return {"h_mean": float(hc.mean()),
            "h_rms": float(np.sqrt(np.mean((hc - hc.mean()) ** 2))),
            "h_max": float(h.max())}

