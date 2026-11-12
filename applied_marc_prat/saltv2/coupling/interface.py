"""Substrate pore network <-> voxel grid: where each surface pore opens into the gas domain.

Regime: evaporation-wicking (every opening is liquid-filled; its concentration comes from the
network).  One object owns the opening geometry, so that the vapour boundary, the per-pore sum of
the base flux, nucleation and the column-to-pore map all use the same footprint:

    (x_i, y_i, r_i)  --PoreOpenings-->  footprint of pore i on the z = 0 plane
    C_i              --driver-->        X_open,i = a_w(C_i) W_sat      (water activity: driver)
    (footprint, X_open, z = 0 state)  --bottom_condition-->  Dirichlet patches for solve_vapor
    base flux j(x, y)                 --per_pore-->          evaporation of pore i  [kg/s / A_cell]

The water-activity law is ``saltv2.materials.brine.a_w`` and the network is ``saltv2.pnm``; this
module has no physics, only geometry and bookkeeping (moved out of the 3-D driver 2026-09-15,
behaviour bit for bit unchanged, tests/test_regression.py).  It is the first piece of the
interface planned in ideas/plan/08_S8_PNM接口.md.

Opening shape [A, decided 2026-09-15].  A surface pore is the QSPNM cube of side 2 r_i that the
pore network uses (``saltv2.pnm.PoreLayer``: V = (2 r_i)^3, A_open = (2 r_i)^2),
so its opening is the SQUARE top face: the footprint of pore i is the set of lateral cells whose
CENTRE lies within |x - x_i| <= r_i and |y - y_i| <= r_i, periodic in both lateral directions.
The vapour boundary and the supply check (demand / A_open) therefore describe the same opening.
Before this date the footprint was the inscribed disk (area pi r_i^2, pi/4 of A_open).

Resolution.  The rasterised area (``voxel_area``) approaches (2 r_i)^2 only as the voxel shrinks
(each side is off by less than one cell, not monotonically).  With the dgf radii (r = 0.36-0.41
pitch) vpp 6 resolves 0.444 pitch^2 of 0.584 (76 %; the same 4 x 4 cells the disk gave, so vpp-6
runs are unchanged), vpp 12 about 0.56, and vpp 3 fills the whole pitch.
"""
from __future__ import annotations

import numpy as np

from ..growth.grid import Grid


class PoreOpenings:
    """Footprints of the surface pores of a periodic square substrate on a 3-D voxel grid.

    ``x, y`` [m] pore centres, ``r`` [m] half the side of the square opening, ``L_box`` [m]
    lateral period.
    ``foot[i]`` is the lateral boolean footprint of pore i; ``pore_id`` maps each lateral cell to
    the pore that opens there (-1: inert substrate between pores; where two footprints overlap the
    higher index wins, as the map is filled in index order).
    """

    def __init__(self, grid: Grid, x, y, r, L_box: float):
        if grid.ndim != 3:
            raise ValueError("PoreOpenings is defined on a 3-D grid (two lateral axes + z)")
        n_x, n_y = grid.shape[:2]
        if n_x != n_y:
            raise ValueError("PoreOpenings expects a square lateral grid")
        x = np.asarray(x, float); y = np.asarray(y, float); r = np.asarray(r, float)
        if not (x.shape == y.shape == r.shape) or x.ndim != 1:
            raise ValueError("x, y and r must be 1-D arrays of the same length")
        self.grid = grid
        self.n_pores = int(len(r))
        self.x = x.copy()
        self.y = y.copy()
        self.r = r.copy()
        self.L_box = float(L_box)
        n_lat = n_x
        a = grid.a
        xc = (np.arange(n_lat) + 0.5) * a
        X, Y = np.meshgrid(xc, xc, indexing="ij")
        foot = np.zeros((self.n_pores, n_lat, n_lat), bool)
        for i in range(self.n_pores):
            dx = np.abs(X - x[i]); dx = np.minimum(dx, L_box - dx)
            dy = np.abs(Y - y[i]); dy = np.minimum(dy, L_box - dy)
            foot[i] = (dx <= r[i]) & (dy <= r[i])
        pore_id = np.full((n_lat, n_lat), -1, int)
        for i in range(self.n_pores):
            pore_id[foot[i]] = i
        self.foot = foot
        self.pore_id = pore_id
        self._on_pore = pore_id >= 0
        self._pid_flat = pore_id[self._on_pore]

    def nucleation_seed(self, i: int) -> np.ndarray:
        """One deterministic bottom-layer cell for a local nucleation seed.

        A surface pore is a coarse PNM object, whereas a first crystal is
        localized.  The driver therefore uses this one-cell mask at nucleation
        rather than turning the full pore footprint into wet salt at once.
        The selected cell is the footprint cell whose centre is closest to the
        physical pore centre, with periodic distances.  It is always a subset
        of ``foot[i]`` and ties are resolved reproducibly by NumPy's row-major
        ordering.
        """
        if not 0 <= int(i) < self.n_pores:
            raise IndexError("pore index out of range")
        i = int(i)
        n_lat = self.grid.shape[0]
        xc = (np.arange(n_lat) + 0.5) * self.grid.a
        X, Y = np.meshgrid(xc, xc, indexing="ij")
        dx = np.abs(X - self.x[i]); dx = np.minimum(dx, self.L_box - dx)
        dy = np.abs(Y - self.y[i]); dy = np.minimum(dy, self.L_box - dy)
        dist2 = dx * dx + dy * dy
        dist2 = np.where(self.foot[i], dist2, np.inf)
        out = np.zeros_like(self.foot[i])
        out[np.unravel_index(np.argmin(dist2), dist2.shape)] = True
        return out

    def per_pore(self, lateral) -> np.ndarray:
        """Sum a lateral field over each pore's footprint -> array over pores."""
        return np.bincount(self._pid_flat, weights=np.asarray(lateral, float)[self._on_pore],
                           minlength=self.n_pores)

    def bottom_condition(self, X_open, blocked):
        """Dirichlet patches of the pore openings for ``solve_vapor(bottom=...)``.

        ``X_open`` [-] (one value per pore) is the vapour mass fraction the brine in that opening
        imposes; ``blocked`` (lateral bool) marks the cells whose opening does NOT evaporate
        through the bottom.  The E1 driver blocks openings under a WET z = 0 voxel: that brine is
        liquid-connected to the structure and evaporates through the structure's own faces.  An
        opening under a SEALED voxel stays open; solve_vapor gives it a base flux only where that
        voxel conducts.
        """
        X_open = np.asarray(X_open)
        if X_open.shape != (self.n_pores,):
            raise ValueError("X_open must hold one value per pore")
        lateral = self.foot.shape[1:]
        bmask = np.zeros(lateral, bool); bval = np.zeros(lateral)
        for i in range(self.n_pores):
            bmask |= self.foot[i]; bval[self.foot[i]] = X_open[i]
        bmask &= ~np.asarray(blocked, bool)
        return bmask, bval

    def voxel_area(self) -> np.ndarray:
        """Area [m^2] of each rasterised footprint (what the vapour solve actually sees)."""
        return self.foot.sum(axis=(1, 2)) * self.grid.face_area
