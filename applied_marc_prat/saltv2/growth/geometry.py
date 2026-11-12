"""V0 geometry engine: occupancy growth of a salt set with a prescribed
normal velocity on every exposed face (T section 18).

No physics lives here. The caller provides ``V_n`` per exposed face (m/s);
the engine converts it into solid-volume increments, accumulates an
occupancy fraction f in gas cells adjacent to salt, flips cells whose f
reaches 1, and redistributes any overflow to neighbouring candidate cells so
that the injected volume is conserved to round-off.

Time stepping is a fixed (or adaptively capped) dt, not the one-cell event
step of salt_deterministic.py: V2 couples this to a transient substrate
transport problem, which needs a common dt (plan/11 section 2.5).

P1a of plan/17 changed three things in this file; the acceptance tests are in
tests/test_geometry_multiface.py:

1. The growth step limit is the total occupancy gain of a candidate CELL,
   ``dt_grow = f_max / max_i (Vdot_i / V_cell)``, not ``cfl * a / max_face
   V_n``.  A cell fed by several faces gains occupancy faster than any one of
   its faces implies; the old face-based limit let it overshoot.  Use
   ``dt_occupancy`` to obtain it before combining with the other limiters.
2. Overflow of a cell that flips is redistributed to THAT cell's own gas
   neighbours only.  It used to be pooled over all flipped cells and shared
   equally over all receivers, which transferred solid volume between
   crystal clusters that are not in contact.
3. ``_flip_and_redistribute`` raises instead of silently leaving cells above
   f = 1 when the pass limit is reached.

The redistribution is the only lateral coupling in the lattice growth and its
magnitude scales with dt, so it acted as a dt-dependent regulariser (README
E1 section 8).  It is not a substitute for the explicit regularisation length
of P1b -- see ``regularize.smooth_on_mask`` and ``volume_rate(smooth_length=...)``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import GAS, SALT_SEALED, SALT_WET, Grid
from .regularize import smooth_on_mask


@dataclass
class StepReport:
    """What one geometry step did.  ``float(report)`` is the dt used, so
    callers that only want the time step keep working."""

    dt: float
    flipped: np.ndarray            # cells that became SALT_WET in this step
    volume_injected: float         # solid volume handed in during this step
    volume_banked: float           # part of it that found no receiving cell
    passes: int                    # flip/redistribute passes used

    def __float__(self) -> float:
        return float(self.dt)


class GeometryEngine:
    """Occupancy-based interface advance on a Grid.

    Attributes
    ----------
    state : int8 array, GAS / SALT_WET (SALT_SEALED
            untouched here: they are neither growth candidates nor sources)
    f     : float array, occupancy of candidate cells in [0, 1)
    t     : elapsed time [s]
    volume_injected : total solid volume handed to the engine [m^ndim]
    volume_placed   : solid volume represented on the lattice (flips + f)
    """

    def __init__(self, grid: Grid, state: np.ndarray | None = None):
        self.grid = grid
        if state is None:
            state = np.full(grid.shape, GAS, dtype=np.int8)
        self.state = state.astype(np.int8, copy=True)
        self.f = np.zeros(grid.shape)
        self.t = 0.0
        self.volume_injected = 0.0
        self.bank = 0.0            # overflow that found no candidate cell
        # Active window (plan/17 P1a speed-up, README E1 section 9.11): every per-step
        # operation runs on the lowest ``nz_active`` layers only.  The caller sets it a safe
        # margin above the highest salt voxel (salt advances at most one layer per step);
        # None = the whole box.  ``state`` and ``f`` stay full-size, the window is a view.
        self.bank_carry = 0.0    # carry that went to the bank with its volume (see _flip_and_redistribute)
        self.nz_active: int | None = None
        self._grid_w = grid

    # ------------------------------------------------------------------ window
    def set_window(self, nz_active: int | None):
        """Restrict per-step work to the lowest ``nz_active`` layers (None = all)."""
        g = self.grid
        if nz_active is None or nz_active >= g.n_z:
            self.nz_active = None; self._grid_w = g
        else:
            if nz_active < 2:
                raise ValueError("window must hold at least two layers")
            self.nz_active = int(nz_active)
            self._grid_w = Grid(g.shape[:-1] + (self.nz_active,), g.a)
        return self._grid_w

    @property
    def grid_w(self) -> Grid:
        return self._grid_w

    @property
    def state_w(self) -> np.ndarray:
        return self.state if self.nz_active is None else self.state[..., :self.nz_active]

    @property
    def f_w(self) -> np.ndarray:
        return self.f if self.nz_active is None else self.f[..., :self.nz_active]

    def to_full(self, arr_w: np.ndarray, fill=0) -> np.ndarray:
        """Embed a window-shaped array into a full-grid array (``fill`` above the window)."""
        if self.nz_active is None:
            return arr_w
        full = np.full(self.grid.shape, fill, dtype=arr_w.dtype)
        full[..., :self.nz_active] = arr_w
        return full

    # ------------------------------------------------------------------
    def seed(self, lateral_index, height_cells: int = 1, width_cells: int = 1):
        """Place a salt block on the substrate: lateral index (int or tuple),
        ``width_cells`` wide along every lateral axis, ``height_cells`` tall."""
        g = self.grid
        idx = (lateral_index,) if np.isscalar(lateral_index) else tuple(lateral_index)
        sl = []
        for ax, i in enumerate(idx):
            half = width_cells // 2
            cells = [(i - half + w) % g.shape[ax] for w in range(width_cells)]
            sl.append(cells)
        mesh = np.meshgrid(*sl, indexing="ij")
        for k in range(height_cells):
            self.state[tuple(m for m in mesh) + (k,)] = SALT_WET

    def salt_mask(self) -> np.ndarray:
        """Wet salt only: the growing, evaporating set (window-shaped)."""
        return self.state_w == SALT_WET

    def solid_mask(self) -> np.ndarray:
        """Every salt voxel regardless of state (WET / DRY_OPEN / SEALED), window-shaped.
        Drying and sealing move solid between states without changing its
        volume, so the volume ledger must count all three."""
        st = self.state_w
        return (st == SALT_WET) | (st == SALT_SEALED)

    def solid_mask_full(self) -> np.ndarray:
        st = self.state
        return (st == SALT_WET) | (st == SALT_SEALED)

    def exposed_faces(self):
        """List of (axis, step, mask): mask marks GAS cells that have a salt
        neighbour ``step`` cells along ``axis`` -- i.e. the face between them
        is an exposed salt face and the gas cell is a growth candidate."""
        g = self.grid_w
        salt = self.salt_mask()
        gas = self.state_w == GAS
        out = []
        for ax, s in g.face_directions():
            nb_is_salt = g.shift(salt, ax, s, False)
            out.append((ax, s, gas & nb_is_salt))
        return out

    # ------------------------------------------------------------------
    def volume_rate(self, Vn, smooth_length: float = 0.0,
                    smooth_mask: np.ndarray | None = None) -> np.ndarray:
        """Solid volume rate [m^ndim / s] entering each candidate gas cell.

        ``Vn`` is either a scalar, a full array (velocity attributed to the
        gas cell for every exposed face it owns), or a callable
        ``Vn(axis, step, mask) -> array`` returning per-cell velocities for
        that face direction (lets V1 pass face fluxes).

        ``smooth_length`` ell [m] > 0 applies the P1b regularisation: the
        deposition rate is spread over ell along the front, conserving the
        total.  ``smooth_mask`` is the set the rate may live on; it defaults
        to the cells that currently receive a positive rate, i.e. cells whose
        faces are all frozen do not receive smoothed volume.
        """
        g = self.grid_w
        rate = np.zeros(g.shape)
        contributions = []
        for ax, s, mask in self.exposed_faces():
            if callable(Vn):
                v = Vn(ax, s, mask)
            else:
                v = Vn
            face_rate = np.where(mask, v, 0.0) * g.face_area
            rate += face_rate
            if smooth_length > 0.0:
                contributions.append((ax, s, face_rate))
        if smooth_length > 0.0:
            mask = (rate > 0.0) if smooth_mask is None else np.asarray(smooth_mask, bool)
            # Even adjacent growth candidates can be fed by DIFFERENT liquid
            # clusters. Smooth each source cluster separately before summing.
            from .morphology import label_clusters
            # The kernel itself stays the Euclidean Gaussian of regularize.py (it must
            # smooth across the convex corners of a tip; a graph kernel on the
            # face-connected candidate set cannot, see README E1 section 9.13), so
            # the separation of sources is done here: cluster c may only spread onto
            # cells that own a face to c.
            labels, n = label_clusters(self.salt_mask(), g)
            rate = np.zeros(g.shape)
            for c in range(1, n + 1):
                own = np.zeros(g.shape)
                adjacent = np.zeros(g.shape, dtype=bool)
                for ax, s, fr in contributions:
                    touches = g.shift(labels, ax, s, 0) == c
                    own += np.where(touches, fr, 0.0)
                    adjacent |= touches
                rate += smooth_on_mask(own, mask & adjacent, g, smooth_length)
        return rate

    def dt_occupancy(self, rate: np.ndarray, f_max: float = 0.5) -> float:
        """Growth step limit from the CELL occupancy gain (plan/17 P1a item 1):
        ``dt = f_max / max_i (Vdot_i / V_cell)``.  ``inf`` if nothing grows."""
        peak = float(np.max(rate)) if rate.size else 0.0
        if peak <= 0.0:
            return np.inf
        return f_max * self.grid.cell_volume / peak

    def advance(self, rate: np.ndarray, dt: float, tol: float = 0.0, carry=None) -> StepReport:
        """Apply a precomputed per-cell volume rate for ``dt``.  Use this when
        the caller needs the rate (and hence ``dt_occupancy``) before deciding
        the step; ``step`` is the convenience wrapper that does both.

        ``carry`` (window-shaped, modified in place) travels with the overflow volume, so a
        companion quantity such as the salt accumulated in each candidate cell stays
        proportional to the volume it belongs to."""
        g = self.grid_w
        cand = rate > 0.0
        if dt <= 0.0 or not cand.any():
            self.t += max(dt, 0.0)
            return StepReport(max(dt, 0.0), np.zeros(g.shape, bool), 0.0, 0.0, 0)
        f = self.f_w
        f[cand] += rate[cand] * dt / g.cell_volume
        injected = float(rate[cand].sum()) * dt
        self.volume_injected += injected
        self.t += dt
        bank0 = self.bank
        flipped, passes = self._flip_and_redistribute(tol, carry=carry)
        return StepReport(dt, flipped, injected, self.bank - bank0, passes)

    def _flip_and_redistribute(self, tol: float, max_passes: int = 20, carry=None):
        """Flip full cells and hand each one's overflow to its OWN gas
        neighbours (plan/17 P1a item 2).  Volume is conserved exactly: what a
        flipped cell gives away equals what its neighbours receive, plus the
        part banked when a flipped cell has no gas neighbour left.

        ``carry`` is an optional companion field on the same cells (the driver
        uses it for the salt accumulated in each candidate cell).  It travels
        with the volume: a cell that gives away ``overflow`` out of ``1 +
        overflow`` gives the same FRACTION of its carry, so the carry stays
        proportional to the volume and a voxel's porosity does not depend on
        how the occupancy scheme happened to split the last step."""
        g = self.grid_w
        st = self.state_w; f = self.f_w
        if carry is not None and carry.shape != g.shape:
            raise ValueError("carry must have the window shape")
        flipped_all = np.zeros(g.shape, dtype=bool)
        for p in range(1, max_passes + 1):
            full = (st == GAS) & (f >= 1.0 - tol)
            if not full.any():
                return flipped_all, p - 1
            overflow = np.where(full, f - 1.0, 0.0)   # may be < 0 if tol > 0
            carry_out = None
            if carry is not None:
                frac = np.zeros(g.shape)
                np.divide(overflow, f, out=frac, where=full & (f > 0.0))
                carry_out = np.where(full, carry * np.clip(frac, 0.0, 1.0), 0.0)
                carry -= carry_out
            f[full] = 0.0
            st[full] = SALT_WET
            flipped_all |= full
            if not np.any(overflow):
                continue
            gas = st == GAS
            # per flipped cell: how many of its own neighbours can receive
            n_nb = np.zeros(g.shape)
            nb_gas = {}
            for ax, s in g.face_directions():
                nb_gas[(ax, s)] = g.shift(gas, ax, s, False)
                n_nb += nb_gas[(ax, s)]
            has_nb = full & (n_nb > 0)
            self.bank += float(overflow[full & (n_nb == 0)].sum()) * g.cell_volume
            share = np.zeros(g.shape)
            np.divide(overflow, n_nb, out=share, where=has_nb)
            share = np.where(has_nb, share, 0.0)
            share_c = None
            if carry is not None:
                self.bank_carry += float(carry_out[full & (n_nb == 0)].sum())
                share_c = np.zeros(g.shape)
                np.divide(carry_out, n_nb, out=share_c, where=has_nb)
                share_c = np.where(has_nb, share_c, 0.0)
            for ax, s in g.face_directions():
                # cell i gives share[i] to its neighbour at i+s when that is gas;
                # the receiver at index j reads the giver at j-s, i.e. shift(., -s)
                given = np.where(nb_gas[(ax, s)], share, 0.0)
                f += g.shift(given, ax, -s, 0.0)
                if carry is not None:
                    carry += g.shift(np.where(nb_gas[(ax, s)], share_c, 0.0), ax, -s, 0.0)
        if ((st == GAS) & (f >= 1.0 - tol)).any():
            raise RuntimeError(
                f"geometry overflow did not settle in {max_passes} passes: "
                f"{int(((st == GAS) & (f >= 1.0 - tol)).sum())} cells still full. "
                "Reduce f_max / dt (plan/17 P1a item 3)")
        return flipped_all, max_passes

    # ------------------------------------------------------------------
    @property
    def volume_placed(self) -> float:
        g = self.grid
        return float(self.solid_mask_full().sum()) * g.cell_volume + float(self.f.sum()) * g.cell_volume

    def conservation_residual(self, volume_initial: float) -> float:
        """(placed + banked - initial - injected) / max(injected, tiny)."""
        inj = self.volume_injected
        return (self.volume_placed + self.bank - volume_initial - inj) / max(inj, 1e-300)
