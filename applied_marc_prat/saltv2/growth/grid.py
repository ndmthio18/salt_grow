"""Lattice and voxel states, dimension-agnostic (2-D or 3-D).

Conventions
-----------
* Arrays have shape ``(n_1, ..., n_lat, n_z)``: the LAST axis is vertical
  (index 0 = the row sitting on the substrate surface z = 0), all other axes
  are lateral and periodic.
* Voxel size ``a`` [m] is uniform. Cell centre of index k along z is
  ``(k + 0.5) a``.
* States are small integers so that masks are cheap.
"""
from __future__ import annotations

import numpy as np

GAS = 0          # air (vapour diffuses)
SALT_WET = 1     # porous salt with mobile brine; inventory owns its saturation
SALT_SEALED = 4  # reaction band or bulk densified to eps_min: no liquid supply, no growth;
                 # conducts vapour through its closed band (crust.sealed_conductance).
                 # The value 4 is kept: saved records and report_3d.py rely on it (the
                 # intermediate "dry open" state 2 was removed 2026-09-14)


class Grid:
    """Uniform lattice with periodic lateral axes and a vertical last axis."""

    def __init__(self, shape, a: float):
        self.shape = tuple(int(n) for n in shape)
        self.ndim = len(self.shape)
        if self.ndim not in (2, 3):
            raise ValueError("Grid supports 2-D or 3-D lattices")
        self.a = float(a)
        self.z_axis = self.ndim - 1
        self.lateral_axes = tuple(range(self.ndim - 1))
        self.n_z = self.shape[-1]
        self.n_lateral = int(np.prod(self.shape[:-1]))
        # face area and cell volume per unit depth (2-D) or true (3-D)
        self.face_area = self.a ** (self.ndim - 1)
        self.cell_volume = self.a ** self.ndim

    # ------------------------------------------------------------------
    def z_centres(self) -> np.ndarray:
        """Cell-centre heights [m] along the vertical axis."""
        return (np.arange(self.n_z) + 0.5) * self.a

    def shift(self, arr: np.ndarray, axis: int, step: int, fill):
        """Neighbour lookup: value of the cell ``step`` cells further along
        ``axis``. Lateral axes wrap (periodic); the vertical axis is padded
        with ``fill`` beyond the domain."""
        if axis in self.lateral_axes:
            return np.roll(arr, -step, axis=axis)
        out = np.full_like(arr, fill)
        n = arr.shape[axis]
        src = [slice(None)] * arr.ndim
        dst = [slice(None)] * arr.ndim
        if step > 0:
            src[axis] = slice(step, n)
            dst[axis] = slice(0, n - step)
        else:
            src[axis] = slice(0, n + step)
            dst[axis] = slice(-step, n)
        out[tuple(dst)] = arr[tuple(src)]
        return out

    def face_directions(self):
        """All (axis, step) pairs: 2*ndim face directions."""
        return [(ax, s) for ax in range(self.ndim) for s in (+1, -1)]

    def column_max_height(self, salt_mask: np.ndarray) -> np.ndarray:
        """Height [m] of the topmost salt voxel in each lateral column
        (0 where the column holds no salt)."""
        k = np.arange(1, self.n_z + 1)
        shape = [1] * self.ndim
        shape[-1] = self.n_z
        top = np.max(np.where(salt_mask, k.reshape(shape), 0), axis=self.z_axis)
        return top * self.a
