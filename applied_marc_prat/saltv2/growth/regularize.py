"""Regularisation length for the lattice growth (plan/17 section P1b).

Why this exists
---------------
Pure Laplacian growth has no fastest wavelength: every tip is unstable down
to the smallest resolved scale, so the selected morphology is set by whatever
lateral smoothing the discretisation happens to provide.  In
``GeometryEngine`` that smoothing was the overflow redistribution, whose
magnitude scales with ``dt``; the 3-D ``dt`` scan of 2026-09 showed the columns narrowing monotonically as ``dt`` decreases,
i.e. no convergence.  Removing the ``dt``-dependent smoothing (P1a) leaves the
voxel size ``a`` as the only cut-off, so an explicit physical length must be
put into the growth law instead.

The regularisation
------------------
The solid-volume deposition rate is redistributed over a length ``ell`` along
the growth front before it drives the occupancy.  In Fourier space this
multiplies the growth rate of a perturbation of wavenumber ``k`` by
``exp(-k^2 ell^2 / 2)``, which turns the Laplacian instability into one with a
fastest-growing wavelength of order ``ell``.

Physical reading (**[A], assumption**): ``ell`` is the lateral redistribution
length of the liquid inside the top layer of the salt structure.  Candidate
values are the first layer thickness ``h_1 = c_h d_sub`` or the first pore
scale ``d_1 = 0.14 d_sub`` of the plan/13 stack.  ``ell`` must NOT be tuned on
the 2014 evaporation plateau or on the crusty/patchy transition position
(plan/17 P1b item 4); it is reported as a sensitivity.

Numerics
--------
``smooth_on_mask`` is a normalised (Nadaraya-Watson) Gaussian spread, written
so that the total over the mask is conserved exactly rather than only in the
interior:

    norm_i = sum_k K(i-k) mask_k                (kernel weight seen by source i)
    out_j  = mask_j * sum_i K(j-i) field_i / norm_i

Summing ``out`` over ``j`` gives back ``sum_i field_i`` for any mask and any
boundary treatment that extends the arrays consistently, because the same
truncated kernel enters ``norm_i`` and the sum over ``j``.  Lateral axes use
``wrap`` (the lattice is periodic); the vertical axis uses ``constant`` with
``cval = 0`` so that nothing is created or destroyed at the box top/bottom.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from .grid import Grid


def smooth_on_mask(field: np.ndarray, mask: np.ndarray, grid: Grid,
                   length: float) -> np.ndarray:
    """Spread ``field`` over the physical length ``length`` [m], staying inside
    ``mask`` and conserving ``field[mask].sum()``.

    ``length <= 0`` returns ``field`` restricted to the mask unchanged (no
    regularisation), so the caller can switch the length off without a branch.
    """
    field = np.asarray(field, float)
    mask = np.asarray(mask, bool)
    if field.shape != grid.shape or mask.shape != grid.shape:
        raise ValueError("field and mask must have the shape of the grid")
    inside = np.where(mask, field, 0.0)
    if not (length > 0.0) or not mask.any():
        return inside
    sigma = float(length) / grid.a
    truncate = 4.0                       # kernel support in sigmas (scipy default)
    modes = ["wrap"] * len(grid.lateral_axes) + ["constant"]
    norm = ndimage.gaussian_filter(mask.astype(float), sigma, mode=modes,
                                   cval=0.0, truncate=truncate)
    # norm >= K(0) > 0 wherever mask is True; the division is only used there.
    src = np.zeros(grid.shape)
    np.divide(inside, norm, out=src, where=mask)
    out = ndimage.gaussian_filter(src, sigma, mode=modes, cval=0.0,
                                  truncate=truncate)
    return np.where(mask, out, 0.0)
