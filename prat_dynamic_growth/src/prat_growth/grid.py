from __future__ import annotations

import numpy as np


def top_refined_grid(n_cells: int, power: float = 2.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return normalized ALE faces, centers, and widths on [0, 1]."""
    if n_cells < 4:
        raise ValueError("n_cells must be at least 4")
    if power < 1.0:
        raise ValueError("grid power must be at least 1")
    s = np.linspace(0.0, 1.0, n_cells + 1)
    faces = 1.0 - (1.0 - s) ** power
    centers = 0.5 * (faces[:-1] + faces[1:])
    widths = np.diff(faces)
    return faces, centers, widths


def cells_across_top_layer(faces: np.ndarray, height: float, layer_thickness: float) -> int:
    """Count complete/partial cells intersecting a layer measured from the top."""
    if layer_thickness <= 0.0:
        return 0
    physical_depths = height * (1.0 - faces)
    # A cell intersects the top layer when its top/right face lies within it.
    return int(np.count_nonzero(physical_depths[1:] < layer_thickness))
