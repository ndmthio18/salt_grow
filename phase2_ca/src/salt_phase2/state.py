"""Physical state representation for the minimal Phase 2A salt-growth CA."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Sequence

import numpy as np
from numpy.typing import NDArray


class Material(IntEnum):
    """Material identity stored in the CA lattice."""

    EMPTY = 0
    SUBSTRATE = 1
    SALT = 2


Coordinate = tuple[int, int, int]


def _validate_shape(shape: Sequence[int]) -> tuple[int, int, int]:
    if len(shape) != 3:
        raise ValueError("shape must contain exactly three dimensions")
    checked: list[int] = []
    for value in shape:
        if isinstance(value, bool) or int(value) != value or value <= 0:
            raise ValueError("shape dimensions must be positive integers")
        checked.append(int(value))
    return checked[0], checked[1], checked[2]


def _validate_coordinate(coordinate: Sequence[int], shape: tuple[int, int, int]) -> Coordinate:
    if len(coordinate) != 3:
        raise ValueError("each seed coordinate must have three entries")
    checked: list[int] = []
    for value, upper in zip(coordinate, shape, strict=True):
        if isinstance(value, bool) or int(value) != value:
            raise ValueError("seed coordinates must be integers")
        integer = int(value)
        if not 0 <= integer < upper:
            raise ValueError("seed coordinate lies outside the domain")
        checked.append(integer)
    return checked[0], checked[1], checked[2]


@dataclass
class PhysicalState:
    """Independent material, liquid-accessibility, and seed arrays.

    ``liquid_accessible`` is not an occupancy flag.  On a salt voxel it marks a
    direct liquid-film supply anchor.  Film access on other salt voxels is a
    derived connectivity property, recomputed from the solid salt network.
    Empty target voxels may therefore be air-filled and still be eligible when
    they touch supply-connected salt.
    """

    material: NDArray[np.uint8]
    liquid_accessible: NDArray[np.bool_]
    seed_mask: NDArray[np.bool_]
    substrate_layers: int

    def validate(self) -> None:
        if self.material.ndim != 3:
            raise ValueError("material must be a three-dimensional array")
        if self.material.dtype != np.uint8:
            raise ValueError("material must use uint8 storage")
        if self.liquid_accessible.dtype != np.bool_:
            raise ValueError("liquid_accessible must use boolean storage")
        if self.seed_mask.dtype != np.bool_:
            raise ValueError("seed_mask must use boolean storage")
        if self.liquid_accessible.shape != self.material.shape:
            raise ValueError("liquid_accessible must match material shape")
        if self.seed_mask.shape != self.material.shape:
            raise ValueError("seed_mask must match material shape")
        allowed = np.isin(
            self.material,
            np.array([Material.EMPTY, Material.SUBSTRATE, Material.SALT], dtype=np.uint8),
        )
        if not bool(np.all(allowed)):
            raise ValueError("material contains an unknown state")
        if np.any(self.seed_mask & (self.material != Material.SALT)):
            raise ValueError("every seed must remain a salt voxel")
        if not 1 <= self.substrate_layers < self.material.shape[2]:
            raise ValueError("substrate_layers must leave space above the substrate")
        if not np.all(self.material[:, :, : self.substrate_layers] == Material.SUBSTRATE):
            raise ValueError("the prescribed bottom substrate is incomplete")

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.material.shape

    @property
    def initial_seed_voxel_count(self) -> int:
        return int(np.count_nonzero(self.seed_mask))


def initialize_physical_state(
    *,
    shape: Sequence[int],
    substrate_layers: int,
    seed_coordinates: Sequence[Sequence[int]],
    wetted_seed_indices: Sequence[int] | None = None,
) -> PhysicalState:
    """Create a bottom substrate and salt seeds immediately above its surface."""

    checked_shape = _validate_shape(shape)
    if isinstance(substrate_layers, bool) or int(substrate_layers) != substrate_layers:
        raise ValueError("substrate_layers must be an integer")
    substrate_layers = int(substrate_layers)
    if not 1 <= substrate_layers < checked_shape[2]:
        raise ValueError("substrate_layers must leave space above the substrate")
    if len(seed_coordinates) == 0:
        raise ValueError("at least one seed coordinate is required")

    coordinates = [
        _validate_coordinate(coordinate, checked_shape) for coordinate in seed_coordinates
    ]
    if len(set(coordinates)) != len(coordinates):
        raise ValueError("seed coordinates must be unique")
    if any(coordinate[2] != substrate_layers for coordinate in coordinates):
        raise ValueError("seeds must lie immediately above the substrate interface")

    if wetted_seed_indices is None:
        wet_indices = set(range(len(coordinates)))
    else:
        wet_indices: set[int] = set()
        for index in wetted_seed_indices:
            if isinstance(index, bool) or int(index) != index:
                raise ValueError("wetted_seed_indices must contain integers")
            integer = int(index)
            if not 0 <= integer < len(coordinates):
                raise ValueError("wetted seed index is out of range")
            wet_indices.add(integer)

    material = np.full(checked_shape, Material.EMPTY, dtype=np.uint8)
    material[:, :, :substrate_layers] = Material.SUBSTRATE
    seed_mask = np.zeros(checked_shape, dtype=bool)
    liquid_accessible = np.zeros(checked_shape, dtype=bool)
    for index, coordinate in enumerate(coordinates):
        material[coordinate] = Material.SALT
        seed_mask[coordinate] = True
        if index in wet_indices:
            liquid_accessible[coordinate] = True

    state = PhysicalState(
        material=material,
        liquid_accessible=liquid_accessible,
        seed_mask=seed_mask,
        substrate_layers=substrate_layers,
    )
    state.validate()
    return state
