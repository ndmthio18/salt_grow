"""Connectivity and void-morphology diagnostics for Phase 2A.1."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .state import Coordinate, Material, PhysicalState


_FACE_OFFSETS: tuple[Coordinate, ...] = (
    (-1, 0, 0),
    (1, 0, 0),
    (0, -1, 0),
    (0, 1, 0),
    (0, 0, -1),
    (0, 0, 1),
)


@dataclass(frozen=True)
class EnvelopeVoidStatistics:
    name: str
    envelope_voxels: int
    salt_voxels: int
    void_voxels: int
    salt_volume_fraction: float | None
    void_component_count: int
    void_component_sizes: tuple[int, ...]
    external_air_connected_void_voxels: int
    enclosed_void_voxels: int
    enclosed_void_component_count: int
    enclosed_void_component_sizes: tuple[int, ...]


@dataclass(frozen=True)
class MorphologyAnalysis:
    salt_cluster_count: int
    salt_cluster_sizes: tuple[int, ...]
    seed_occupied_cluster_count: int
    external_air_voxels: int
    supply_connected_salt_voxels: int
    bounding_box: EnvelopeVoidStatistics
    column_envelope: EnvelopeVoidStatistics


def connected_components(
    mask: NDArray[np.bool_],
) -> tuple[NDArray[np.int32], tuple[int, ...]]:
    """Label non-periodic six-face components and return descending sizes."""

    if mask.ndim != 3:
        raise ValueError("component mask must be three-dimensional")
    labels = np.zeros(mask.shape, dtype=np.int32)
    sizes: list[int] = []
    shape = mask.shape
    next_label = 0
    for start_array in np.argwhere(mask & (labels == 0)):
        start = tuple(int(value) for value in start_array)
        if labels[start] != 0:
            continue
        next_label += 1
        labels[start] = next_label
        queue: deque[Coordinate] = deque([start])
        size = 0
        while queue:
            x, y, z = queue.popleft()
            size += 1
            for dx, dy, dz in _FACE_OFFSETS:
                neighbor = x + dx, y + dy, z + dz
                if not (
                    0 <= neighbor[0] < shape[0]
                    and 0 <= neighbor[1] < shape[1]
                    and 0 <= neighbor[2] < shape[2]
                ):
                    continue
                if mask[neighbor] and labels[neighbor] == 0:
                    labels[neighbor] = next_label
                    queue.append(neighbor)
        sizes.append(size)
    return labels, tuple(sorted(sizes, reverse=True))


def external_air_connected_mask(state: PhysicalState) -> NDArray[np.bool_]:
    """Return empty voxels connected to the top or lateral domain boundaries."""

    empty = state.material == Material.EMPTY
    boundary = np.zeros(state.shape, dtype=bool)
    boundary[0, :, :] = True
    boundary[-1, :, :] = True
    boundary[:, 0, :] = True
    boundary[:, -1, :] = True
    boundary[:, :, -1] = True
    labels, _ = connected_components(empty)
    external_labels = np.unique(labels[boundary & empty])
    external_labels = external_labels[external_labels > 0]
    return np.isin(labels, external_labels)


def _bounding_box_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    salt_coordinates = np.argwhere(state.material == Material.SALT)
    envelope = np.zeros(state.shape, dtype=bool)
    if salt_coordinates.size == 0:
        return envelope
    lower = np.min(salt_coordinates, axis=0)
    upper = np.max(salt_coordinates, axis=0)
    lower[2] = max(int(lower[2]), state.substrate_layers)
    envelope[
        int(lower[0]) : int(upper[0]) + 1,
        int(lower[1]) : int(upper[1]) + 1,
        int(lower[2]) : int(upper[2]) + 1,
    ] = True
    return envelope


def _column_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    salt = state.material == Material.SALT
    envelope = np.zeros(state.shape, dtype=bool)
    for x, y in np.argwhere(np.any(salt, axis=2)):
        occupied_z = np.flatnonzero(salt[int(x), int(y), :])
        top = int(np.max(occupied_z))
        envelope[int(x), int(y), state.substrate_layers : top + 1] = True
    return envelope


def _envelope_statistics(
    *,
    name: str,
    envelope: NDArray[np.bool_],
    salt: NDArray[np.bool_],
    external_air: NDArray[np.bool_],
) -> EnvelopeVoidStatistics:
    void = envelope & ~salt
    _, void_sizes = connected_components(void)
    enclosed = void & ~external_air
    _, enclosed_sizes = connected_components(enclosed)
    envelope_count = int(np.count_nonzero(envelope))
    salt_count = int(np.count_nonzero(envelope & salt))
    void_count = int(np.count_nonzero(void))
    return EnvelopeVoidStatistics(
        name=name,
        envelope_voxels=envelope_count,
        salt_voxels=salt_count,
        void_voxels=void_count,
        salt_volume_fraction=(salt_count / envelope_count if envelope_count else None),
        void_component_count=len(void_sizes),
        void_component_sizes=void_sizes,
        external_air_connected_void_voxels=int(np.count_nonzero(void & external_air)),
        enclosed_void_voxels=int(np.count_nonzero(enclosed)),
        enclosed_void_component_count=len(enclosed_sizes),
        enclosed_void_component_sizes=enclosed_sizes,
    )


def analyze_morphology(
    state: PhysicalState,
    *,
    supply_connected_salt: NDArray[np.bool_] | None = None,
) -> MorphologyAnalysis:
    """Analyze solid clusters and voids under two explicit envelope assumptions."""

    state.validate()
    salt = state.material == Material.SALT
    salt_labels, salt_sizes = connected_components(salt)
    seed_labels = np.unique(salt_labels[state.seed_mask])
    seed_labels = seed_labels[seed_labels > 0]
    external_air = external_air_connected_mask(state)
    if supply_connected_salt is None:
        supply_connected_count = 0
    else:
        if supply_connected_salt.shape != state.shape:
            raise ValueError("supply_connected_salt must match state shape")
        supply_connected_count = int(np.count_nonzero(supply_connected_salt & salt))
    return MorphologyAnalysis(
        salt_cluster_count=len(salt_sizes),
        salt_cluster_sizes=salt_sizes,
        seed_occupied_cluster_count=len(seed_labels),
        external_air_voxels=int(np.count_nonzero(external_air)),
        supply_connected_salt_voxels=supply_connected_count,
        bounding_box=_envelope_statistics(
            name="axis_aligned_salt_bounding_box",
            envelope=_bounding_box_envelope(state),
            salt=salt,
            external_air=external_air,
        ),
        column_envelope=_envelope_statistics(
            name="footprint_columns_to_local_salt_top",
            envelope=_column_envelope(state),
            salt=salt,
            external_air=external_air,
        ),
    )
