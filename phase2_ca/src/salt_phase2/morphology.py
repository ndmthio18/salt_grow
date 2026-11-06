"""Connectivity, envelope, and geometric diagnostics for Phase 2A."""

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
    """Void diagnostics inside one explicitly named analysis mask."""

    name: str
    envelope_voxels: int
    salt_voxels: int
    void_voxels: int
    salt_volume_fraction: float | None
    void_volume_fraction: float | None
    void_component_count: int
    void_component_sizes: tuple[int, ...]
    external_air_connected_void_voxels: int
    enclosed_void_voxels: int
    enclosed_void_component_count: int
    enclosed_void_component_sizes: tuple[int, ...]
    bottom_connected_void_voxels: int
    bottom_to_top_void_percolates: bool
    bottom_to_top_pathway_voxels: int
    bottom_to_lateral_void_percolates: bool
    bottom_to_lateral_pathway_voxels: int
    accessible_pore_volume_m3: float | None
    enclosed_void_volume_m3: float | None


@dataclass(frozen=True)
class MorphologyAnalysis:
    salt_cluster_count: int
    salt_cluster_sizes: tuple[int, ...]
    seed_occupied_cluster_count: int
    external_air_voxels: int
    supply_connected_salt_voxels: int
    computational_box: EnvelopeVoidStatistics
    bounding_box: EnvelopeVoidStatistics
    column_envelope: EnvelopeVoidStatistics
    face_shell_envelope: EnvelopeVoidStatistics
    salt_void_surface_faces: int
    salt_surface_area_m2: float | None
    salt_volume_m3: float | None
    surface_to_volume_ratio_m_inv: float | None
    characteristic_solid_size_m: float | None
    tested_envelope_void_fraction_range: float
    tested_envelope_pathway_agreement: bool
    intrinsic_porosity_identifiable: bool
    intrinsic_porosity_flag: str


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
    """Return empty voxels connected to top or lateral domain boundaries."""

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


def computational_box_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    """Return the whole computational region above the prescribed substrate."""

    envelope = np.zeros(state.shape, dtype=bool)
    envelope[:, :, state.substrate_layers :] = True
    return envelope


def bounding_box_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    """Return the tight axis-aligned bounding box containing all salt."""

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


def column_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    """Fill each occupied x-y column from substrate surface to local salt top."""

    salt = state.material == Material.SALT
    envelope = np.zeros(state.shape, dtype=bool)
    for x, y in np.argwhere(np.any(salt, axis=2)):
        occupied_z = np.flatnonzero(salt[int(x), int(y), :])
        top = int(np.max(occupied_z))
        envelope[int(x), int(y), state.substrate_layers : top + 1] = True
    return envelope


def face_shell_envelope(state: PhysicalState) -> NDArray[np.bool_]:
    """Return salt plus its one-voxel six-face shell, clipped above substrate."""

    salt = state.material == Material.SALT
    shell = salt.copy()
    shell[1:, :, :] |= salt[:-1, :, :]
    shell[:-1, :, :] |= salt[1:, :, :]
    shell[:, 1:, :] |= salt[:, :-1, :]
    shell[:, :-1, :] |= salt[:, 1:, :]
    shell[:, :, 1:] |= salt[:, :, :-1]
    shell[:, :, :-1] |= salt[:, :, 1:]
    shell[:, :, : state.substrate_layers] = False
    return shell


def _boundary_masks(envelope: NDArray[np.bool_]) -> tuple[NDArray[np.bool_], ...]:
    bottom = np.zeros(envelope.shape, dtype=bool)
    top = np.zeros(envelope.shape, dtype=bool)
    lateral = np.zeros(envelope.shape, dtype=bool)
    coordinates = np.argwhere(envelope)
    if coordinates.size == 0:
        return bottom, top, lateral
    lower = np.min(coordinates, axis=0)
    upper = np.max(coordinates, axis=0)
    bottom[:, :, int(lower[2])] = envelope[:, :, int(lower[2])]
    top[:, :, int(upper[2])] = envelope[:, :, int(upper[2])]
    lateral[int(lower[0]), :, :] |= envelope[int(lower[0]), :, :]
    lateral[int(upper[0]), :, :] |= envelope[int(upper[0]), :, :]
    lateral[:, int(lower[1]), :] |= envelope[:, int(lower[1]), :]
    lateral[:, int(upper[1]), :] |= envelope[:, int(upper[1]), :]
    return bottom, top, lateral


def _labels_touching(
    labels: NDArray[np.int32], boundary: NDArray[np.bool_]
) -> NDArray[np.int32]:
    values = np.unique(labels[boundary & (labels > 0)])
    return values[values > 0]


def _envelope_statistics(
    *,
    name: str,
    envelope: NDArray[np.bool_],
    salt: NDArray[np.bool_],
    external_air: NDArray[np.bool_],
    voxel_size_m: float | None,
) -> EnvelopeVoidStatistics:
    void = envelope & ~salt
    labels, void_sizes = connected_components(void)
    enclosed = void & ~external_air
    _, enclosed_sizes = connected_components(enclosed)
    bottom, top, lateral = _boundary_masks(envelope)
    bottom_labels = _labels_touching(labels, bottom)
    top_labels = _labels_touching(labels, top)
    lateral_labels = _labels_touching(labels, lateral)
    bottom_top_labels = np.intersect1d(bottom_labels, top_labels)
    bottom_lateral_labels = np.intersect1d(bottom_labels, lateral_labels)
    bottom_connected = np.isin(labels, bottom_labels)
    bottom_top = np.isin(labels, bottom_top_labels)
    bottom_lateral = np.isin(labels, bottom_lateral_labels)
    envelope_count = int(np.count_nonzero(envelope))
    salt_count = int(np.count_nonzero(envelope & salt))
    void_count = int(np.count_nonzero(void))
    accessible_count = int(np.count_nonzero(void & external_air))
    enclosed_count = int(np.count_nonzero(enclosed))
    voxel_volume = None if voxel_size_m is None else voxel_size_m**3
    return EnvelopeVoidStatistics(
        name=name,
        envelope_voxels=envelope_count,
        salt_voxels=salt_count,
        void_voxels=void_count,
        salt_volume_fraction=(salt_count / envelope_count if envelope_count else None),
        void_volume_fraction=(void_count / envelope_count if envelope_count else None),
        void_component_count=len(void_sizes),
        void_component_sizes=void_sizes,
        external_air_connected_void_voxels=accessible_count,
        enclosed_void_voxels=enclosed_count,
        enclosed_void_component_count=len(enclosed_sizes),
        enclosed_void_component_sizes=enclosed_sizes,
        bottom_connected_void_voxels=int(np.count_nonzero(bottom_connected)),
        bottom_to_top_void_percolates=bool(bottom_top_labels.size),
        bottom_to_top_pathway_voxels=int(np.count_nonzero(bottom_top)),
        bottom_to_lateral_void_percolates=bool(bottom_lateral_labels.size),
        bottom_to_lateral_pathway_voxels=int(np.count_nonzero(bottom_lateral)),
        accessible_pore_volume_m3=(
            None if voxel_volume is None else accessible_count * voxel_volume
        ),
        enclosed_void_volume_m3=(
            None if voxel_volume is None else enclosed_count * voxel_volume
        ),
    )


def _salt_void_surface_faces(state: PhysicalState) -> int:
    salt = state.material == Material.SALT
    empty = state.material == Material.EMPTY
    faces = 0
    faces += int(np.count_nonzero(salt[1:, :, :] & empty[:-1, :, :]))
    faces += int(np.count_nonzero(salt[:-1, :, :] & empty[1:, :, :]))
    faces += int(np.count_nonzero(salt[:, 1:, :] & empty[:, :-1, :]))
    faces += int(np.count_nonzero(salt[:, :-1, :] & empty[:, 1:, :]))
    faces += int(np.count_nonzero(salt[:, :, 1:] & empty[:, :, :-1]))
    faces += int(np.count_nonzero(salt[:, :, :-1] & empty[:, :, 1:]))
    return faces


def analyze_morphology(
    state: PhysicalState,
    *,
    supply_connected_salt: NDArray[np.bool_] | None = None,
    voxel_size_m: float | None = None,
) -> MorphologyAnalysis:
    """Analyze solid and void connectivity under explicit, non-equivalent masks."""

    state.validate()
    if voxel_size_m is not None and (
        not np.isfinite(voxel_size_m) or voxel_size_m <= 0
    ):
        raise ValueError("voxel_size_m must be finite and positive")
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

    masks = {
        "computational_box_above_substrate": computational_box_envelope(state),
        "axis_aligned_salt_bounding_box": bounding_box_envelope(state),
        "footprint_columns_to_local_salt_top": column_envelope(state),
        "salt_plus_one_voxel_face_shell": face_shell_envelope(state),
    }
    statistics = {
        name: _envelope_statistics(
            name=name,
            envelope=envelope,
            salt=salt,
            external_air=external_air,
            voxel_size_m=voxel_size_m,
        )
        for name, envelope in masks.items()
    }
    tested_names = tuple(
        name for name in masks if name != "computational_box_above_substrate"
    )
    tested = [statistics[name] for name in tested_names]
    fractions = [
        item.void_volume_fraction
        for item in tested
        if item.void_volume_fraction is not None
    ]
    fraction_range = float(max(fractions) - min(fractions)) if fractions else 0.0
    signatures = {
        (
            item.bottom_to_top_void_percolates,
            item.bottom_to_lateral_void_percolates,
        )
        for item in tested
    }
    pathway_agreement = len(signatures) <= 1
    masks_equal = all(
        np.array_equal(masks[tested_names[0]], masks[name])
        for name in tested_names[1:]
    )
    identifiable = bool(masks_equal)
    flag = (
        "tested envelope masks coincide"
        if identifiable
        else "tested envelope masks differ; no unique intrinsic pore volume"
    )

    surface_faces = _salt_void_surface_faces(state)
    salt_count = int(np.count_nonzero(salt))
    if voxel_size_m is None:
        surface_area = volume = ratio = characteristic_size = None
    else:
        surface_area = surface_faces * voxel_size_m**2
        volume = salt_count * voxel_size_m**3
        ratio = surface_area / volume if volume else None
        characteristic_size = 6.0 * volume / surface_area if surface_area else None

    return MorphologyAnalysis(
        salt_cluster_count=len(salt_sizes),
        salt_cluster_sizes=salt_sizes,
        seed_occupied_cluster_count=len(seed_labels),
        external_air_voxels=int(np.count_nonzero(external_air)),
        supply_connected_salt_voxels=supply_connected_count,
        computational_box=statistics["computational_box_above_substrate"],
        bounding_box=statistics["axis_aligned_salt_bounding_box"],
        column_envelope=statistics["footprint_columns_to_local_salt_top"],
        face_shell_envelope=statistics["salt_plus_one_voxel_face_shell"],
        salt_void_surface_faces=surface_faces,
        salt_surface_area_m2=surface_area,
        salt_volume_m3=volume,
        surface_to_volume_ratio_m_inv=ratio,
        characteristic_solid_size_m=characteristic_size,
        tested_envelope_void_fraction_range=fraction_range,
        tested_envelope_pathway_agreement=pathway_agreement,
        intrinsic_porosity_identifiable=identifiable,
        intrinsic_porosity_flag=flag,
    )
