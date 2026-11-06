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
    boundary_definition: str


@dataclass(frozen=True)
class PhysicalBoundarySurfaces:
    """Explicit inlet and outlet surfaces independent of an analysis envelope."""

    name: str
    inlet: NDArray[np.bool_]
    top_outlet: NDArray[np.bool_]
    lateral_outlet: NDArray[np.bool_]

    def validate(self, shape: tuple[int, int, int]) -> None:
        for name, mask in (
            ("inlet", self.inlet),
            ("top_outlet", self.top_outlet),
            ("lateral_outlet", self.lateral_outlet),
        ):
            if mask.shape != shape or mask.dtype != np.bool_:
                raise ValueError(f"{name} must be a boolean mask matching the domain")
        if np.any(self.inlet & self.top_outlet):
            raise ValueError("inlet and top_outlet must not overlap")


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
    fixed_window_solid_volume_fraction: float | None
    intrinsic_porosity: float | None
    accessible_intergranular_void_voxels: int | None
    accessible_intergranular_void_volume_m3: float | None
    ambiguous_exterior_void_voxels: int
    ambiguous_exterior_void_volume_m3: float | None
    physical_shell_thickness_m: float | None


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


def rectangular_window_mask(
    shape: tuple[int, int, int],
    *,
    lower: Coordinate,
    upper_exclusive: Coordinate,
) -> NDArray[np.bool_]:
    """Return a fixed rectangular analysis window with explicit coordinates."""

    if any(
        lower[axis] < 0
        or upper_exclusive[axis] > shape[axis]
        or lower[axis] >= upper_exclusive[axis]
        for axis in range(3)
    ):
        raise ValueError("window bounds must define a non-empty box inside the domain")
    mask = np.zeros(shape, dtype=bool)
    mask[
        lower[0] : upper_exclusive[0],
        lower[1] : upper_exclusive[1],
        lower[2] : upper_exclusive[2],
    ] = True
    return mask


def rectangular_boundary_surfaces(
    shape: tuple[int, int, int],
    *,
    lower: Coordinate,
    upper_exclusive: Coordinate,
    name: str,
) -> PhysicalBoundarySurfaces:
    """Create fixed planar inlet/top/lateral surfaces for a physical window."""

    window = rectangular_window_mask(
        shape, lower=lower, upper_exclusive=upper_exclusive
    )
    inlet = np.zeros(shape, dtype=bool)
    top = np.zeros(shape, dtype=bool)
    lateral = np.zeros(shape, dtype=bool)
    inlet[
        lower[0] : upper_exclusive[0],
        lower[1] : upper_exclusive[1],
        lower[2],
    ] = True
    top[
        lower[0] : upper_exclusive[0],
        lower[1] : upper_exclusive[1],
        upper_exclusive[2] - 1,
    ] = True
    lateral[lower[0], :, :] |= window[lower[0], :, :]
    lateral[upper_exclusive[0] - 1, :, :] |= window[upper_exclusive[0] - 1, :, :]
    lateral[:, lower[1], :] |= window[:, lower[1], :]
    lateral[:, upper_exclusive[1] - 1, :] |= window[:, upper_exclusive[1] - 1, :]
    boundaries = PhysicalBoundarySurfaces(
        name=name,
        inlet=inlet,
        top_outlet=top,
        lateral_outlet=lateral,
    )
    boundaries.validate(shape)
    return boundaries


def domain_boundary_surfaces(state: PhysicalState) -> PhysicalBoundarySurfaces:
    """Return substrate, domain-top, and domain-lateral physical surfaces."""

    return rectangular_boundary_surfaces(
        state.shape,
        lower=(0, 0, state.substrate_layers),
        upper_exclusive=state.shape,
        name="computational_domain_surfaces",
    )


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


def face_shell_envelope(
    state: PhysicalState,
    *,
    voxel_size_m: float | None = None,
    shell_thickness_m: float | None = None,
) -> NDArray[np.bool_]:
    """Return salt plus a six-face shell of explicit voxel or physical width."""

    if (voxel_size_m is None) != (shell_thickness_m is None):
        raise ValueError(
            "voxel_size_m and shell_thickness_m must be supplied together"
        )
    if voxel_size_m is None:
        shell_steps = 1
    else:
        if (
            not np.isfinite(voxel_size_m)
            or voxel_size_m <= 0.0
            or not np.isfinite(shell_thickness_m)
            or shell_thickness_m <= 0.0
        ):
            raise ValueError("voxel size and shell thickness must be positive")
        ratio = float(shell_thickness_m / voxel_size_m)
        shell_steps = int(round(ratio))
        if shell_steps < 1 or not np.isclose(
            ratio, shell_steps, rtol=0.0, atol=1.0e-10
        ):
            raise ValueError("shell thickness must be an integer number of voxels")

    salt = state.material == Material.SALT
    shell = salt.copy()
    for _ in range(shell_steps):
        previous = shell.copy()
        shell[1:, :, :] |= previous[:-1, :, :]
        shell[:-1, :, :] |= previous[1:, :, :]
        shell[:, 1:, :] |= previous[:, :-1, :]
        shell[:, :-1, :] |= previous[:, 1:, :]
        shell[:, :, 1:] |= previous[:, :, :-1]
        shell[:, :, :-1] |= previous[:, :, 1:]
    shell[:, :, : state.substrate_layers] = False
    return shell


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
    boundaries: PhysicalBoundarySurfaces,
    voxel_size_m: float | None,
) -> EnvelopeVoidStatistics:
    boundaries.validate(envelope.shape)
    void = envelope & ~salt
    labels, void_sizes = connected_components(void)
    enclosed = void & ~external_air
    _, enclosed_sizes = connected_components(enclosed)
    bottom_labels = _labels_touching(labels, boundaries.inlet)
    top_labels = _labels_touching(labels, boundaries.top_outlet)
    lateral_labels = _labels_touching(labels, boundaries.lateral_outlet)
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
        boundary_definition=boundaries.name,
    )


def analyze_void_space(
    *,
    name: str,
    analysis_mask: NDArray[np.bool_],
    salt_mask: NDArray[np.bool_],
    external_air_mask: NDArray[np.bool_],
    boundaries: PhysicalBoundarySurfaces,
    voxel_size_m: float | None = None,
) -> EnvelopeVoidStatistics:
    """Analyze void paths using caller-defined material and boundary masks."""

    if (
        analysis_mask.shape != salt_mask.shape
        or analysis_mask.shape != external_air_mask.shape
    ):
        raise ValueError("analysis, salt, and external-air masks must match")
    for mask in (analysis_mask, salt_mask, external_air_mask):
        if mask.dtype != np.bool_:
            raise ValueError("analysis masks must be boolean")
    return _envelope_statistics(
        name=name,
        envelope=analysis_mask,
        salt=salt_mask,
        external_air=external_air_mask,
        boundaries=boundaries,
        voxel_size_m=voxel_size_m,
    )


def face_adjacent_intergranular_gap_mask(
    state: PhysicalState,
    grain_ids: NDArray[np.integer],
    *,
    analysis_mask: NDArray[np.bool_] | None = None,
) -> NDArray[np.bool_]:
    """Return empty cells touching at least two distinct positive grain IDs."""

    if grain_ids.shape != state.shape or not np.issubdtype(grain_ids.dtype, np.integer):
        raise ValueError("grain_ids must be an integer array matching state shape")
    if np.any((state.material == Material.SALT) & (grain_ids <= 0)):
        raise ValueError("every salt voxel must have a positive grain ID")
    if np.any((state.material != Material.SALT) & (grain_ids != 0)):
        raise ValueError("non-salt voxels must have grain ID zero")
    if analysis_mask is None:
        analysis = computational_box_envelope(state)
    else:
        if analysis_mask.shape != state.shape or analysis_mask.dtype != np.bool_:
            raise ValueError("analysis_mask must be boolean and match state shape")
        analysis = analysis_mask
    candidate = analysis & (state.material == Material.EMPTY)
    gap = np.zeros(state.shape, dtype=bool)
    shape = state.shape
    for coordinate_array in np.argwhere(candidate):
        coordinate = tuple(int(value) for value in coordinate_array)
        x, y, z = coordinate
        neighbors: set[int] = set()
        for dx, dy, dz in _FACE_OFFSETS:
            neighbor = x + dx, y + dy, z + dz
            if (
                0 <= neighbor[0] < shape[0]
                and 0 <= neighbor[1] < shape[1]
                and 0 <= neighbor[2] < shape[2]
            ):
                grain_id = int(grain_ids[neighbor])
                if grain_id > 0:
                    neighbors.add(grain_id)
        if len(neighbors) >= 2:
            gap[coordinate] = True
    return gap


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
    boundary_surfaces: PhysicalBoundarySurfaces | None = None,
    fixed_window_mask: NDArray[np.bool_] | None = None,
    physical_shell_thickness_m: float | None = None,
    intrinsic_material_mask: NDArray[np.bool_] | None = None,
    intrinsic_boundary_source: str | None = None,
    grain_ids: NDArray[np.integer] | None = None,
) -> MorphologyAnalysis:
    """Analyze solid and void connectivity under explicit, non-equivalent masks."""

    state.validate()
    if voxel_size_m is not None and (
        not np.isfinite(voxel_size_m) or voxel_size_m <= 0
    ):
        raise ValueError("voxel_size_m must be finite and positive")
    if physical_shell_thickness_m is not None and voxel_size_m is None:
        raise ValueError("a physical shell thickness requires voxel_size_m")
    boundaries = boundary_surfaces or domain_boundary_surfaces(state)
    boundaries.validate(state.shape)
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

    computational_mask = computational_box_envelope(state)
    if fixed_window_mask is None:
        fixed_window = computational_mask
    else:
        if (
            fixed_window_mask.shape != state.shape
            or fixed_window_mask.dtype != np.bool_
        ):
            raise ValueError("fixed_window_mask must be boolean and match state shape")
        fixed_window = fixed_window_mask
    if physical_shell_thickness_m is None:
        shell = face_shell_envelope(state)
        shell_name = "salt_plus_one_voxel_face_shell"
        recorded_shell_thickness = voxel_size_m
    else:
        shell = face_shell_envelope(
            state,
            voxel_size_m=voxel_size_m,
            shell_thickness_m=physical_shell_thickness_m,
        )
        shell_name = (
            f"salt_plus_{physical_shell_thickness_m:.9g}_m_face_shell"
        )
        recorded_shell_thickness = physical_shell_thickness_m

    masks = {
        "computational_box_above_substrate": computational_mask,
        "axis_aligned_salt_bounding_box": bounding_box_envelope(state),
        "footprint_columns_to_local_salt_top": column_envelope(state),
        shell_name: shell,
    }
    statistics = {
        name: _envelope_statistics(
            name=name,
            envelope=envelope,
            salt=salt,
            external_air=external_air,
            boundaries=boundaries,
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
    fixed_window_count = int(np.count_nonzero(fixed_window))
    fixed_window_solid_fraction = (
        int(np.count_nonzero(fixed_window & salt)) / fixed_window_count
        if fixed_window_count
        else None
    )
    if intrinsic_material_mask is None:
        identifiable = False
        intrinsic_porosity = None
        flag = "undefined material boundary; intrinsic porosity not reported"
        ambiguous_exterior = statistics[
            "axis_aligned_salt_bounding_box"
        ].external_air_connected_void_voxels
    else:
        if (
            intrinsic_material_mask.shape != state.shape
            or intrinsic_material_mask.dtype != np.bool_
        ):
            raise ValueError(
                "intrinsic_material_mask must be boolean and match state shape"
            )
        if np.any(
            intrinsic_material_mask & (state.material == Material.SUBSTRATE)
        ):
            raise ValueError("intrinsic material boundary cannot include substrate")
        if not intrinsic_boundary_source:
            raise ValueError(
                "intrinsic_boundary_source is required with intrinsic_material_mask"
            )
        material_count = int(np.count_nonzero(intrinsic_material_mask))
        material_void = int(
            np.count_nonzero(
                intrinsic_material_mask & (state.material == Material.EMPTY)
            )
        )
        identifiable = True
        intrinsic_porosity = material_void / material_count if material_count else None
        flag = f"reported within explicit boundary: {intrinsic_boundary_source}"
        ambiguous_exterior = 0

    if grain_ids is None:
        accessible_intergranular = None
    else:
        gap_mask = face_adjacent_intergranular_gap_mask(
            state, grain_ids, analysis_mask=fixed_window
        )
        accessible_intergranular = int(np.count_nonzero(gap_mask & external_air))

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
        face_shell_envelope=statistics[shell_name],
        salt_void_surface_faces=surface_faces,
        salt_surface_area_m2=surface_area,
        salt_volume_m3=volume,
        surface_to_volume_ratio_m_inv=ratio,
        characteristic_solid_size_m=characteristic_size,
        tested_envelope_void_fraction_range=fraction_range,
        tested_envelope_pathway_agreement=pathway_agreement,
        intrinsic_porosity_identifiable=identifiable,
        intrinsic_porosity_flag=flag,
        fixed_window_solid_volume_fraction=fixed_window_solid_fraction,
        intrinsic_porosity=intrinsic_porosity,
        accessible_intergranular_void_voxels=accessible_intergranular,
        accessible_intergranular_void_volume_m3=(
            None
            if accessible_intergranular is None or voxel_size_m is None
            else accessible_intergranular * voxel_size_m**3
        ),
        ambiguous_exterior_void_voxels=ambiguous_exterior,
        ambiguous_exterior_void_volume_m3=(
            None
            if voxel_size_m is None
            else ambiguous_exterior * voxel_size_m**3
        ),
        physical_shell_thickness_m=recorded_shell_thickness,
    )
