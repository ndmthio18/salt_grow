"""Synthetic voxel geometries with analytical void-connectivity answers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .morphology import (
    PhysicalBoundarySurfaces,
    analyze_void_space,
    external_air_connected_mask,
    face_adjacent_intergranular_gap_mask,
    rectangular_boundary_surfaces,
    rectangular_window_mask,
)
from .state import Material, PhysicalState


SHAPE = (9, 9, 9)
LOWER = (1, 1, 1)
UPPER_EXCLUSIVE = (8, 8, 8)


@dataclass(frozen=True)
class GeometryBenchmark:
    name: str
    state: PhysicalState
    fixed_window: NDArray[np.bool_]
    intended_material_mask: NDArray[np.bool_]
    boundaries: PhysicalBoundarySurfaces
    grain_ids: NDArray[np.int32] | None
    expected_window_bottom_to_top: bool
    expected_material_bottom_to_top: bool
    expected_material_enclosed_void_voxels: int
    expected_material_external_void_voxels: int
    expected_accessible_intergranular_void_voxels: int | None


def _state_from_salt(
    salt_mask: NDArray[np.bool_],
    *,
    seed_coordinate: tuple[int, int, int],
) -> PhysicalState:
    material = np.full(SHAPE, Material.EMPTY, dtype=np.uint8)
    material[:, :, 0] = Material.SUBSTRATE
    material[salt_mask] = Material.SALT
    seed_mask = np.zeros(SHAPE, dtype=bool)
    seed_mask[seed_coordinate] = True
    liquid_accessible = np.zeros(SHAPE, dtype=bool)
    liquid_accessible[seed_coordinate] = True
    state = PhysicalState(
        material=material,
        liquid_accessible=liquid_accessible,
        seed_mask=seed_mask,
        substrate_layers=1,
    )
    state.validate()
    return state


def _common_masks() -> tuple[NDArray[np.bool_], PhysicalBoundarySurfaces]:
    window = rectangular_window_mask(
        SHAPE, lower=LOWER, upper_exclusive=UPPER_EXCLUSIVE
    )
    boundaries = rectangular_boundary_surfaces(
        SHAPE,
        lower=LOWER,
        upper_exclusive=UPPER_EXCLUSIVE,
        name="benchmark_fixed_window_surfaces",
    )
    return window, boundaries


def build_geometry_benchmarks() -> tuple[GeometryBenchmark, ...]:
    """Return five independent geometries with prespecified reference answers."""

    window, boundaries = _common_masks()
    cases: list[GeometryBenchmark] = []

    pillar = np.zeros(SHAPE, dtype=bool)
    pillar[4, 4, 1:8] = True
    cases.append(
        GeometryBenchmark(
            name="solid_pillar_exterior_bypass",
            state=_state_from_salt(pillar, seed_coordinate=(4, 4, 1)),
            fixed_window=window,
            intended_material_mask=pillar.copy(),
            boundaries=boundaries,
            grain_ids=None,
            expected_window_bottom_to_top=True,
            expected_material_bottom_to_top=False,
            expected_material_enclosed_void_voxels=0,
            expected_material_external_void_voxels=0,
            expected_accessible_intergranular_void_voxels=None,
        )
    )

    slab_mask = rectangular_window_mask(
        SHAPE, lower=(2, 2, 1), upper_exclusive=(7, 7, 8)
    )
    porous_slab = slab_mask.copy()
    porous_slab[4, 4, 1:8] = False
    cases.append(
        GeometryBenchmark(
            name="porous_slab_through_channel",
            state=_state_from_salt(porous_slab, seed_coordinate=(2, 2, 1)),
            fixed_window=window,
            intended_material_mask=slab_mask,
            boundaries=boundaries,
            grain_ids=None,
            expected_window_bottom_to_top=True,
            expected_material_bottom_to_top=True,
            expected_material_enclosed_void_voxels=0,
            expected_material_external_void_voxels=7,
            expected_accessible_intergranular_void_voxels=None,
        )
    )

    two_grain_material = rectangular_window_mask(
        SHAPE, lower=(2, 2, 1), upper_exclusive=(7, 7, 8)
    )
    grain_one = np.zeros(SHAPE, dtype=bool)
    grain_two = np.zeros(SHAPE, dtype=bool)
    grain_one[2:4, 2:7, 1:8] = True
    grain_two[5:7, 2:7, 1:8] = True
    two_grain_salt = grain_one | grain_two
    grain_ids = np.zeros(SHAPE, dtype=np.int32)
    grain_ids[grain_one] = 1
    grain_ids[grain_two] = 2
    cases.append(
        GeometryBenchmark(
            name="two_grains_connected_intergranular_gap",
            state=_state_from_salt(two_grain_salt, seed_coordinate=(2, 2, 1)),
            fixed_window=window,
            intended_material_mask=two_grain_material,
            boundaries=boundaries,
            grain_ids=grain_ids,
            expected_window_bottom_to_top=True,
            expected_material_bottom_to_top=True,
            expected_material_enclosed_void_voxels=0,
            expected_material_external_void_voxels=35,
            expected_accessible_intergranular_void_voxels=35,
        )
    )

    closed_material = rectangular_window_mask(
        SHAPE, lower=(2, 2, 1), upper_exclusive=(7, 7, 7)
    )
    closed_salt = closed_material.copy()
    closed_salt[4, 4, 4] = False
    cases.append(
        GeometryBenchmark(
            name="closed_internal_cavity",
            state=_state_from_salt(closed_salt, seed_coordinate=(2, 2, 1)),
            fixed_window=window,
            intended_material_mask=closed_material,
            boundaries=boundaries,
            grain_ids=None,
            expected_window_bottom_to_top=True,
            expected_material_bottom_to_top=False,
            expected_material_enclosed_void_voxels=1,
            expected_material_external_void_voxels=0,
            expected_accessible_intergranular_void_voxels=None,
        )
    )

    blind_material = closed_material.copy()
    blind_salt = blind_material.copy()
    blind_salt[4, 4, 4:7] = False
    cases.append(
        GeometryBenchmark(
            name="externally_connected_blind_pore",
            state=_state_from_salt(blind_salt, seed_coordinate=(2, 2, 1)),
            fixed_window=window,
            intended_material_mask=blind_material,
            boundaries=boundaries,
            grain_ids=None,
            expected_window_bottom_to_top=True,
            expected_material_bottom_to_top=False,
            expected_material_enclosed_void_voxels=0,
            expected_material_external_void_voxels=3,
            expected_accessible_intergranular_void_voxels=None,
        )
    )
    return tuple(cases)


def evaluate_geometry_benchmark(case: GeometryBenchmark) -> dict[str, object]:
    """Evaluate one benchmark without substituting exterior access for a path."""

    salt = case.state.material == Material.SALT
    external_air = external_air_connected_mask(case.state)
    window = analyze_void_space(
        name="fixed_window",
        analysis_mask=case.fixed_window,
        salt_mask=salt,
        external_air_mask=external_air,
        boundaries=case.boundaries,
    )
    material = analyze_void_space(
        name="intended_material",
        analysis_mask=case.intended_material_mask,
        salt_mask=salt,
        external_air_mask=external_air,
        boundaries=case.boundaries,
    )
    if case.grain_ids is None:
        intergranular_count = None
    else:
        intergranular = face_adjacent_intergranular_gap_mask(
            case.state,
            case.grain_ids,
            analysis_mask=case.intended_material_mask,
        )
        intergranular_count = int(np.count_nonzero(intergranular & external_air))
    return {
        "benchmark": case.name,
        "window_bottom_to_top": window.bottom_to_top_void_percolates,
        "material_bottom_to_top": material.bottom_to_top_void_percolates,
        "exterior_bypass_detected": (
            window.bottom_to_top_void_percolates
            and not material.bottom_to_top_void_percolates
        ),
        "material_external_void_voxels": (
            material.external_air_connected_void_voxels
        ),
        "material_enclosed_void_voxels": material.enclosed_void_voxels,
        "accessible_intergranular_void_voxels": intergranular_count,
        "expected_window_bottom_to_top": case.expected_window_bottom_to_top,
        "expected_material_bottom_to_top": case.expected_material_bottom_to_top,
        "expected_material_external_void_voxels": (
            case.expected_material_external_void_voxels
        ),
        "expected_material_enclosed_void_voxels": (
            case.expected_material_enclosed_void_voxels
        ),
        "expected_accessible_intergranular_void_voxels": (
            case.expected_accessible_intergranular_void_voxels
        ),
    }
