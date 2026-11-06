"""Controlled Phase 2A.1 simulation and comparison helpers."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from time import perf_counter
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from salt_phase1 import EvaporationParameters, compute_mass_budget

from .ca import AttachmentWeights, SaltGrowthCA
from .morphology import MorphologyAnalysis, analyze_morphology
from .state import Coordinate, Material, initialize_physical_state


@dataclass(frozen=True)
class CaseConfiguration:
    case_id: str
    rule_name: str
    shape: tuple[int, int, int]
    substrate_layers: int
    seed_coordinates: tuple[Coordinate, ...]
    wetted_seed_indices: tuple[int, ...] | None
    enforce_wetting: bool
    growth_mode: str
    update_mode: str
    attachment_weights: AttachmentWeights
    rng_seed: int
    timestep_s: float
    total_duration_s: float
    evaporation_area_scale: float = 1.0
    voxel_size_m: float = 10.0e-6
    solid_density_kg_m3: float = 2160.0


@dataclass
class CaseResult:
    configuration: CaseConfiguration
    summary: dict[str, object]
    history: list[dict[str, object]]
    final_salt_mask: NDArray[np.bool_]
    seed_mask: NDArray[np.bool_]
    snapshots: dict[int, NDArray[np.bool_]]
    final_analysis: MorphologyAnalysis


def _mass_parameters(config: CaseConfiguration) -> EvaporationParameters:
    area = (
        config.shape[0]
        * config.voxel_size_m
        * config.shape[1]
        * config.voxel_size_m
        * config.evaporation_area_scale
    )
    return EvaporationParameters(
        temperature_k=298.15,
        water_activity=0.753,
        relative_humidity=0.30,
        mass_transfer_coefficient_m_s=0.005,
        evaporation_area_m2=area,
        timestep_s=config.timestep_s,
        salt_mass_fraction=0.264,
        solid_density_kg_m3=config.solid_density_kg_m3,
    )


def run_case(config: CaseConfiguration) -> CaseResult:
    """Run one auditable case and retain conservation and morphology histories."""

    step_ratio = config.total_duration_s / config.timestep_s
    if not math.isclose(step_ratio, round(step_ratio), rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("total_duration_s must be divisible by timestep_s")
    steps = int(round(step_ratio))
    state = initialize_physical_state(
        shape=config.shape,
        substrate_layers=config.substrate_layers,
        seed_coordinates=config.seed_coordinates,
        wetted_seed_indices=config.wetted_seed_indices,
    )
    ca = SaltGrowthCA(
        state,
        solid_density_kg_m3=config.solid_density_kg_m3,
        voxel_size_m=config.voxel_size_m,
        enforce_wetting=config.enforce_wetting,
        growth_mode=config.growth_mode,
        update_mode=config.update_mode,
        attachment_weights=config.attachment_weights,
        rng_seed=config.rng_seed,
        tie_break_seed=config.rng_seed,
    )
    mass_budget = compute_mass_budget(_mass_parameters(config))
    snapshot_steps = {0, max(1, steps // 2), steps}
    snapshots: dict[int, NDArray[np.bool_]] = {
        0: (state.material == Material.SALT).copy()
    }
    history: list[dict[str, object]] = []
    initial_analysis = analyze_morphology(
        state, supply_connected_salt=ca.supply_connected_salt_mask()
    )
    previous_cluster_count = initial_analysis.salt_cluster_count
    coalescence_events = 0
    coalescence_steps: list[int] = []
    started = perf_counter()
    for step in range(1, steps + 1):
        result = ca.grow(mass_budget.precipitated_salt_mass_kg)
        analysis = analyze_morphology(
            state, supply_connected_salt=ca.supply_connected_salt_mask()
        )
        decrease = max(previous_cluster_count - analysis.salt_cluster_count, 0)
        if decrease:
            coalescence_events += decrease
            coalescence_steps.append(step)
        previous_cluster_count = analysis.salt_cluster_count
        direct_new_voxels = int(np.count_nonzero(state.material == Material.SALT)) - len(
            config.seed_coordinates
        )
        if direct_new_voxels != ca.allocator.cumulative_realized_solid_voxels:
            raise RuntimeError("direct geometric inventory disagrees with mass ledger")
        history.append(
            {
                "case_id": config.case_id,
                "step": step,
                "time_s": step * config.timestep_s,
                "incoming_mass_kg": mass_budget.precipitated_salt_mass_kg,
                "cumulative_supplied_mass_kg": ca.allocator.cumulative_input_mass_kg,
                "committed_new_salt_mass_kg": (
                    ca.allocator.cumulative_committed_solid_mass_kg
                ),
                "retained_mass_kg": (
                    ca.allocator.deferred_precipitation_mass_kg
                    + ca.allocator.mass_remaining_in_solution_kg
                ),
                "mass_closure_error_kg": ca.allocator.closure_error_kg(),
                "requested_voxels": result.transaction.requested_solid_voxels,
                "realized_voxels": result.transaction.realized_solid_voxels,
                "cumulative_new_voxels": direct_new_voxels,
                "candidate_voxels_at_step_start": result.candidate_voxel_count,
                "maximum_height_voxels": (
                    result.morphology.maximum_height_voxels_above_substrate
                ),
                "footprint_voxels": result.morphology.footprint_voxels,
                "salt_cluster_count": analysis.salt_cluster_count,
                "supply_connected_salt_voxels": (
                    analysis.supply_connected_salt_voxels
                ),
                "transaction_status": result.transaction.status,
            }
        )
        if step in snapshot_steps:
            snapshots[step] = (state.material == Material.SALT).copy()
    runtime_s = perf_counter() - started
    ca.assert_invariants()
    final_connected = ca.supply_connected_salt_mask()
    final_analysis = analyze_morphology(
        state, supply_connected_salt=final_connected
    )
    morphology = ca.morphology_diagnostics()
    retained_mass = (
        ca.allocator.deferred_precipitation_mass_kg
        + ca.allocator.mass_remaining_in_solution_kg
    )
    unrealized_whole_voxels = int(retained_mass // ca.voxel_mass_kg)
    total_salt = morphology.total_salt_voxels
    summary: dict[str, object] = {
        "case_id": config.case_id,
        "rule_name": config.rule_name,
        "growth_mode": config.growth_mode,
        "update_mode": config.update_mode,
        "rng_seed": config.rng_seed,
        "w1": config.attachment_weights.w1,
        "w2": config.attachment_weights.w2,
        "w3": config.attachment_weights.w3,
        "w4": config.attachment_weights.w4,
        "timestep_s": config.timestep_s,
        "total_duration_s": config.total_duration_s,
        "evaporation_area_scale": config.evaporation_area_scale,
        "initial_seed_count": len(config.seed_coordinates),
        "enforce_wetting": config.enforce_wetting,
        "runtime_s": runtime_s,
        "total_supplied_mass_kg": ca.allocator.cumulative_input_mass_kg,
        "committed_new_salt_mass_kg": (
            ca.allocator.cumulative_committed_solid_mass_kg
        ),
        "retained_mass_kg": retained_mass,
        "mass_closure_error_kg": ca.allocator.closure_error_kg(),
        "new_salt_voxels": morphology.new_salt_voxels,
        "total_salt_voxels": total_salt,
        "new_salt_volume_m3": morphology.new_salt_voxels * config.voxel_size_m**3,
        "total_salt_volume_m3": total_salt * config.voxel_size_m**3,
        "maximum_height_voxels": morphology.maximum_height_voxels_above_substrate,
        "footprint_voxels": morphology.footprint_voxels,
        "salt_cluster_count": final_analysis.salt_cluster_count,
        "salt_cluster_sizes": json.dumps(final_analysis.salt_cluster_sizes),
        "seed_occupied_cluster_count": final_analysis.seed_occupied_cluster_count,
        "coalescence_events": coalescence_events,
        "coalescence_steps": json.dumps(coalescence_steps),
        "supply_connected_salt_voxels": (
            final_analysis.supply_connected_salt_voxels
        ),
        "supply_connected_fraction": (
            final_analysis.supply_connected_salt_voxels / total_salt
            if total_salt
            else 0.0
        ),
        "bbox_salt_fraction": final_analysis.bounding_box.salt_volume_fraction,
        "bbox_void_components": final_analysis.bounding_box.void_component_count,
        "bbox_enclosed_void_voxels": (
            final_analysis.bounding_box.enclosed_void_voxels
        ),
        "bbox_external_void_voxels": (
            final_analysis.bounding_box.external_air_connected_void_voxels
        ),
        "column_salt_fraction": final_analysis.column_envelope.salt_volume_fraction,
        "column_void_components": final_analysis.column_envelope.void_component_count,
        "column_enclosed_void_voxels": (
            final_analysis.column_envelope.enclosed_void_voxels
        ),
        "column_external_void_voxels": (
            final_analysis.column_envelope.external_air_connected_void_voxels
        ),
        "unrealized_whole_voxels": unrealized_whole_voxels,
        "full_whole_voxel_realization": unrealized_whole_voxels == 0,
    }
    return CaseResult(
        configuration=config,
        summary=summary,
        history=history,
        final_salt_mask=(state.material == Material.SALT).copy(),
        seed_mask=state.seed_mask.copy(),
        snapshots=snapshots,
        final_analysis=final_analysis,
    )


def mask_jaccard(first: NDArray[np.bool_], second: NDArray[np.bool_]) -> float:
    if first.shape != second.shape:
        raise ValueError("masks must have matching shapes")
    union = int(np.count_nonzero(first | second))
    if union == 0:
        return 1.0
    return int(np.count_nonzero(first & second)) / union


def mask_difference_fraction(
    first: NDArray[np.bool_], second: NDArray[np.bool_]
) -> float:
    if first.shape != second.shape:
        raise ValueError("masks must have matching shapes")
    return float(np.mean(first != second))
