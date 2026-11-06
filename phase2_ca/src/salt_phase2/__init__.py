"""Independent CPU Phase 2A salt-growth cellular automaton."""

from .ca import AttachmentWeights, GrowthStepResult, MorphologyDiagnostics, SaltGrowthCA
from .morphology import (
    EnvelopeVoidStatistics,
    MorphologyAnalysis,
    PhysicalBoundarySurfaces,
    analyze_morphology,
    analyze_void_space,
    bounding_box_envelope,
    column_envelope,
    computational_box_envelope,
    connected_components,
    domain_boundary_surfaces,
    external_air_connected_mask,
    face_adjacent_intergranular_gap_mask,
    face_shell_envelope,
    rectangular_boundary_surfaces,
    rectangular_window_mask,
)
from .state import Coordinate, Material, PhysicalState, initialize_physical_state
from .geometry_benchmarks import (
    GeometryBenchmark,
    build_geometry_benchmarks,
    evaluate_geometry_benchmark,
)
from .validation import (
    CaseConfiguration,
    CaseResult,
    mask_difference_fraction,
    mask_jaccard,
    run_case,
)

__all__ = [
    "AttachmentWeights",
    "Coordinate",
    "CaseConfiguration",
    "CaseResult",
    "GrowthStepResult",
    "GeometryBenchmark",
    "EnvelopeVoidStatistics",
    "Material",
    "MorphologyDiagnostics",
    "MorphologyAnalysis",
    "PhysicalBoundarySurfaces",
    "PhysicalState",
    "SaltGrowthCA",
    "analyze_morphology",
    "analyze_void_space",
    "build_geometry_benchmarks",
    "bounding_box_envelope",
    "column_envelope",
    "computational_box_envelope",
    "connected_components",
    "domain_boundary_surfaces",
    "external_air_connected_mask",
    "evaluate_geometry_benchmark",
    "face_adjacent_intergranular_gap_mask",
    "face_shell_envelope",
    "rectangular_boundary_surfaces",
    "rectangular_window_mask",
    "initialize_physical_state",
    "mask_difference_fraction",
    "mask_jaccard",
    "run_case",
]
