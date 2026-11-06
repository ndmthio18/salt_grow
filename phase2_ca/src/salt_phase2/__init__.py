"""Independent CPU Phase 2A salt-growth cellular automaton."""

from .ca import AttachmentWeights, GrowthStepResult, MorphologyDiagnostics, SaltGrowthCA
from .morphology import (
    EnvelopeVoidStatistics,
    MorphologyAnalysis,
    analyze_morphology,
    bounding_box_envelope,
    column_envelope,
    computational_box_envelope,
    connected_components,
    external_air_connected_mask,
    face_shell_envelope,
)
from .state import Coordinate, Material, PhysicalState, initialize_physical_state
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
    "EnvelopeVoidStatistics",
    "Material",
    "MorphologyDiagnostics",
    "MorphologyAnalysis",
    "PhysicalState",
    "SaltGrowthCA",
    "analyze_morphology",
    "bounding_box_envelope",
    "column_envelope",
    "computational_box_envelope",
    "connected_components",
    "external_air_connected_mask",
    "face_shell_envelope",
    "initialize_physical_state",
    "mask_difference_fraction",
    "mask_jaccard",
    "run_case",
]
