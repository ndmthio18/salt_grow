# Progress

## 2026-11-06 — Phase 1 complete

- Audited the pinned external CA reference without copying its unlicensed
  source into this repository.
- Recorded the unavailable NVIDIA/CUDA/PyCUDA baseline honestly; no original
  GPU simulation or morphology result was claimed.
- Implemented the independent SI-unit evaporation, NaCl precipitation, and
  mass-to-voxel accounting modules.
- Verified equilibrium/high-RH clipping, analytical balances, timestep
  partition invariance, input validation, sub-voxel carry-forward, multiple
  voxel sizes, initial-seed exclusion, and cumulative mass closure.
- Generated and visually checked the deterministic numerical tables and
  RH/flux/rate/cumulative-mass figure.
- Passed 34 Phase 1 tests and 50 combined repository tests.

Status: Phase 1 accounting module GO. Phase 2 remains on hold for review and
the explicit licensing, CUDA-baseline, calibration, and spatial-allocation
gates in `docs/phase1_feasibility_and_validation.md`.
