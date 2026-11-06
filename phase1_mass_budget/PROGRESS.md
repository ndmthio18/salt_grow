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

## 2026-11-06 — Phase 1.1 complete

- Added scalar domain-capacity diagnostics with configurable porosity reserve,
  seed accounting, overflow status, and constant-rate capacity time.
- Verified that the original 64-cubed, 10-um, RH=0.30, 24 h case requires
  `116.924%` of the domain and reaches full-solid capacity at `20.526 h`.
- Added a two-stage mass allocator: geometry requests do not become committed
  solid until the realized voxel count is confirmed; all other mass remains
  explicitly deferred or in solution.
- Removed the hard-coded 1441-point time grid and made SVG identifiers and
  metadata deterministic within the documented environment.
- Added the Phase 2 material/phase/eligibility design without implementing the
  CA, transport, permeability, or vapor solver.
- Passed 53 Phase 1/1.1 tests and all 69 combined repository tests.
- Confirmed remote NVIDIA hardware on frodo, but the exact reference-CA
  environment remains unavailable because CUDA Toolkit and PyCUDA are absent.

Status: conditional GO for a bounded independent CPU Phase 2 prototype after
this task; NO-GO for the original 64-cubed 24 h geometry and for physical
validation claims.
