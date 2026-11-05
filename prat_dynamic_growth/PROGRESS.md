# Progress

## 2026-11-05

- Located and verified the Prat (2026) source PDF; recorded its SHA-256.
- Implemented Stage A fixed-property finite-volume verification and passed the
  `1e-5` velocity-error gate.
- Implemented conservative ALE Model B with independent kinetic interface
  velocity, birth porosity, B1/B2 closures, signed precipitation/dissolution,
  and separate mass inventories.
- Passed automated tests (`5 passed`) and global conservation gates.
- Completed baseline, grid/time refinement, and the smallest sensitivity matrix
  covering `d_sub`, evaporation, `k_r`, closure choice, and initial height.
- Accepted the 160-cell `grid_power=2.5` baseline after comparison with the
  more strongly clustered grid; retained stopped/retried runs in
  `NUMERICAL_PROVENANCE.md`.
- Generated all six required figures in PNG and PDF and visually reviewed the
  final versions.
- Completed `README.md`, `MODEL_ASSUMPTIONS.md`, `RESULTS.md`, and
  `GO_NO_GO.md`.

Current status: prototype complete. Physical validation is not claimed.

## 2026-11-05 — M2 critical review

- Audited the existing ALE equations, top boundary, inventories, campaign
  controls, and verification tests without changing solver code.
- Identified constant birth porosity as the self-replication closure and showed
  quantitatively why wake densification has weak direct leverage on velocity.
- Defined the Q0--Q6 frozen-snapshot attribution matrix and T0--T3 initial-state
  controls required to isolate the approximately 3% deviation.
- Proposed a conservative morphology-state front closure and a minimal
  permeability/capillary supply constraint.
- Specified sparse-solver and independent local/interface residual acceptance
  gates in `M2_CRITICAL_REVIEW.md`.

Current M2 status: review complete; additional models are intentionally not yet
implemented pending the governing-equation and acceptance-test gate.
