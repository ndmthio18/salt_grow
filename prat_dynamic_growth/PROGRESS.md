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

## 2026-11-05 — M3 controlled numerical attribution

- Implemented Q0--Q6 at 5 min, 1 h, 3 h and 6 h, plus sink-only and
  boundary-only nonlinear controls; retained precipitation-induced pore-volume
  rates in the frozen concentration equation.
- Corrected the inconsistent ALE face porosity in salt versus liquid storage;
  replaced the unphysical-root fallback with failure and a stable quadratic.
- Implemented T0--T3, resolved startup through 1 s, independent local/interface
  residuals, separate trajectory-derivative audits, and configuration/source/data
  hash validation before cache reuse.
- Completed 160/240/320/480/640-cell and 120/60/30-s convergence, eleven distinct
  6 h runs, all finishing without stop events. The final T0 BDF subtraction
  warning is retained in the manifest; the other warning records are empty.
- Revised the baseline difference from +3.325% to +2.226%. At 6 h, frozen
  porosity contributes +2.612 points, diffusivity +1.160, area -1.588, and
  their interaction +0.031. Initialization does not explain the persistent gap.
- Passed 13 automated tests and the successive-grid, algebraic local/interface,
  and global balance gates. The additional trajectory-differentiation audit
  does not uniformly meet the stricter 1e-8 local / 1e-10 interface criteria;
  conditioning and unresolved tiny effects are explicitly recorded in RESULTS.
- Archived data and six PNG/SVG figure pairs under `outputs/m3/`; updated
  RESULTS, M3_METHODS, README and GO_NO_GO. No morphology/hydraulic extension
  or physical-validation claim was made.

Current status: requested numerical attribution complete; percent-level
closure effects resolved, strict trajectory-level residual qualification
remains limited as documented.
