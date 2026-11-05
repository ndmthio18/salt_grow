# Numerical provenance

## Retained development events

- The first environment-check command failed before execution because of an
  empty PowerShell pipe element. No files or model state were changed.
- The first pytest invocation used an incorrect nested `PYTHONPATH`; collection
  failed before tests ran. The corrected invocation passed all tests.
- A 240-cell, 6 h B2 refinement was started after the 80-, 120-, and 160-cell
  runs completed. It was manually stopped before producing a registered result
  because SciPy's finite-difference BDF Jacobian became a costly dense solve and
  triggered high multithreaded BLAS use. No partial `.npz` or summary row was
  accepted. The replacement spatial check keeps 160 cells but increases the ALE
  top-clustering power from 2.0 to 2.5; this directly refines the only sharp
  gradient while retaining the full-domain conservation equations.
- A subsequent clustering-power 3.5 run was manually stopped after SciPy's
  finite-difference Jacobian step-size estimator emitted overflow warnings.
  The completed power-2.5 and power-3.0 runs already changed terminal velocity
  by less than 1%, so the more extreme clustering was neither registered nor
  needed for the stated acceptance gate.
- The first `k_r = 2.3e-4 m/s` sensitivity attempt used the accepted baseline's
  cubic clustering and was stopped after repeated SciPy finite-difference
  Jacobian step-size overflow warnings. Because lowering `k_r` thickens the
  analytical reaction layer, the S3 low-`k_r` cases were rerun on the retained
  160-cell quadratic grid and were accepted only if their own resolution and
  conservation diagnostics passed.
- Initial-height sensitivities on `grid_power = 3.0` were stopped after the
  same Jacobian step-size warnings. The smallest cells were then far below both
  the modeled reaction scale and the estimated crystal scale. The accepted
  baseline was set to the coarsest mesh passing the predefined observable and
  resolution gates (`power = 2.5`), and S6 was run with quadratic clustering,
  subject to each run's own `cells_across_delta` check.

The stopped 240-cell attempt is a performance limitation, not evidence of
numerical failure or physical behavior. It must not be presented as a completed
convergence point.
