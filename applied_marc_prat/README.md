# Applied Marc–Prat snapshot from liqin

This directory is a curated, provenance-preserving copy of the useful Marc
Prat salt-growth work developed on `liqin`. It is retained as a historical
reference implementation and validation chain. It has **not** been merged into
the repository's current `prat_dynamic_growth` model.

## What is here

- `prat2026/`: analytical equations, parameters, digitised data, figure
  reproductions, and formula-level tests for Prat (2026).
- `saltv2/`: the applied package containing material properties, pore-network
  transport, vapour/growth kernels, and 2-D/3-D coupling drivers.
- `tests/`: the corresponding unit, conservation, interface, regression, and
  Step 0/1/2A verification tests.
- `runs/S_surface/`: the surface-first verification chain, including the
  Step 2A independent one-dimensional reaction–transport reference.
- `knowledge/` and `ideas/`: assumptions, model design, and earlier critical
  analyses needed to interpret the code.
- `reviews/Prat生长律适用性审查_2026-10-25.md`: a later applicability
  review copied from the main liqin development branch because it materially
  limits the scientific claims that can be made from `saltv2`.
- `branch_archive/`: complete 14-branch Git bundle, browsable deltas for the
  seven independent branch tips, and the uncommitted hierarchical overlay.
- `MIGRATION_ANALYSIS.md`: current assessment of what is reusable and what is
  not validated.
- `PROVENANCE.md`: exact source commits, integrity checks, exclusions, and
  local verification results.

## Recommended entry points

1. Read `MIGRATION_ANALYSIS.md` before reusing a closure or result.
2. Read `runs/S_surface/step2a_prat/README.md` for the verified one-dimensional
   Prat problem and its assumptions.
3. Read `reviews/Prat生长律适用性审查_2026-10-25.md` before interpreting the
   3-D morphology or patchy–crusty behaviour.
4. Use `prat2026/core.py` and `saltv2/growth/prat.py` as the analytical-law
   entry points; use `saltv2/coupling/driver2d.py` and `driver3d.py` only after
   reviewing their additional closures.
5. Read `branch_archive/BRANCH_INDEX.md` for the competing structure,
   audit, layered-reference, beta/hierarchy, and Stage-T routes.

## Minimal local verification

On Windows, force UTF-8 because three architecture tests call
`Path.read_text()` without specifying the source encoding:

```powershell
cd "D:\3OutputFile\1Research\Article\23salt grow\applied_marc_prat"
$env:PYTHONUTF8 = '1'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python -m pytest -q tests/test_layering.py tests/test_step2a_prat.py prat2026/tests
```

On 2026-11-12 this gave `48 passed in 2.58 s` with Python 3.12.4 on
Windows. The complete suite gave `269 passed, 1 skipped, 6 failed`; the six
failures and their interpretation are recorded in `PROVENANCE.md`.

## Scientific scope

The strongest result in this snapshot is numerical verification of Prat's
analytical growth law against an independent one-dimensional reference under
Prat's own linearised, saturated-supply assumptions. It is not experimental
validation of the three-dimensional local law, nucleation, capillary supply,
crust formation, pore-scale morphology, or a patchy–crusty transition.

