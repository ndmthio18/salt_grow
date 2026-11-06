# Phase 1 evaporation-driven salt mass budget

This directory is an independent CPU implementation of the Phase 1 physical
accounting model. It does not import, copy, or modify the reference CUDA CA.

Run the deterministic tests from the repository root:

```powershell
python -m pytest phase1_mass_budget/tests -q
```

Regenerate the reference tables and plots:

```powershell
python phase1_mass_budget/run_reference_case.py
```

The model uses SI units throughout. Existing seed voxels are morphology
metadata and are never counted as newly precipitated salt. Spatial placement of
new voxels is intentionally deferred to Phase 2.

Phase 1.1 adds a scalar domain-capacity diagnostic and a transactional
mass-allocation interface for a future geometry layer. It still contains no CA
growth implementation. See `docs/phase1_1_spatial_growth_readiness.md` for the
reference capacity audit, proposed Phase 2 state representation, and bounded
GO/NO-GO decision.
