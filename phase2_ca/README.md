# Phase 2A minimal 3D salt-growth CA

This directory contains an independent CPU prototype that couples the verified
Phase 1 evaporation mass supply and transactional voxel ledger to a minimal 3D
geometry update.

Run its tests from the repository root:

```powershell
python -m pytest phase2_ca/tests -q
```

Run the deterministic demonstration:

```powershell
python phase2_ca/run_demo.py
```

Run the Phase 2A.1 stochastic rule, timestep, nucleation, wetting, and void
morphology validation matrix:

```powershell
python phase2_ca/run_phase2a1_validation.py
```

The campaign outputs raw conservation histories, morphology tables, 3D Boolean
masks, figures, provenance, and hashes under `outputs/phase2a1/`. See
`docs/phase2a_1_growth_rule_validation.md` for the evidence-bounded assessment.

The material and liquid-accessibility arrays are separate. Salt-film supply is
represented only by six-face connectivity to a wetted seed. An empty air voxel
can grow when it is adjacent to supply-connected salt; the target itself need
not contain liquid.

Candidates are ranked by local six-face salt coordination. Ties use a pinned,
mixed-coordinate priority to avoid a systematic low-height axis bias while
remaining reproducible. This priority is an algorithmic prototype rule, not a
validated NaCl crystal-growth law.

Phase 2A.1 also provides weighted stochastic selection and sequential
event-based updates. Relative weights change where fixed-budget salt is placed;
they are not calibrated attachment probabilities and are not equivalent to the
reference Redkov Monte Carlo process.

The wetting filter can be disabled as an explicitly nonphysical control. This
prototype contains no pore-network flow, film hydraulics, vapor diffusion,
permeability, morphology-dependent evaporation, dissolution, or experimental
calibration. Numerical mass closure does not validate the growth physics.

Run the Phase 2A.3a surface-normal and crystallographic-orientation verification
in its pre-registered order with:

```powershell
python phase2_ca/run_phase2a3a_verification.py --stage development --config phase2_ca/configs/phase2a3a_method_development.json --output phase2_ca/outputs/phase2a3a_development_reproduced
python phase2_ca/run_phase2a3a_verification.py --stage locked --config phase2_ca/configs/phase2a3a_locked_protocol.json --output phase2_ca/outputs/phase2a3a_locked_reproduced
```

The locked parameters and thresholds are immutable. The evidence-bounded
assessment is `docs/phase2a_3a_normals_orientation_verification.md`. The result
is conditional GO for a resolved single-grain orientation pilot only; one- and
two-voxel seeds, multi-grain growth, and transport coupling remain NO-GO.

Run the Phase 2A.3b-0 candidate-normal and anisotropic-growth interface
validation with:

```powershell
python phase2_ca/run_phase2a3b0_validation.py --config phase2_ca/configs/phase2a3b0_protocol.json --output phase2_ca/outputs/phase2a3b0_reproduced
```

The campaign compares two predeclared edge/corner policies, physical-scale
transfer, full/local/cached normal updates, alpha-zero equivalence, pinned RNG
reproducibility, and joint-rotation covariance. See
`docs/phase2a_3b0_growth_interface_validation.md`. The GO is limited to a
resolved single-grain numerical pilot with both edge policies retained.

Run the Phase 2A.3b-1 method-development gate and locked anisotropic morphology
verification with:

```powershell
python phase2_ca/run_phase2a3b1_validation.py --stage development --config phase2_ca/configs/phase2a3b1_method_development.json --output phase2_ca/outputs/phase2a3b1_development_reproduced
python phase2_ca/run_phase2a3b1_validation.py --stage locked --config phase2_ca/configs/phase2a3b1_locked_verification.json --output phase2_ca/outputs/phase2a3b1_locked_reproduced
```

See `docs/phase2a_3b1_anisotropic_growth_mechanism.md`. Fixed-front probability,
cache correctness, mass closure, rotation covariance, and a separate anchored
supply control pass. The registered 10% growth-only morphology-signal gate
fails, so the complete single-grain morphology model and all multigrain or
transport extensions remain NO-GO.

Run the Phase 2A.3b-2 weak-anisotropy diagnosis in its gated order with:

```powershell
python phase2_ca/run_phase2a3b2_diagnosis.py --stage exploratory --config phase2_ca/configs/phase2a3b2_exploratory.json --output phase2_ca/outputs/phase2a3b2_exploratory_reproduced
python phase2_ca/run_phase2a3b2_diagnosis.py --stage controlled --config phase2_ca/configs/phase2a3b2_controlled.json --output phase2_ca/outputs/phase2a3b2_controlled_reproduced
```

See `docs/phase2a_3b2_anisotropy_diagnosis.md`. The diagnosis preserves the
Phase 2A.3b-1 failure, quantifies sub-voxel growth depth, audits the missing
probability-to-velocity closure, and tests one mean voxel layer with paired
seeds. The current site-weighted CA remains NO-GO as a physical morphology law.

Run the Phase 2A.3b-3 conservative interface-area development and locked
verification with:

```powershell
python phase2_ca/run_phase2a3b3_validation.py --config phase2_ca/configs/phase2a3b3_method_development.json --output phase2_ca/outputs/phase2a3b3_development_reproduced
python phase2_ca/run_phase2a3b3_validation.py --config phase2_ca/configs/phase2a3b3_locked_verification.json --output phase2_ca/outputs/phase2a3b3_locked_reproduced
```

See `docs/phase2a_3b3_interface_flux_closure.md`. The independent area-times-
kinetic-factor allocation, transactional patch ledger, and frozen-front
realization pass their locked gates. The GO is limited to this numerical
closure: the small evolving-grain comparison does not validate anisotropic NaCl
morphology, and multigrain or transport extensions remain out of scope.

Run the Phase 2A.2 rollback, grid-sensitivity, pore-envelope, and morphology
feasibility audit with:

    python phase2_ca/run_phase2a2_validation.py

Outputs are under outputs/phase2a2/. The evidence-bounded assessment is
docs/phase2a_2_morphology_feasibility.md.

Run the Phase 2A.2R corrected pore-connectivity benchmarks and selected-mask
reanalysis with:

    python phase2_ca/run_phase2a2r_validation.py

Outputs are under outputs/phase2a2r/. The corrected scientific assessment and
minimal grain-aware architecture are in
docs/phase2a_2r_pore_diagnostics_and_grain_design.md.
