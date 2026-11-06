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

The material and liquid-accessibility arrays are separate. Salt-film supply is
represented only by six-face connectivity to a wetted seed. An empty air voxel
can grow when it is adjacent to supply-connected salt; the target itself need
not contain liquid.

Candidates are ranked by local six-face salt coordination. Ties use a pinned,
mixed-coordinate priority to avoid a systematic low-height axis bias while
remaining reproducible. This priority is an algorithmic prototype rule, not a
validated NaCl crystal-growth law.

The wetting filter can be disabled as an explicitly nonphysical control. This
prototype contains no pore-network flow, film hydraulics, vapor diffusion,
permeability, morphology-dependent evaporation, dissolution, or experimental
calibration. Numerical mass closure does not validate the growth physics.
