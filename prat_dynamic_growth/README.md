# Dynamic internal precipitation versus Prat (2026)

This directory contains a reproducible one-dimensional prototype comparing
Prat's analytical growth law with a conservative transient ALE finite-volume
model in which porosity evolves through internal precipitation.

The implementation uses the moving-interface solution and salt balances from
Prat (2026), Eqs. (6)-(7), and computes the interface velocity independently
from the kinetic relation in Eq. (20). It never substitutes the analytical
growth velocity into Model B.

## Source and scope

- Paper: Marc Prat, *The Evaporation-Driven Growth Law of Porous Salt
  Structures*, Langmuir 42 (2026), 21269-21282,
  DOI `10.1021/acs.langmuir.6c02211`.
- Local source: `C:/Users/12608/Downloads/Prat2026.pdf`.
- Recorded SHA-256:
  `82E4F450481F7185AFAF3F1B1C2A7A60C25C16FF04B0A570FCB2C989EB23AB3A`.
- `rho_solution = 1200 kg/m3` is an explicit illustrative assumption, not a
  value attributed to Prat (2026).

## Reproduction on PowerShell

```powershell
cd "D:\3OutputFile\1Research\Article\23salt grow\prat_dynamic_growth"
$env:PYTHONPATH = (Join-Path (Resolve-Path '.').Path 'src')
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'

python -m prat_growth.attribution
python -m pytest -q
```

This runs the current M3 Q0--Q6/T0--T3 study, including spatial/time refinement,
startup controls, local/interface audits, and PNG/SVG figures in `outputs/m3/`.
Cache reuse requires matching configuration, code and data hashes. See
[`M3_METHODS.md`](M3_METHODS.md) for the frozen-snapshot equations and attribution
definitions, and [`RESULTS.md`](RESULTS.md) for the corrected conclusions.

The fixed-material M4 diagnostic is reproduced separately with:

```powershell
python -m prat_growth.m4_analysis
```

Its interpretation is in [`M4_MATERIAL_HISTORY.md`](M4_MATERIAL_HISTORY.md),
and its CSV/NPZ/PNG/SVG evidence is under `outputs/m4/`. It does not overwrite
the archived M3 evidence.

To run the older broad campaign with the **current corrected solver**, use a
new output folder so that historical M1 evidence remains intact:

```powershell
python -m prat_growth.campaign all --config configs/baseline.json --output outputs/rerun_current
```

The historical baseline campaign used an equilibrium concentration initial condition,
160 ALE cells, `grid_power = 2.5`, a 6 h duration, and SciPy's implicit BDF
integrator. The time step is adaptive; `max_step = 60 s` is only an upper
bound.

## Numerical formulation

- Cell-integrated pore volume and dissolved-salt mass are advanced on a
  normalized moving mesh.
- Mesh-motion fluxes are included explicitly; new volume enters at the top
  with `epsilon_birth = 0.43`.
- The liquid flux is reconstructed from liquid mass conservation and the
  reaction-induced divergence. The bottom reservoir flux is therefore an
  output, not a second imposed boundary condition.
- The top concentration is obtained from the moving-interface salt balance and
  kinetic velocity relation.
- No porosity or concentration clipping is used. Terminal events stop a run at
  stated physical bounds.
- B1 evolves porosity but freezes reference `a_v` and `D_eff`; B2 updates
  `a_v = 36(1-epsilon)/d_sub` and `D_eff = sqrt(epsilon) D_s`.

## Outputs

Current M3 results are in `outputs/m3/`: `contributions.csv`, `frozen.csv`,
`matched.csv`, `transients.csv`, `startup.csv`, `convergence.csv`,
`local_audits.csv`, `trajectory_residuals.csv`, full run arrays and
`manifest.json`. The manifest records solver warnings and the code hash.

The following are archived M1 outputs produced before the M3 correction:

- `outputs/stage_a_convergence.csv`: independent analytical-equation
  verification.
- `outputs/run_time_series.csv`: every recorded time point.
- `outputs/run_summary.csv`: terminal row of every registered run, including
  stopped runs and their event messages.
- `outputs/run_assessment.csv`: acceptance interpretation for each run.
- `outputs/runs/*.npz` and `*.json`: complete fields and provenance per run.
- `outputs/figures/`: required figures in PNG and PDF.

Interpretation is in `RESULTS.md`; limitations are in
`MODEL_ASSUMPTIONS.md`; the extension decision is in `GO_NO_GO.md`. Failed and
replaced development attempts are retained in `NUMERICAL_PROVENANCE.md`.
The M2 audit of the interface closure, controlled attribution matrix, hydraulic
feedback, and sparse-solver acceptance gates is in `M2_CRITICAL_REVIEW.md`.
