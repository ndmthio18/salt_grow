# Implementation plan

1. Reproduce Prat (2026), Eqs. (18)-(29), with an independent finite-volume
   boundary-value solve. Verify the concentration profile, the general
   `tanh(sqrt(Da))` growth law, its large-Da limit, and spatial convergence.
2. Implement Model B with cell-integrated pore volume and dissolved-salt mass
   on a top-refined ALE mesh. Use the original moving-interface solution and
   salt balances, the independent kinetic interface law, and a prescribed
   evaporation flux. Assign newly created volume `epsilon_birth = 0.43`.
3. Run frozen-coefficient (B1) and porosity-dependent (B2) closures. Audit
   total salt and total solution-plus-solid mass from boundary fluxes, while
   retaining separate top-growth, internal-precipitation, dissolved, bottom
   transfer, and evaporation inventories.
4. Establish a small serial baseline, then refine grid and maximum time step
   independently. Only after the conservation and refinement gates pass, run
   the smallest sensitivity matrix requested in `instruction/1.md`.
5. Save run-level CSV data and PNG/PDF figures, then write evidence-bounded
   `RESULTS.md` and `GO_NO_GO.md`. Numerical agreement will be labelled
   verification, not physical validation.

## Baseline acceptance gates

- Stage A relative growth-velocity error: below `1e-5` on the accepted mesh.
- Model B normalized global salt and bulk-mass residuals: below `1e-6`.
- Reported height and terminal velocity: less than 1% change under the accepted
  grid and time-step refinements.
- At least four finite-volume cells across the instantaneous analytical
  concentration scale `delta_C` near the moving interface.
- Stop rather than clip if porosity or concentration leaves physical bounds.

## Initial bounded run plan

- Serial Python/SciPy BDF integration; no MPI or GPU.
- Baseline duration: 6 h; output interval: 5 min.
- Accepted baseline mesh after refinement: 160 ALE cells with
  `grid_power = 2.5`. It resolves the final `delta_C` with seven cells and
  differs from the more strongly clustered `power = 3.0` result by less than
  1% in terminal height and velocity. The quadratic and cubic meshes are
  retained in the convergence records.
- Refinement checks: 80, 120, and 160 cells; maximum BDF steps 120, 60, and
  30 s where needed.
- Expected resources: under 2 CPU cores, 4 GB RAM, 1 GB output, and 30 min wall
  time for the initial campaign. These are hard campaign bounds, not estimates
  of numerical correctness.
