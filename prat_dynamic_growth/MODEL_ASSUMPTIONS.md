# Model assumptions and limitations

## Paper-grounded elements

The prototype follows Prat (2026) for the volume-averaged dissolved-salt,
solution-mass, and porosity equations; the moving-interface solution and salt
balances; `beta = 2` for the cubic-crystal construction; the reference
`k_r = 2.3e-3 m/s`; and the definitions of `Pe`, `Da`, `Da_top`, and `Da*`.

## Added assumptions

1. The solution density is constant at `1200 kg/m3`. This is an illustrative
   saturated-brine value supplied by the task, not a paper measurement.
2. Vapor density in the top solution balance is neglected, as allowed in the
   task and justified by its smallness relative to liquid and crystal density.
3. The initial concentration is uniform at `C_eq`; the first few minutes are a
   numerical/physical startup transient, not the paper's quasi-steady profile.
4. New volume enters with `epsilon_birth = 0.43`. This birth porosity is used in
   the top moving-boundary inventory, while the adjacent existing top cell can
   subsequently densify through internal precipitation.
5. B2's `a_v = 36(1-epsilon)/d_sub` is an exploratory closure inherited from a
   cube-packing estimate. It does not represent crystal coarsening, changing
   morphology, or surface shielding.
6. `D_eff = sqrt(epsilon) D_s` changes with porosity in B2, but no dispersion,
   activity-coefficient, concentration-dependent molecular diffusion, or
   anisotropy is included.
7. The finite-volume advective flux is first-order upwind. Baseline transport is
   diffusion/reaction dominated (`Pe << Da`), but numerical diffusion remains a
   discretization limitation.
8. The prescribed evaporation flux remains identical in all models. There is no
   permeability law, capillary pressure, liquid disconnection, gas transport,
   or evaporation-area feedback.
9. A unit cross-sectional area and one-dimensional geometry are assumed. Patchy
   coverage, lateral growth, branching, mechanics, and detachment are absent.
10. Dissolution is allowed by the signed linear reaction law. No arbitrary
    reaction or porosity clipping is applied. Runs stop at `epsilon = 0.02`
    because the connected-continuum assumptions are no longer credible.

## Conservation accounting

The salt inventory is

`integral[rho epsilon C + rho_crystal (1-epsilon)] dz`.

Its change is compared with the full bottom advective-plus-diffusive salt flux.
The total solution-plus-solid inventory is compared with bottom solution influx
minus evaporation. Separate integrated states record new-layer solid, internal
solid, bottom salt transfer, bottom solution transfer, and evaporation.

## Continuum limitation

For the accepted baseline, `delta_C = 0.511 micrometres` and
`L95 = 1.24 micrometres`, whereas the paper's crystal estimate gives
`d_c approximately 0.14 d_sub = 3.5 micrometres`. Thus the modeled active layer
is smaller than one representative crystal. Grid convergence verifies the
discrete continuum equations; it does not validate their physical applicability
at this scale.

## Meaning of the sensitivity stops

The `k_r = 2.3e-6` and `2.3e-5 m/s` runs reach the porosity stop within minutes.
This is not evidence that real salt structures necessarily clog in that manner.
It is evidence that, under fixed evaporation and the selected closures, the
initial model becomes internally inconsistent before 6 h. A permeability and
capillary-supply feedback would be required to continue physically.
