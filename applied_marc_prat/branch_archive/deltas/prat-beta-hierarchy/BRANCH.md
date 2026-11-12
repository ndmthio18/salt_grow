# `review/prat-beta-hierarchy`

- Tip: `a43a0a80ca07c6ce140432a3c2d6521fe507c2f2`
- Unique commits relative to Step 2A: 2
- Exported changed files: 49
- Role: beta geometry, hierarchy sensitivity, and solver audit.
- Local restore test on 2026-11-12: `11 passed in 87.76 s`.

The route contains the layered-validation branch plus beta/hierarchy sweeps.
It shows that `beta = 2` follows from Prat's particular free-cube growth cell,
not from a universal crystal constant; contact, partial wetting, habit, or
face-dependent kinetics change beta. The predicted closure layer shifts with
beta and the imposed hierarchy.

The solver audit also finds that testing closure on every Picard iterate can
close one layer early at large time steps. Testing only the converged step
removes that time-step-dependent layer jump in the audited cases.

Scientific status: useful numerical and causal audit; not a production merge
or an experimental identification of beta/hierarchy.

