# `review/prat-pnm-audit`

- Tip: `facdfff86e171def76ffc845839911388c887254`
- Unique commits relative to Step 2A: 3
- Exported changed files: 11
- Role: implementation and physical-consistency audit.

The audit verifies that the `(1-F)` salt share is not double counted and that
the ledgers close for the implemented equations. It also finds material
problems: unsourced salt can be deposited by disconnected clusters;
`CrustField` hydraulic fields are written but not used by the supply gate;
hierarchical layers continue to use first-layer Prat quantities; two salt
diffusivities coexist; and a face declared supply-limited can continue to
evaporate and precipitate.

Scientific status: strong negative/diagnostic evidence. No production repair
or experimental validation is contained in this branch.

