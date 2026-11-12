# liqin Marc–Prat work: migration and scientific assessment

## Bottom line

The liqin work contains two genuinely useful assets:

1. a well-tested analytical reproduction of Prat (2026), including explicit
   parameter provenance and an independent one-dimensional reaction–transport
   check; and
2. a modular applied framework connecting vapour transport, pore-network
   supply, salt growth, and conservation ledgers in 2-D/3-D.

The first asset is suitable as a formula-level reference. The second is useful
as an implementation library and hypothesis generator, but its morphology and
patchy–crusty predictions do not meet a top-journal validation standard. Those
predictions depend on additional closures that are neither identified by the
Prat data nor independently validated.

This snapshot therefore complements, rather than replaces, the current local
`prat_dynamic_growth` model. The local model contains the later conservative
ALE treatment, controlled attribution, and fixed-material histories; the
liqin snapshot contributes the analytical figure reproduction, gas/PNM/3-D
coupling structure, and the earlier surface-first verification ladder.

## Evidence that is strong enough to reuse

### Analytical law and numerical reference

`runs/S_surface/step2a_prat/` verifies the large-Da Prat growth law against an
independent moving-coordinate finite-volume reference. In the stated case,
the reported growth velocity, front/internal precipitation partition,
concentration profile, and porosity closure converge to the analytical result
at roughly `10^-5` relative level. Salt and solid residuals pass the registered
acceptance gates.

The verified statement is narrow: under saturated supply, constant density,
one-dimensional geometry, Prat's linearisation, and fixed transport properties
in the reaction band, equations (1)–(3), (6), (7), and (20) produce the
analytical relations used in equations (29), (31), and (32).

`prat2026/` is also useful because it records otherwise hidden or inferred
inputs: solution density `1200 kg/m3` is back-calibrated, the Fig. 5 height is
inferred, and the Fig. 4 `Da*` curve appears more consistent with `beta` near
1.8 than the stated value 2. These are reproducibility facts, not independent
physical validation.

### Software structure and conservation accounting

`saltv2/` separates materials, PNM transport, growth, and coupling layers and
contains explicit liquid, salt, and solid ledgers. Its tests cover interface
conservation, pore precipitation, vapour transport, geometry, and 2-D/3-D
drivers. This architecture is worth mining when integrating a local growth
law with a pore-network or gas-phase solver.

## Evidence–claim mismatches that must be preserved

### One-dimensional verification does not validate the local 3-D law

Prat's experimental comparisons use one-dimensional or area-averaged
observables. The local normal-growth extension used by `saltv2` has no direct
local experimental validation. A successful Step 2A reference test therefore
cannot support a claim that the model predicts 3-D morphology.

### `F`, porosity, and kinetic parameters are not separately identified

The height and salt-mass observables constrain the combination of growth
partition and solid fraction, not `F`, `k_r`, `a_v`, `beta`, and newborn
porosity separately. Treating `epsilon = 0.43` as a newborn-layer porosity is
an extra model decision. If the measured 0.43 instead represents the final
salt structure, the inferred particle-size dependence changes materially.

### Patchy–crusty behaviour comes from added closures

The transition is governed mainly by choices outside the verified Prat law:

- how birth and final porosity are related;
- whether `Da*` uses the substrate scale or the current hierarchical layer;
- whether the sub-micrometre reaction length `delta_C` is treated as a
  physical sealing thickness;
- the regularisation length and nucleation placement;
- the permeability/capillary gate and its calibrated factor; and
- the absence of a dissolution branch when supplied brine is undersaturated.

These assumptions can generate plausible trends, but they do not constitute
an identified mechanism. The model may reproduce a pattern for the wrong
causal reason.

### Constant-density closure is internally conservative but physically idealised

The Step 2A equations remove a brine volume `S/rho` during internal
precipitation. The associated water fraction would physically remain as more
dilute liquid. The reference solver correctly closes the equations it was
given; that does not remove the constitutive idealisation.

### Continuum reaction-band closure is below its likely representative scale

The inferred reaction thickness is often only a small fraction of the crystal
or pore-throat scale. Treating that band as a homogeneous continuum sealing
layer is therefore the weakest part of the crust mechanism and requires a
separate pore/crystal-scale validation or a sensitivity envelope.

## Recommended use in the current project

- Reuse `prat2026/core.py`, its formula tests, and the Step 2A independent
  reference as cross-checks for the current analytical/ALE implementation.
- Reuse the `saltv2` module boundaries and conservation-ledger patterns when a
  gas or PNM coupling is needed.
- Do not merge historical porosity, hierarchy, sealing, or hydraulic-gate
  closures without an explicit model-choice table and sensitivity study.
- Do not treat stored `.npz` equality as portable reproducibility. The current
  regression instrument intentionally requires bit-for-bit identity and fails
  on a different Python/NumPy/platform at floating-point roundoff scale.
- Before a publication claim about morphology or transition, require data that
  jointly constrain salt height, salt mass, and porosity, plus a validation of
  the chosen sealing and hydraulic length scales.

## Deliberately excluded from the migration

Git metadata, virtual environments, caches, paper PDFs, Claude/Codex session
state, and large E1 HTML/image reports were not copied. They add storage or
execution-history noise without strengthening the reusable source or its
scientific audit trail.

