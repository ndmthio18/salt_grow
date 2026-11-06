# Phase 2A.1 stochastic growth-rule and morphology validation

## Decision

**GO for the stochastic/sequential CA as a reproducible numerical test bed;
NO-GO for advancing to Phase 2B transport coupling or for claiming a robust
porous-salt morphology.** The sequential event algorithm removes the tested
30--120 s timestep-grouping dependence, and all matched cases conserve the
prescribed mass. However, final morphology remains strongly dependent on the
growth rule and random seed, while inferred salt fraction depends strongly on
the assumed analysis envelope. Enclosed internal voids are rare and tiny.

Adding PNM transport, vapor diffusion, or permeability now would increase
model complexity without resolving the dominant morphology-rule uncertainty.

## Evidence boundary

The implemented weights are relative sampling weights under a fixed Phase 1
precipitation mass budget. They decide *which* eligible sites receive already
prescribed salt mass; they do not change cumulative salt supply.

This is not the Redkov Monte Carlo process. The audited reference CA applies
attachment/detachment probabilities and diffusion to a closed particle
inventory during algorithmic iterations. Phase 2A.1 has no mobile particles,
detachment, diffusion sweep, calibrated Monte Carlo time, or fitted NaCl
attachment probability. Numerical similarity of a weight pattern would not
establish physical or code equivalence.

## Implemented algorithms

### Growth modes

- `deterministic`: preserved legacy control. Candidates are ranked by local
  six-face NaCl coordination and a deterministic mixed-coordinate tie-break.
- `stochastic`: candidates are sampled with a NumPy generator seeded once at
  CA construction. The generator state persists across physical timesteps.

The stochastic coordination classes use weights `w1`, `w2`, `w3`, and `w4`
for one, two, three, and four-or-more crystalline neighbors. All weights must be
finite and positive.

### Update modes

- `synchronous`: compute one frozen candidate/coordination list at the start of
  a physical step and select the requested batch from it.
- `sequential`: add one realized voxel at a time, then recompute wetting
  connectivity, eligibility, and coordination before the next event.

The sequential implementation is intentionally event based. When no other
within-step physics is present, its final state should depend on cumulative
deposition events and the persistent RNG sequence, not how those events are
grouped into physical timesteps.

### Transaction failure handling

`TransactionalMassToVoxelAllocator.abort` now clears an uncommitted proposal
without altering mass pools. Candidate-selection or commit exceptions cause
the CA to restore every tentatively added voxel, abort the pending proposal,
restore the pre-step RNG state, and recheck geometry/ledger invariants.
Injected-failure tests verify that no pending transaction, mass input,
random-sequence advance, or artificial salt survives rollback.

## Controlled matrix

Common rule-sensitivity cases use a `20 x 20 x 18` domain, `10 um` voxels, one
substrate layer, four wetted seeds, 600 s duration, and the verified Phase 1
RH=0.30 evaporation supply. Each complete case realizes 207 new voxels.

| Rule | Weights `(w1,w2,w3,w4)` | Update | Replicates |
|---|---|---|---:|
| Legacy deterministic | not applicable | synchronous | 1 |
| Uniform stochastic | `(1,1,1,1)` | sequential | 5 paired seeds |
| Strong high coordination | `(1,4,16,64)` | sequential | 5 paired seeds |
| Mild high coordination | `(1,1.5,2,3)` | sequential | 5 paired seeds |

Random seeds are `101, 202, 303, 404, 505`. These are simulation realizations,
not independent experimental replicates. Analyses therefore report raw values,
means, sample SDs, ranges, paired mask Jaccard indices, and voxel-difference
fractions. No inferential p-values are reported.

The timestep study uses nine wetted seeds, paired RNG seeds `101, 202, 303`,
an evaporation-area scale of 0.5, identical 600 s duration, and timesteps 30,
60, and 120 s. All paired cases realize 103 new voxels, so morphology is
compared at matched inventory.

## Numerical correctness

Across 44 campaign cases:

- maximum absolute mass-closure residual: `1.55096e-25 kg`;
- 43 cases realized every requested whole voxel;
- the deliberate one-step synchronous capacity stress realized 5 of 41
  requested voxels and retained 36 whole-voxel equivalents;
- the paired sequential stress realized all 41 voxels;
- summed recorded per-case simulation runtime in the final run: `31.40 s`;
- all original substrate voxels remained unchanged.

The capacity-stress pair is explicitly unmatched and is not used for morphology
comparison. It demonstrates that a frozen synchronous list can exhaust its
initial eligible sites, while sequential events expose new candidates inside
the same physical step.

## Timestep sensitivity

All sequential paired cases were exactly equal for 30, 60, and 120 s:

- mask Jaccard: `1.0`;
- differing-domain fraction: `0`;
- height and footprint differences: `0`.

This confirms the expected event-sequence invariance for the present model,
which has no other explicit within-step physics.

The matched-inventory synchronous cases were strongly timestep dependent:

- paired mask Jaccard relative to 30 s: `0.191--0.295`;
- differing-domain fraction: `0.0169--0.0211`;
- height differences: `-2` to `0` voxels;
- footprint differences: `-1` to `+12` voxels.

The synchronous effect is algorithmic: coarse batches sample from a frozen
front and do not expose sites created earlier in the same step.

## Growth-rule sensitivity

All rule cases contain exactly 207 newly committed voxels.

| Rule | Height, mean ± SD [vox] | Footprint, mean ± SD [vox] | Bounding-box salt fraction | Column-envelope salt fraction |
|---|---:|---:|---:|---:|
| Legacy deterministic (`n=1`) | `3` | `82` | `0.451` | `0.972` |
| Uniform (`n=5`) | `6.0 ± 1.22` | `112.0 ± 8.86` | `0.146 ± 0.026` | `0.744 ± 0.036` |
| Strong high coordination (`n=5`) | `5.8 ± 0.84` | `86.2 ± 9.93` | `0.184 ± 0.051` | `0.909 ± 0.022` |
| Mild high coordination (`n=5`) | `6.0 ± 2.00` | `105.0 ± 8.00` | `0.159 ± 0.045` | `0.824 ± 0.049` |

Strong high-coordination weighting consistently reduces footprint relative to
the paired uniform case (`-18` to `-34` voxels) and increases the local column
salt fraction (`+0.124` to `+0.253`). It does not consistently increase height:
paired differences range from `-2` to `+2` voxels. Thus the main detected
effect is compaction/local filling, not a robust vertical-growth enhancement.

Paired rule masks remain substantially different even with the same RNG seed
and mass. Strong-versus-uniform Jaccard values are `0.230--0.357`; mild-versus-
uniform values are `0.252--0.397`.

## Nucleation and connectivity

At matched 207-voxel precipitation inventory:

- single-seed cases remain one salt cluster, reach height 7 voxels, and occupy
  footprints of 79--80 voxels;
- four-seed cases reach heights 5--6 voxels and footprints 111--116 voxels;
- four-seed cases undergo two cluster coalescences but finish with two clusters
  for all three paired seeds.

In the dry-seed control with enforced wetting, cluster sizes are `[104,1]`:
the disconnected dry seed remains exactly one voxel. With wetting disabled,
the same mass produces clusters `[80,25]`, so the otherwise dry seed grows.
Only 80 of 105 salt voxels in that nonphysical control remain connected to the
actually wetted seed. This separates solid-network connectivity from assumed
liquid-film supply connectivity.

## Void morphology and envelope sensitivity

Two explicit analysis windows are retained:

1. the tight axis-aligned bounding box of all salt;
2. footprint columns from the substrate surface to each column's local salt top.

The bounding box includes large regions of exterior air between separated
clusters; the column envelope excludes those unsupported lateral gaps. Their
salt fractions therefore differ strongly. For example, the stochastic rule
means span `0.146--0.184` in the bounding box but `0.744--0.909` in the column
envelope. Reporting `1 - salt fraction` as one intrinsic porosity would be
misleading.

Most void voxels within the bounding box remain connected to the external air
boundary. Enclosed voids are sparse:

- uniform: mean `1.8`, range `0--4` voxels;
- mild high coordination: mean `1.0`, range `0--3` voxels;
- strong high coordination and deterministic control: `0` enclosed voxels.

The generated structures range from compact crystalline clusters to open,
branched aggregates, but they do not exhibit a robust, internally connected
porous-salt void network at this size and duration. “Potentially porous” is
therefore only a geometric possibility under the loose bounding-box envelope,
not a supported intrinsic-material conclusion.

## Figures and underlying data

- `rule_morphologies.png`: one 3D rendering plus horizontal and vertical
  sections for each rule at matched supplied mass.
- `rule_metric_distributions.png`: all stochastic seed observations with mean
  ± sample SD; the deterministic control is a single explicitly marked point.
- `timestep_synchronous.png` and `timestep_sequential.png`: paired seed 101 at
  30, 60, and 120 s.
- `nucleation_morphologies.png`: single- versus four-seed structures.
- `wetting_control.png`: enforced-wetting and nonphysical bypass comparison.
- `capacity_stress.png`: unmatched one-step failure/control pair.

Every morphology row includes a true 3D voxel rendering and two orthogonal
sections. Slice indices are recorded in `figure_slice_indices.csv`. Salt is
blue and initial seeds orange; rule groups also use marker-shape redundancy.
PNG files are 200 dpi, RGB, opaque white. Palette screening found adequate
contrast against white; several color pairs have weak grayscale separation,
so colors are never used as the sole encoding.

Machine-readable evidence includes:

- `cases.csv` and `conservation_history.csv`;
- rule descriptive and paired-difference tables;
- paired timestep comparisons;
- 44 final 3D Boolean `.npy` masks;
- representative initial/mid/final masks;
- configuration/provenance manifest and SHA-256 inventory.

## Tests executed

```text
python -m pytest phase2_ca/tests -q
30 passed

python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests phase2_ca/tests -q
99 passed
```

Tests cover stochastic reproducibility, weight-dependent geometry at fixed
mass, sequential timestep partition invariance, synchronous timestep/capacity
sensitivity, wetting connectivity, direct geometry/ledger agreement,
insufficient capacity, component and void classification, and injected
selection/commit rollback failures.

## Remaining scientific limitations

- Coordination weights are uncalibrated and omit crystallographic orientation,
  attachment/detachment kinetics, mobile particles, and surface energy.
- Five stochastic seeds quantify only small-ensemble algorithmic spread.
- Domain, duration, substrate geometry, and seed placement remain synthetic.
- The continuous-film rule is connectivity logic, not a liquid-film flow model.
- The column envelope is an analysis assumption and not a measured salt-body
  boundary.
- No comparison with experimental 3D NaCl morphology or a licensed executable
  reference CA is available.

## Phase 2B recommendation

**NO-GO at present.** Advance only after at least one of the following gates is
met: morphology-rule calibration against independent 3D evidence; robustness
to domain/seed/envelope choices; or a justified physical attachment model that
reduces the current algorithmic degrees of freedom. Until then, transport
coupling would produce precise fields on an underdetermined geometry.
