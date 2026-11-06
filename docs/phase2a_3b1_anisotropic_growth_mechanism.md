# Phase 2A.3b-1 anisotropic-growth mechanism and morphology verification

## Decision

**NO-GO for a complete single-grain anisotropic morphology model.** The
orientation-dependent selector is numerically correct and produces a measurable
event-level statistical response, but the locked 10% growth campaign does not
show a growth-only morphology difference beyond the predeclared alpha-zero grid
and stochastic baseline. This result does not support multigrain growth or
transport coupling.

The bounded conclusions are:

- **GO for numerical mechanism correctness:** normalized fixed-front sampling,
  conservative local cache updates, mass accounting, and the separate anchored
  supply control pass their locked gates.
- **GO for an event-level statistical response:** changing alpha changes the
  crystallographic-family distribution of selected new voxels.
- **NO-GO for morphology verification:** the event-level response does not
  produce a locked growth-only shape signal above the registered baseline.
- **NO physical validation:** alpha is uncalibrated, edge policy is not selected,
  and neither facet velocities nor liquid-film hydraulics are solved.

## Scope, provenance, and locked protocol

The work started from commit
`3a9a7b2575f0c1cb5e0eef2d540ed26257e2a52e` on the dedicated branch
`codex/phase2a3b1-anisotropic-growth`. Method development and locked verification
were separated:

| Item | Value |
|---|---:|
| Development configuration SHA-256 | `c9cc1d9abf19c98ae3c21de2c3347ec148e7b8e1b663c0c470e9c5bb7d6beeee` |
| Development gate SHA-256 | `5462eaa3da893a2696e8698f1cc74c9c7242cdc346e6fd65da217afdfb17c80f` |
| Locked configuration SHA-256 | `162dabd1956c04b3f9a37ea825235cfc3df5e7fe51f359ad8406aabff268c3a3` |
| Voxel size / nominal grain size | 10 / 160 micrometres |
| Normal-fit radius | 0.25 grain size |
| Alpha | 0 and plus/minus 1.03972077084, i.e. plus/minus (3/2) ln(2) |
| Edge policies | `smooth_average`, `isotropic_ambiguous` |
| Non-grid-aligned rotations | (17, 23, 11) and (-19, 31, 7) degrees |
| Paired random seeds | 101, 202, 303 |
| Growth stages | 41, 205, 410 new voxels, approximately 1%, 5%, 10% |

Development required 100-event full/local identity and fixed-front probability
agreement before locking. It also covered face-like, edge-like, a true
three-face corner-like insertion, and a domain-boundary insertion. All passed.
The locked thresholds were then frozen at absolute probability error at most
0.005, standardized probability error at most 4.5, cache weight error at most
1e-12, nonzero-alpha rotation excess at most 1.5 voxels, and morphology
signal/baseline ratio at least 1.25.

The first locked execution exposed a diagnostic bug: selected candidates with
an unresolved normal were passed to a routine that required finite nonzero
directions. Such events are isotropic by construction and are now retained in
the explicit `other` class with an undefined angular error. A regression test
and a 41-event end-to-end check were added. No threshold or controlled
configuration changed, and the locked run was restarted from the beginning.

## Fixed-front probability verification

Frozen candidate sets were built from a digital sphere and the rotated cube.
They contained candidates near `<100>`, `<110>`, and `<111>` as well as an
explicit `other` bin. For every candidate the output records the reconstructed
lab normal, crystal-frame normal, q4, weight, represented face area, expected
normalized probability, and empirical count.

There were eight independent fixed-front distributions: two geometries, two
edge policies, and two nonzero alpha values. Each was sampled 100,000 times
without evolving the front. Across all family bins:

| Locked statistic | Result | Threshold |
|---|---:|---:|
| Maximum absolute probability error | 0.002227 | 0.005 |
| Maximum absolute standardized error | 2.073 | 4.5 |
| Probability gate | PASS | both limits |

The weights are distinguishable before any front feedback. For example, on the
sphere with `smooth_average` and positive alpha, the family-mean weights were
1.868 (`<100>`), 1.220 (`<110>`), and 1.060 (`<111>`). Negative alpha reversed
the preference, with means 0.536, 0.822, and 0.944. Empirical frequencies agree
with the probabilities obtained directly by normalizing all candidate weights.
These are algorithmic selection probabilities on a frozen digital surface;
they are **not** physical facet growth-velocity ratios.

## Face, edge, corner, and refinement attribution

The rotated cube was digitized at 16, 24, and 32 cells across the nominal grain.
The fraction of candidates assigned to face/edge/corner classes did not converge
toward a purely face-interior population:

| Cells across | Face-like | Edge-like | Corner-like |
|---:|---:|---:|---:|
| 16 | 0.638 | 0.226 | 0.136 |
| 24 | 0.633 | 0.217 | 0.150 |
| 32 | 0.622 | 0.223 | 0.156 |

At 32 cells across, face-like candidates contributed 63.0-63.4% of the total
probability shift under `smooth_average`; edge-like sites contributed about
22.2-22.3% and corner-like sites 14.4-14.7%. Under
`isotropic_ambiguous`, face-like candidates contributed 54.4-57.8% and the
combined edge/corner share was 42.2-45.6%. In the latter policy, ambiguous
sites retain unit local weight, but their normalized selection shares still
change indirectly when face weights change.

Thus, the selector contains a genuine face-interior crystallographic component,
but ambiguous digital geometry remains material and does not diminish over this
16-32-cell refinement range. The evidence does not justify calling the observed
response a grid-independent continuum facet competition, and neither edge
policy is preferred on aesthetic grounds.

## Longer single-grain growth

The locked campaign contains 36 trajectories: two rotations by two policies by
three alpha values by three paired seeds. The initial rotated cubes contained
4094 or 4096 voxels. Cumulative deposited mass produced 41, 205, and 410 new
voxels; the final realized increment was 10.010-10.015% of the initial seed.
All primary comparisons use newly deposited voxels in the crystal frame rather
than whole-crystal overlap or gyration.

At 10% growth, median new-event fractions showed a clear statistical response:

| Policy | Alpha | `<100>` | `<110>` | `<111>` | Other |
|---|---:|---:|---:|---:|---:|
| smooth | negative | 0.512 | 0.156 | 0.017 | 0.329 |
| smooth | zero | 0.602 | 0.134 | 0.010 | 0.259 |
| smooth | positive | 0.673 | 0.098 | 0.007 | 0.227 |
| isotropic ambiguous | negative | 0.556 | 0.132 | 0.010 | 0.307 |
| isotropic ambiguous | zero | 0.602 | 0.134 | 0.010 | 0.259 |
| isotropic ambiguous | positive | 0.602 | 0.127 | 0.015 | 0.259 |

The corresponding median surface-area fractions at 10% remained much closer:
`<100>` ranged from 0.570 to 0.595, `<110>` from 0.108 to 0.138, `<111>` from
0.0076 to 0.0099, and `other` from 0.273 to 0.293. Median final
surface-to-volume ratio ranged from 6.87e4 to 7.59e4 per metre. Median total
growth extents across the crystal axes ranged from 3.49 to 5.65 voxels,
depending on condition. Full tables also report area-normalized event counts,
air-exposed area, facet-area distribution, directional positive/negative
extents, and per-stage runtime.

The important evidence split is therefore:

1. alpha changes which orientation classes receive events;
2. that response has not produced a morphology displacement larger than the
   alpha-zero grid/stochastic displacement by 10% added volume.

## Growth-only morphology signal and rotation covariance

The primary signal is the symmetric Chamfer distance between cumulative new
voxel locations for nonzero alpha and the paired alpha-zero case, after
transforming to the crystal frame. The conservative baseline was the larger of
the median alpha-zero rotation discrepancy and median alpha-zero stochastic
discrepancy: 1.1258 voxels.

| Policy | Alpha sign | Median signal (voxels) | Signal / baseline |
|---|---|---:|---:|
| smooth average | negative | 1.0782 | 0.958 |
| smooth average | positive | 1.0720 | 0.952 |
| isotropic ambiguous | negative | 1.0049 | 0.893 |
| isotropic ambiguous | positive | 1.0719 | 0.952 |

All four ratios fail the locked minimum of 1.25. The failure is not evidence
that alpha has exactly zero effect; it means the tested campaign cannot separate
its morphology effect from the registered numerical/stochastic background.

Joint-rotation covariance itself passed. Across all stages, the largest
nonzero-alpha Chamfer excess over the paired alpha-zero rotation error was
0.487 voxels, below the 1.5-voxel limit. At the final stage the excess range was
-0.074 to 0.106 voxels. This bounds orientation implementation error but does
not rescue the failed morphology-signal gate.

## Conservative local cache verification

Eight locked 256-event comparisons covered two rotations, two policies, and two
random seeds. Every selected sequence and final mask matched full recomputation;
the maximum selected-weight error was zero. Periodic full checks at every stage
of all 36 growth trajectories also found zero weight error, zero site-class
mismatches, and zero eligibility mismatches. Targeted corner, concave, and
domain-boundary insertions passed.

The local method gave a median 7.76-fold runtime speedup over full recomputation
(range 7.19-8.09) while full recomputation evaluated a median 14.36 times as many
candidates. This validates the cache numerically for the tested single-grain
matrix only.

## Mass conservation, runtime, and anchored supply control

The 108 recorded free-standing stages had a maximum absolute mass-closure error
of 1.03e-25 kg. Median and maximum cumulative runtime per 10% trajectory were
28.39 and 63.46 s, respectively. The complete locked workflow used six CPU
workers and took 476.77 s on the recorded Windows/Python environment.

The physically distinct control shifted the rotated crystal to touch the
one-layer substrate. It had four substrate-contact voxels and one explicit
liquid-supply anchor at `[19, 12, 1]`. Phase 1 conditions were 298.15 K,
water activity 0.753, RH 0.30, mass-transfer coefficient 0.005 m/s,
evaporation area 2.56e-8 m2, and duration 1800 s. The evaporation module
delivered 8.62059e-10 kg NaCl; 399 voxels realized 8.61840e-10 kg and
2.18670e-13 kg remained deferred, with zero ledger closure error. Initial and
final air-exposed face areas were 2.132e-7 and 3.135e-7 m2.

This is a prescribed single-anchor supply control. It does not establish that
all crystal surfaces remain wetted and it does not solve liquid-film hydraulics.

## Reproducibility

Run the development gate first, then the immutable locked configuration:

```powershell
python phase2_ca/run_phase2a3b1_validation.py --stage development --config phase2_ca/configs/phase2a3b1_method_development.json --output phase2_ca/outputs/phase2a3b1_development_reproduced
python phase2_ca/run_phase2a3b1_validation.py --stage locked --config phase2_ca/configs/phase2a3b1_locked_verification.json --output phase2_ca/outputs/phase2a3b1_locked_reproduced
```

Raw candidate/event tables, stage summaries, masks, figures, manifests, and
SHA-256 inventories are stored under
`phase2_ca/outputs/phase2a3b1_development/` and
`phase2_ca/outputs/phase2a3b1_locked/`. The locked manifest explicitly sets
`physical_validation_claimed` to false.

## Final recommendation

Retain the verified selector and conservative local cache as research code, but
do not label the present result a verified anisotropic crystal-morphology model.
A future study would require an independently motivated extension in growth
volume, resolution, or experimental calibration with a newly locked protocol;
the failed 1.25 signal gate must not be relabelled after inspection. Phase
2A.3b-1 stops here. Multigrain growth, impingement, film flow, permeability, and
PNM coupling remain outside the authorized and evidential scope.
