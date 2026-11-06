# Phase 2A.3b-0 — orientation-to-growth interface validation

## Decision

**Conditional GO for a bounded full Phase 2A.3b single-grain numerical pilot.**
The verified surface-orientation estimator can be coupled transactionally to
the sequential mass-conservative CA, provided that the grain is at least 16
voxels across, the fitting radius remains inside the tested dimensionless
range, and both predeclared edge/corner policies are retained as a numerical
sensitivity pair.

This is not a selection of a physically correct edge/corner rule. Non-zero
anisotropy produced strongly policy-dependent event locations, so a subsequent
campaign may compare the policies but may not choose one because it gives a
more plausible-looking crystal. No NaCl kinetic coefficient, facet velocity,
film-flow law, permeability, vapour transport, or PNM coupling was validated.

## Locked numerical protocol

The campaign configuration is
`phase2_ca/configs/phase2a3b0_protocol.json` (SHA-256
`97372f1d6e077ed492e6876f37d4a355c0411ca31b7ad55eebdfa82d544acbef`).
The independent numerical units were pinned RNG event sequences, not physical
or experimental replicates. Comparisons were paired by continuous seed,
voxelization, orientation, event count, and RNG seed.

The protocol fixed:

- cubic anisotropy magnitude \(|\alpha|=(3/2)\ln 2\);
- edge policies `smooth_average` and `isotropic_ambiguous`;
- minimum normal quality 0.60;
- two non-grid-aligned rotations, `(17°, 23°, 11°)` and
  `(-19°, 31°, 7°)` using \(R_zR_yR_x\);
- three numerical RNG seeds, 101, 202, and 303;
- a 160 µm continuous cube on a 10 µm grid, giving 16 cells across;
- 12 sequential deposition events per growth case; and
- a joint-rotation covariance threshold of 10% for normalized gyration
  eigenvalues.

## Candidate-site normal mapping

The CA grows into empty voxels, whereas Phase 2A.3a estimates normals at
solid/exterior interface-face centres. Each eligible candidate is therefore
mapped to every distinct face shared with an adjacent salt voxel. The verified
120 µm absolute radius was not copied to the smaller grain. Normals are
reconstructed by the weighted plane fit using a physical radius
\(r_{\rm fit}\), then grouped by candidate.

For every evaluated candidate the implementation records:

- lattice coordinate;
- all reconstructed adjacent-face normals;
- number of adjacent salt faces and valid reconstructed normals;
- `face_like`, `edge_like`, or `corner_like` digital classification;
- normalized resultant normal;
- quality, defined as valid-normal fraction times resultant length;
- anisotropy eligibility; and
- final cubic attachment weight.

The classification is geometric and digital. It is not a crystallographic
facet assignment. Two mutually perpendicular salt faces are `edge_like`;
opposed or three-or-more-face configurations are `corner_like` because they do
not define a unique convex edge normal.

The resulting audit table contains 16,704 candidate records in
`candidate_snapshot_diagnostics.csv`.

## Edge and corner policies

Two policies were defined before the campaign:

1. **Smooth reconstructed normal (`smooth_average`).** Valid adjacent-face
   normals are averaged and normalized. The candidate receives anisotropic
   weighting when quality is at least 0.60, irrespective of the digital class.
   This treats a locally rounded edge as a reconstructed smooth surface; it
   does not label it as a facet.
2. **Isotropic ambiguous fallback (`isotropic_ambiguous`).** The same
   reconstruction is retained for audit, but only `face_like` candidates may
   receive anisotropic weighting. Edge- and corner-like candidates receive
   weight one.

Across the controlled growth cases, approximately 59% of candidate-evaluations
were face-like, 22% edge-like, and 19% corner-like. Consequently,
`isotropic_ambiguous` applied anisotropy to approximately 59% of evaluations,
whereas `smooth_average` passed the 0.60 quality criterion for essentially all
evaluations.

This policy choice mattered. At \(\alpha=0\), both policies reproduced the
unmodified isotropic CA exactly for all 12 paired cases. At non-zero alpha, the
two policies had no identical event sequence. Their selected new-voxel sets
overlapped by only 13.9% for negative alpha and 16.7% for positive alpha on
average. Up to 24 final voxels differed, meaning that in some pairs none of the
12 new voxels coincided. Whole-mask Jaccard values near 0.995 are not evidence
of morphological agreement because the common approximately 4,095-voxel seed
dominates that statistic.

## Physical-scale transfer

The transfer matrix used continuous cubes of 320 and 480 µm, 16/24/32 cells
across, two rotations, and \(\lambda=r_{\rm fit}/L_g\) of 0.125 and 0.25.
The physical size and voxel size changed together so that dimensionless
resolution could be separated from absolute length.

| Cells across | Lambda | Mean p90 across sizes/rotations | Maximum p90 | Ambiguous edge/corner fraction |
|---:|---:|---:|---:|---:|
| 16 | 0.125 | 9.14° | 10.73° | 0.42 |
| 16 | 0.25 | 2.12° | 2.58° | 0.74 |
| 24 | 0.125 | 3.56° | 3.97° | 0.43 |
| 24 | 0.25 | 1.10° | 1.29° | 0.75 |
| 32 | 0.125 | 1.91° | 2.58° | 0.43 |
| 32 | 0.25 | 0.48° | 0.57° | 0.75 |

The selected campaign value, \(\lambda=0.25\), passed the 10° face-interior
gate at both physical sizes and all tested resolutions. However, its large
ambiguous fraction shows that the fitting neighbourhood reaches across
multiple facets for much of these small cubes. None of the tested radii spans
the entire grain (\(r_{\rm fit}<L_g\)), but the estimator is not a purely
single-facet local measurement near edges. The low face-interior error and the
high multi-facet exposure are separate observations and must not be conflated.

The near equality between the 320 and 480 µm results at matched \(N_g\) and
\(\lambda\) supports dimensionless transfer over this bounded matrix only. It
does not establish a universal grain-size rule.

## Corrected Phase 2A.3a resolution metadata

The original locked Phase 2A.3a CSV and its SHA-256
`6a00e3db0604c741a25972e05e52fab54e401cfb53b77dcc00af8570a8038b3c`
were preserved unchanged. A separate corrected table now requires evaluated
normals, valid fraction at least 0.90, at least 100 face-interior samples, and
p90 error at most 10° before a row may count as facet-orientation evidence.

The original four-voxel row had `eligible=True` despite an unresolved normal
estimator; its corrected flag is false. The eight-voxel row is also excluded
from facet evidence because p90 is 11.41°. Sixteen and 32 voxels pass. The
single-alignment volume values remain labelled as the original locked
diagnostic; the later translation/rotation robustness audit remains a separate
post-lock result. No threshold was retroactively described as pre-registered.

The Phase 2A.3a runner was corrected so future reruns derive eligibility from
the actual normal evidence rather than `cells_across > 2`.

## Cubic anisotropy integration

For grain-to-laboratory rotation \(R_g\), the experimental CA mode uses

\[
\mathbf n_c=R_g^T\mathbf n_{\rm lab},\qquad
q_4=n_{c,x}^4+n_{c,y}^4+n_{c,z}^4,
\]

\[
w=\exp[\alpha(q_4-1/3)].
\]

For positive alpha the `<100>/<110>/<111>` weights were
2.000/1.189/1.000; for negative alpha they were
0.500/0.841/1.000. The maximum weight discrepancy over cubic-symmetry-
equivalent orientations was \(2.22\times10^{-15}\).

All 12 alpha-zero orientation-mode cases reproduced the unmodified stochastic
CA event sequence, final mask, and committed mass exactly. Positive and
negative alpha changed every paired event sequence relative to alpha zero. The
new-voxel set overlap was only about 0–1.4% on average, demonstrating that the
anisotropy weights are active. This is numerical sensitivity, not evidence that
either sign describes NaCl growth.

The largest paired normalized-gyration difference between the two jointly
rotated physical problems was 0.00490, below the predeclared 0.10 bound. The
two voxelized continuous seeds contained 4,094 and 4,096 voxels, respectively.
This supports bounded numerical covariance of the short pilot; it is not exact
rotation covariance of a binary lattice.

## Transactionality and mass closure

The optional candidate-weight model is included in the CA transaction
snapshot. A failed selection, insertion, or mass commit restores material,
allocator, RNG, candidate cache, selected diagnostics, and performance
counters. New failure-injection tests verify this path.

All 42 campaign cases committed exactly 12 new voxels. The mass-ledger closure
error was exactly zero in every case. Four representative repeated runs with
pinned RNG seeds reproduced selected coordinates and final masks exactly.

## Normal-update performance

The same eight-event case compared three strategies:

| Strategy | Runtime | Recomputed candidates | Agreement with full |
|---|---:|---:|---|
| Full recomputation | 2.241 s | 11,420 | reference |
| Conservative local invalidation | 0.603 s | 2,114 | exact sequence, weights, and mask |
| Adjacency-only cache | 0.433 s | 1,460 | failed: different sequence and mask |

The conservative local strategy was selected. It invalidates cached candidates
within \(r_{\rm fit}+\sqrt{3}\Delta x\) of each insertion and was 3.71 times
faster than full recomputation in the recorded run without measurable weight
or morphology discrepancy. The faster adjacency-only cache was rejected; its
speed does not compensate for stale normals.

Runtime values are implementation and host observations, not portable
performance guarantees.

## Acceptance summary

| Gate | Result |
|---|---|
| Exact mass closure | PASS |
| Alpha-zero exact isotropic reproduction | PASS |
| Pinned-seed reproducibility | PASS |
| Cubic-equivalent orientation weights | PASS |
| Two-size dimensionless scale transfer at selected lambda | PASS |
| Joint-rotation covariance | PASS |
| Four-voxel erroneous eligibility corrected | PASS |
| Conservative local update equals full | PASS |

The numerical interface therefore passes Phase 2A.3b-0. The scientific
decision is narrower: **GO for a small, resolved, single-grain numerical
campaign while retaining both edge policies; NO-GO for selecting a physical
edge rule, calibrating alpha, validating NaCl morphology, or adding multi-grain
and transport physics.**

## Reproduction

From the repository root:

```powershell
python phase2_ca/run_phase2a3b0_validation.py --config phase2_ca/configs/phase2a3b0_protocol.json --output phase2_ca/outputs/phase2a3b0_reproduced
$env:PYTHONPATH=(Resolve-Path 'phase2_ca/src').Path + ';' + (Resolve-Path 'phase1_mass_budget/src').Path + ';' + (Resolve-Path 'prat_dynamic_growth/src').Path
python -m pytest -q
```

The runner refuses implicit overwrite. The committed output directory contains
the protocol manifest, corrected metadata, raw scale/performance/growth/
sensitivity tables, candidate-level diagnostics, figures, representative
masks, gate summary, and SHA-256 inventory.
