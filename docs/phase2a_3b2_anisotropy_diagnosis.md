# Phase 2A.3b-2 weak-anisotropy mechanism diagnosis

## Decision

**NO-GO for continuing the porous-salt model with the current site-weighted
growth law.** Phase 2A.3b-1 remains a failed morphology verification. The new
study shows that its 10% campaign had sub-voxel mean facet advancement, but this
is not the only cause of the weak morphology signal. Candidate-site density,
digital face geometry, edge policy, stochastic variation, and the density
dependence of symmetric Chamfer distance also materially affect the result.

The current law is a valid numerical candidate selector. It is not a physically
interpretable local facet-normal velocity law. A separately specified and
verified area/flux-based kinetic closure is required before multigrain or
porous-salt continuation.

## Provenance and evidence boundaries

This work starts from the immutable Phase 2A.3b-1 commit
`f39fb8bb1e304f50c5b0662cbb5e100f54dd19c4`. Its NO-GO result and 1.25
morphology signal/baseline threshold were preserved without reinterpretation or
retuning.

Exploration and controlled calculation were separated:

| Item | SHA-256 |
|---|---|
| Exploratory configuration | `97c2231df5b67fbf91767e501659bb942b2643ae0f4f94bb4346353ac8fca038` |
| Exploratory gate | `7a52929216d3595250ddd5e58d09c5028a254901c4e10025e1574ef077108dbf` |
| Exploratory inventory | `6668b6d059e66479c561224c861cbe7e798eb17de3e06640ea4c0a17c20c1b19` |
| Controlled configuration | `aadab2f5c6ff4c089785a467d75af0631f930f4f139b766c75edbde4d697a33e` |

The exploratory stage reused the existing Phase 2A.3b-1 events, audited frozen
planes at 5 and 10 micrometres, and timed one 42.38% growth trajectory. The
single feasibility case realized all 1735 requested events, remained off every
domain boundary, passed periodic full/local cache checks and mass closure, and
finished in 145.25 s. Only then was the controlled 15-case protocol frozen.

The controlled study uses one non-grid-aligned continuous 160 micrometre cube,
10 micrometre voxels, a 48 cubed domain, alpha equal to zero and plus/minus
(3/2) ln(2), both edge policies, three paired seeds, and cumulative 1%, 5%,
10%, 20%, and 42.3828% growth stages. Alpha zero is mathematically independent
of edge policy and was therefore run once per seed as the common paired control.
This is a one-resolution, one-orientation mechanism diagnosis, not a convergence
or physical-validation study.

## 1. Geometric growth resolution

For uniform expansion of a cube with initial side length L0, the one-sided
facet displacement is

\[
\Delta h=\frac{L_0}{2}\left[(1+f)^{1/3}-1\right].
\]

For L0 = 160 micrometres and a 10 micrometre voxel:

| Added solid volume | Displacement (micrometres) | Displacement (voxels) |
|---:|---:|---:|
| 1% | 0.266 | 0.0266 |
| 5% | 1.312 | 0.131 |
| 10% | 2.582 | 0.258 |
| 20% | 5.013 | 0.501 |
| 42% | 9.919 | 0.992 |
| 42.3828% | 10.000 | 1.000 |

Thus the previous 10% campaign supplied only about one quarter of an average
voxel layer under uniform-face expansion. It was intrinsically weak for
resolving an average interface displacement. This calculation does not imply
that every CA surface point moves by that amount: stochastic voxel insertion is
patchy and can create multi-voxel protrusions while the area-averaged advance
remains sub-voxel.

## 2. Existing Phase 2A.3b-1 data before new calculation

The saved 14,760 event records were analyzed without rerunning them. Paired
comparisons separated alpha response, alpha-zero seed variation, edge policy,
and joint rotation/grid effects. At 10% growth:

| Effect source | n | Chamfer q25 / median / q75 (voxels) | Event-family L1 median | Facet-area L1 median |
|---|---:|---:|---:|---:|
| Alpha response | 24 | 1.028 / 1.066 / 1.100 | 0.115 | 0.0337 |
| Alpha-zero seed | 12 | 1.085 / 1.098 / 1.111 | 0.034 | 0.0286 |
| Edge policy | 18 | 0.000 / 1.029 / 1.072 | 0.085 | 0.0289 |
| Rotation/grid | 18 | 1.126 / 1.150 / 1.173 | 0.078 | 0.0320 |

The alpha-response Chamfer distribution lies inside the seed and rotation/grid
background. Event-family response is more separated, but exposed facet-area
change is not. Across all 216 paired comparisons, Spearman correlation between
Chamfer and event-family L1 was 0.643, whereas correlation between Chamfer and
facet-area L1 was -0.120. Chamfer is therefore a complementary spatial
dissimilarity measure, not a signed measure of facet-normal evolution.

## 3. Probability-to-velocity audit

For a frozen facet family k with Nk candidate sites, physical area Ak, and mean
weight wk, the implemented law gives

\[
P_k=\frac{N_k\bar w_k}{\sum_j N_j\bar w_j},\qquad
\Delta h_k=\frac{\Delta V P_k}{A_k}.
\]

Consequently,

\[
\frac{\Delta h_k}{\Delta h_l}
=\frac{(N_k/A_k)}{(N_l/A_l)}\frac{\bar w_k}{\bar w_l}.
\]

The attachment-weight ratio equals a normal-displacement ratio only if candidate
site densities are equal. They are not.

For the 10 micrometre, 320 micrometre-wide frozen planes:

| Family | Candidates | Candidates per physical area (m-2) | Digital / physical area |
|---|---:|---:|---:|
| `<100>` | 1024 | 1.000e10 | 1.000 |
| `<110>` | 992 | 6.850e9 | 1.370 |
| `<111>` | 768 | 5.774e9 | 1.696 |

At 5 micrometres the absolute densities increase fourfold, while the ratios to
`<100>` are 0.696 for `<110>` and 0.577 for `<111>`. The orientation-dependent
site-density bias therefore persists with refinement over this audit.

Even at alpha zero, equal weights produce effective site-weighted displacement
ratios 1 : 0.685 : 0.577 at 10 micrometres. For smooth positive alpha, the
implemented ratios are 1 : 0.407 : 0.289, whereas grouping flux by physical
area would give 1 : 0.595 : 0.500, matching the family mean-weight ratios. The
hypothetical area-aware numbers are an audit, not an implemented model.

The `isotropic_ambiguous` policy exposes a more severe interpretation problem.
The `<110>` digital plane is entirely edge-like and the `<111>` plane is 93.75%
corner-like at 10 micrometres; both therefore receive unit site weight rather
than the intended crystallographic factor. A policy intended as numerical
fallback changes the effective kinetic law.

### Required closure before implementation

A physically interpretable alternative would assign each interface patch or
facet a volume flux

\[
\frac{dV_k}{dt}\propto A_k\,\kappa(\mathbf n_k),
\]

normalize the integrated flux to the externally supplied salt-volume budget,
and use voxel insertion only as a conservative geometric realization. Edges and
corners require an explicit surface partition or flux-sharing rule. This closure
was documented before the controlled calculation and was deliberately not
implemented in Phase 2A.3b-2.

## 4. Controlled one-layer-depth campaign

The controlled matrix contains 15 unique trajectories and 75 stage states,
ending at 1735 new voxels, or approximately 42.38% of the initial 4094-voxel
crystal. Across all states:

- maximum absolute mass-closure error: 1.03e-25 kg;
- maximum full/local cache weight error: zero;
- site-class and eligibility mismatches: zero;
- salt voxels on the domain boundary: zero;
- median final per-case runtime: 185.30 s;
- final runtime range: 163.88-221.63 s;
- complete six-worker controlled workflow: 593.79 s.

At final depth, median equivalent displacements by family were:

| Policy | Alpha | `<100>` (um) | `<110>` (um) | `<111>` (um) |
|---|---:|---:|---:|---:|
| Smooth | negative | 3.763 | 5.911 | 5.614 |
| Smooth | zero | 4.295 | 4.008 | 3.143 |
| Smooth | positive | 4.943 | 3.209 | 2.520 |
| Isotropic ambiguous | negative | 4.530 | 5.261 | 5.333 |
| Isotropic ambiguous | positive | 4.258 | 3.544 | 2.695 |

The signed `<100>-<111>` equivalent-displacement contrast relative to paired
alpha zero had the expected sign in all three negative-alpha seeds for both
policies. Positive alpha had the expected sign in only two of three seeds for
each policy. Final paired contrast values (um) were:

| Policy | Negative alpha | Positive alpha |
|---|---|---|
| Smooth | -3.735, -1.928, -1.972 | -0.225, 2.677, 1.413 |
| Isotropic ambiguous | -2.788, -0.249, -2.008 | 0.133, 2.400, -1.202 |

This is a systematic orientation response for negative alpha and a weaker,
seed-sensitive response for positive alpha. Its magnitude is policy-dependent,
so it is not robust enough to select a physical edge rule.

Surface-area fractions changed much less cleanly. At final depth the median
`<100>` fraction change ranged from -0.023 to -0.094 across conditions, while
the `other` fraction increased by 0.064-0.098. Digital non-ideal surface area
therefore grows materially even when the selected-event distribution is
crystallographically biased.

## 5. Why the Chamfer gate remains failed

The final alpha-zero seed-pair Chamfer values were 0.744, 0.771, and 0.774
voxels (median 0.771). Nonzero-alpha results were:

| Policy | Alpha | Paired Chamfer range | Median | Signal / seed baseline |
|---|---:|---:|---:|---:|
| Smooth | negative | 0.729-0.775 | 0.743 | 0.963 |
| Smooth | positive | 0.713-0.747 | 0.718 | 0.931 |
| Isotropic ambiguous | negative | 0.662-0.732 | 0.683 | 0.886 |
| Isotropic ambiguous | positive | 0.701-0.794 | 0.773 | 1.002 |

All fail the unchanged 1.25 threshold. Moreover, paired Chamfer decreases from
roughly 2.3-2.7 voxels at 1% growth to 0.68-0.77 voxels at 42.38%. As point
clouds become denser, nearest-neighbour distances shrink even when their
underlying distributions differ. This density dependence makes symmetric
Chamfer unsuitable as the sole measure of facet evolution across growth depths.
Equivalent volume/area displacement, signed crystallographic contrasts,
support-function advance, and facet-area changes must accompany it.

## 6. Answers to the scientific questions

1. **Is the previous weak response primarily explained by sub-voxel advance?**
   Sub-voxel depth is a major and demonstrable contributor: 10% growth is only
   0.258 average voxel layers, and deeper growth produces clearer signed
   equivalent-displacement contrasts. It cannot be identified as the sole or
   uniquely primary cause because the metric and kinetic mapping are also
   structurally biased.
2. **Does more deposited volume reveal systematic facet response?** Yes at the
   event and area-normalized equivalent-displacement levels, especially for
   negative alpha. No at the unchanged whole spatial-signal gate.
3. **Is it robust to stochastic variability and edge policy?** No. Negative
   alpha is consistent across the three seeds; positive alpha is only two of
   three, and effect magnitude depends materially on edge policy.
4. **Does site weighting map physically to local facet velocity?** No. The
   mapping contains candidate-density and digital-staircase factors and changes
   under the edge fallback policy.
5. **Is a modified law required?** Yes. The next legitimate step is to specify
   and independently verify a mass-conserving facet/interface-patch flux law,
   not to add empirical parameters to the present selector.

## Reproducibility

Run the exploratory gate before the controlled configuration:

```powershell
python phase2_ca/run_phase2a3b2_diagnosis.py --stage exploratory --config phase2_ca/configs/phase2a3b2_exploratory.json --output phase2_ca/outputs/phase2a3b2_exploratory_reproduced
python phase2_ca/run_phase2a3b2_diagnosis.py --stage controlled --config phase2_ca/configs/phase2a3b2_controlled.json --output phase2_ca/outputs/phase2a3b2_controlled_reproduced
```

Raw analytical, planar, paired-event, stage, facet, mask, runtime, gate,
manifest, plot, and SHA-256 files are retained under the corresponding output
directories. The manifests preserve the prior NO-GO and set
`physical_validation_claimed` to false.

## Final recommendation

Retain the conservative cache and numerical diagnostics, but stop the current
site-weighted CA as a physical morphology law. Conditional GO is limited to
formulating and testing the area/flux closure above on analytical planar and
single-crystal benchmarks. Do not proceed to multigrain contact, permeability,
vapour diffusion, liquid-film flow, or PNM coupling.
