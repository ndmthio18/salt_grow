# Phase 2A.3b-3 — conservative interface-area growth closure

## Decision

**GO for the bounded numerical interface-area closure only.** The independent
continuum allocation, selected corrected-digital surface estimator, frozen-front
integer realization, mass ledger, rollback, and local patch-data consistency
passed their locked numerical gates.

**NO claim is made that an evolving anisotropic NaCl morphology has been
validated.** The minimal single-grain comparison did not establish robust
directional improvement, and the new method has not been integrated with a
production per-patch reservoir update for an evolving grain. Multigrain growth,
film hydraulics, PNM, permeability, and vapour transport remain out of scope.

The reviewed base is
`58c8db851d8c5de91c6499e7b255ef37256bb530`; work was performed on
`codex/phase2a3b3-interface-flux` in the main project folder. The failed
Phase 2A.3b-1 and 3b-2 evidence was not altered.

## 1. Evidence boundary

This phase tests a numerical closure, not a calibrated material law.

- `kappa(n)` is a positive, prescribed, dimensionless relative factor.
- Analytical plane tests verify the allocation algebra independently of any CA.
- Surface tests quantify deterministic discretization error against known
  plane, sphere, and rotated-cube areas and normals.
- Frozen-front tests quantify integer-voxel error at fixed physical geometry
  and deposited volume.
- The three RNG seeds in the single-grain demonstration are numerical
  sensitivity replicates, not experimental replicates and not a basis for
  material inference.
- Reported area, normal, and advance errors are numerical errors relative to
  exact constructed geometries. They are not measurement uncertainties.

## 2. Continuum allocation law

For prescribed incoming salt mass `Delta M_s` and solid density `rho_s`,

\[
\Delta V_s = \frac{\Delta M_s}{\rho_s}.
\]

The interface is partitioned into non-overlapping patches with physical area
`A_i`, outward unit normal `n_i`, and positive factor `kappa_i`. Define

\[
W = \sum_j A_j\kappa_j,
\qquad
\Delta V_i = \Delta V_s\frac{A_i\kappa_i}{W},
\qquad
\Delta h_i = \frac{\Delta V_i}{A_i}
            = \Delta V_s\frac{\kappa_i}{W}.
\]

Consequently,

\[
\frac{\Delta h_i}{\Delta h_k} = \frac{\kappa_i}{\kappa_k},
\]

independently of candidate-site count and patch area. The implementation uses
`math.fsum` and assigns only the final floating residual to the largest-weight
patch, so the returned floating-point patch volumes sum exactly to the supplied
floating-point volume. This residual correction is bookkeeping; it is not an
unreported physical redistribution.

For verification only, the existing cubic numerical factor was retained:

\[
q_4(\mathbf n_c)=n_{c,x}^4+n_{c,y}^4+n_{c,z}^4,
\qquad
\kappa=\exp\left[\alpha(q_4-1/3)\right].
\]

At `alpha = +(3/2) ln 2`, the exact `<100>/<110>/<111>` factors are
`2`, `2^(1/4)`, and `1`; the negative-alpha values are their reciprocals.
These are numerical ratios, not NaCl kinetic constants.

## 3. Surface-patch definitions

### 3.1 Reconstructed triangulated interface

The first method builds an approximate signed-distance field from the binary
mask and reconstructs its zero isosurface with a fixed six-tetrahedra split of
each active cell cube. Each triangle is non-overlapping within that partition,
has an oriented area vector, and is linked to the nearest face-connected empty
candidate voxel. Domain-boundary surfaces are not invented for the half-space
benchmark.

This implementation was retained as a failed development alternative. On the
10 micrometre development grid, its absolute area errors were 6.28% (plane),
11.78% (sphere), and 10.80% (rotated cube). It failed both the predeclared 10%
maximum and the requirement that all three endpoint errors improve with grid
refinement. Its P90 normal errors also remained 27.6–47.1 degrees. It was not
used for locked allocation.

### 3.2 Corrected digital interface area

The selected method enumerates each internal salt/empty voxel face exactly
once. The face is linked to its unique empty-side candidate. A physical normal
is reconstructed by the already verified fixed-physical-radius local plane fit.
For normal `n`, the physical area represented by a digital face is

\[
A_{\mathrm{face}}^{\mathrm{corr}}
= \frac{\Delta x^2}{|n_x|+|n_y|+|n_z|}.
\]

Faces linked to the same candidate may then be bundled. Bundling preserves the
sum of face areas exactly and uses the area-weighted resultant normal. A
candidate touching multiple faces is explicitly labelled as a numerical
multi-face bundle; it is not called a physical crystallographic corner.

The locked 10 micrometre absolute area errors were 1.058% (rotated plane),
0.546% (shifted sphere), and 0.815% (rotated cube). Corresponding P90 normal
errors were 0.828, 1.319, and 1.262 degrees. The maximum 1.058% error passed the
predeclared 8% locked limit.

### 3.3 Numerical staircases versus physical edges

Two rules were defined before the locked morphology was inspected:

1. evaluate `kappa` on every unique corrected face and then sum fluxes;
2. bundle all faces owned by one candidate, construct one area-weighted normal,
   and then evaluate `kappa`.

On the locked rotated cube, the total `A*kappa` sensitivity between these rules
fell from 0.546% at 40 micrometres to 0.0258% at 10 micrometres. This bounded
numerical sensitivity does not identify a physical edge law. The analytical
cube supplies known edge/corner locations for the comparison; a digital
multi-face count itself is never used as proof of a real physical edge. Because
the 60 micrometre edge-exclusion band covers much of the small reference cube,
the edge-law result must not be extrapolated to unresolved experimental grains.

## 4. Conservative discrete realization

The continuum allocation is completed before any voxel is selected. Each patch
has its own mass reservoir and an ordered candidate queue. A voxel can be
created only if that patch reservoir contains at least one voxel mass,

\[
m_{\mathrm{vox}}=\rho_s\Delta x^3.
\]

After every transaction the ledger enforces

\[
M_{\mathrm{in}} = M_{\mathrm{voxel}} + M_{\mathrm{subvoxel\ reservoirs}} + M_{\mathrm{retained\ by\ constraints}}.
\]

`M_allocated = M_in` is checked separately. Candidate identifiers are globally
unique in the ledger: if two patches list the same candidate, at most one may
consume it. Mass left after candidate exhaustion is retained explicitly and is
not transferred to another patch. Geometry, reservoirs, consumed-candidate
state, totals, and transaction ID are restored together after any exception.

The locked injected-failure case rolled back exactly, consumed three distinct
candidates for three realized voxels despite an overlapping candidate list, and
closed with zero reported mass residual.

The current integer representation necessarily retains fractional patch mass.
It therefore approximates normal advance with a staircase. A subvoxel volume
fraction or level-set front would represent this fraction geometrically, but
this phase did not replace the CA architecture.

## 5. Predeclared development and locked protocols

Method development used
`phase2_ca/configs/phase2a3b3_method_development.json` (SHA-256
`0f7ab34869474a345e259c760e8f5330e89a343610ff4fa46028d03ecaae4348`).
The corrected-digital method was selected only after applying the written rule:
smallest maximum finest-grid error among methods that improve from the coarsest
to finest grid for all three geometries and remain below 10%.

The independent locked configuration is
`phase2_ca/configs/phase2a3b3_locked_verification.json` (SHA-256
`3a28d335b76972f4c03f005779b218466c2b6a6db74ca76cbad206ab97d48ad0`).
Before locked outcomes it fixed:

- analytical ratio error at or below `2e-14`;
- absolute volume closure at or below `1e-27 m3`;
- absolute mass closure at or below `1e-24 kg`;
- finest selected-area error at or below 8%;
- finest frozen-front advance error at or below 0.5%;
- finest/coarsest maximum-error ratio below 0.25 and fitted refinement order at
  least one;
- no candidate-capacity-limited patch;
- exact candidate keys, patch areas, and normals after conservative local
  invalidation, within `1e-24 m2` and `1e-12` respectively;
- exact transactional rollback and no double candidate consumption.

The local-update check used a fresh full reconstruction as the validation
oracle and replaced only records inside the conservative invalidation radius.
It proves patch-data consistency (129 affected candidates, zero difference),
not a production local-reconstruction speedup.

## 6. Verification results

### 6.1 Analytical planes

All `alpha = 0` advances were equal for `<100>`, `<110>`, and `<111>` despite
unequal areas. Nonzero-alpha advance ratios matched the exact prescribed
`kappa` ratios. Jointly rotating the physical planes and crystal frame preserved
the results. Repeating the calculations with candidate-count metadata of 1,
29, and 997 changed no allocation. The maximum ratio error was
`5.28e-16`, and allocated volume closed exactly in the stored rows.

### 6.2 Frozen-front discrete refinement

Physically identical 320 micrometre patches and the same total deposited volume
were evaluated at 20, 10, and 5 micrometres. The maximum new-method advance
errors were:

| alpha | 20 micrometres | 10 micrometres | 5 micrometres | fitted order |
| ---: | ---: | ---: | ---: | ---: |
| `-(3/2) ln 2` | 11.980% | 0.978% | 0.192% | 2.98 |
| `0` | 6.250% | 0.563% | 0.0452% | 3.56 |
| `+(3/2) ln 2` | 3.195% | 0.445% | 0.0440% | 3.09 |

The approximately third-order trend is expected for a fixed physical patch:
the unresolved remainder is less than one voxel volume, so its equivalent
advance error scales approximately as `Delta x^3/A`. It is not evidence that
the evolving CA geometry is a third-order front method.

At 5 micrometres, no patch was capacity limited and all volume identities
closed within the locked bound. Nonzero subvoxel reservoirs were retained and
reported rather than redistributed.

The old site-weighted law remained candidate-density dependent. At `alpha = 0`
its 5 micrometre advances were 1.358, 0.945, and 0.784 micrometres rather than
the exact common 1 micrometre advance; relative errors reached 35.8%. Across
the nonzero-alpha cases its maximum finest-grid error reached 43.8%. The new
area allocation reduced the corresponding maximum to 0.192%. This verifies
removal of candidate-count bias in the **allocation and frozen-front
realization**, not in all subsequent evolving morphology.

### 6.3 Minimal evolving single grain

The gate passed before this demonstration was launched. Old and new CA adapters
then inserted exactly 410 voxels (`4.10e-13 m3`) for each of three matched RNG
seeds. All six mass ledgers closed exactly. The new adapter uses corrected
candidate-owned area times `kappa`; the old adapter uses one site weight per
candidate.

Mean positive crystal-frame support advances were:

| method | `<100>` | `<110>` | `<111>` | mean runtime/trajectory |
| --- | ---: | ---: | ---: | ---: |
| old site weighted | 18.63 micrometres | 6.91 micrometres | 5.98 micrometres | 52.0 s |
| new area weighted | 15.58 micrometres | 8.63 micrometres | 0.00 micrometres | 63.6 s |

The new adapter was about 22% slower. More importantly, `<111>` support did not
advance for any of its three seeds, while other directions were seed sensitive.
Because support is voxel-quantized and this is one small realization, these
values neither prove improvement nor prove a physical failure. They do show
that removal of the planar candidate-density bias is insufficient to validate
an evolving anisotropic morphology. No morphology-level GO is issued.

## 7. Computational cost and reproducibility

The resource audit found 22 visible logical CPUs, no local NVIDIA/AMD GPU, and
unknown host memory availability. Calculations were therefore CPU-only with
BLAS threads limited to one and no more than six planned workers. The runner is
currently mostly serial. Corrected-digital locked surface cases required about
0.01–0.70 s each. The six sequential single-grain trajectories dominated cost
at about 52–68 s each.

Reproduction commands from the repository root are:

```powershell
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
python .\phase2_ca\run_phase2a3b3_validation.py `
  --config .\phase2_ca\configs\phase2a3b3_method_development.json `
  --output .\phase2_ca\outputs\phase2a3b3_development
python .\phase2_ca\run_phase2a3b3_validation.py `
  --config .\phase2_ca\configs\phase2a3b3_locked_verification.json `
  --output .\phase2_ca\outputs\phase2a3b3_locked
```

Both output directories contain manifests and verified `SHA256SUMS.txt`
inventories. The combined Phase 1, Phase 2, and M4 test suite passed 145 tests
when the three `src` directories were placed on `PYTHONPATH`. Static unit and
uncertainty audits reported zero findings in the new model and runner.

## 8. Limitations and next decision gate

1. The selected digital correction relies on resolved local normals and was
   verified only on constructed single-interface geometries over 10–40
   micrometre grids.
2. The triangulated alternative did not pass and cannot be used as independent
   confirmation of surface area.
3. Physical edge/corner kinetics remain unidentified. The two plausible
   partition rules were retained as a sensitivity, not selected visually.
4. The evolving CA adapter samples candidate probabilities proportional to
   corrected area times `kappa`; it does not yet carry the independent patch
   reservoirs through every evolving geometry update.
5. Integer voxels leave quantized support advances and subvoxel mass. A
   subvoxel/level-set representation is a conceptual alternative, not an
   implemented replacement.
6. `kappa` has not been calibrated against NaCl experiments, and no experimental
   morphology, local growth rate, or transport observable was validated.

Accordingly, Phase 2A.3b-3 stops here. The result permits further work only on a
separately authorized, reservoir-coupled evolving single-grain numerical gate.
It does not authorize multigrain growth or transport coupling.
