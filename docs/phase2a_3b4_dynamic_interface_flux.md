# Phase 2A.3b-4 — dynamic conservative interface-flux CA

## Decision

**NO-GO for the end-to-end dynamic planar-front closure.**

The new event-driven scheme passes global/local conservation, persistent mass
ownership, timestep partition, interface-cache consistency, candidate-collision,
domain-boundary, and complete rollback checks. It also removes the mathematical
candidate-bundle approximation by integrating `A*kappa` over distinct interface
faces before assigning volume to a candidate reservoir.

However, the locked growing-plane tests fail four predeclared front-accuracy
gates at 5 micrometres:

- alpha-zero actual orientation spread: 20.97% (limit 10%);
- alpha-zero partial-volume orientation spread: 12.16% (limit 5%);
- maximum actual reconstructed-front error: 25.43% (limit 15%);
- maximum partial-volume equivalent error: 13.14% (limit 5%).

The actual error improves with refinement, but improvement alone is not
acceptance. Per the registered stopping rule, the single-grain comparison was
not run, and no further model complexity was added.

The reviewed base is
`920aed2b8933801eb44ac65efada513e8cbc243f`. Work was performed on
`codex/phase2a3b4-dynamic-flux` in the main project directory. The previous
Phase 2A.3b-1 and 3b-2 NO-GO findings remain unchanged.

## 1. Scientific scope

This phase tests a geometric/kinetic numerical closure only.

- `kappa(n)` is a prescribed positive relative numerical factor.
- The total salt mass is externally supplied.
- Local supersaturation, liquid-film flow, vapour transport, dissolution,
  permeability, and pore-network coupling are not solved.
- Numerical RNG seeds and grid refinements are computational checks, not
  experimental replicates.
- Numerical closure does not constitute physical NaCl validation.

## 2. Why the Phase 2A.3b-3 adapter is not the required closure

For all reconstructed face elements owned by candidate `i`, the required
integrated weight is

\[
W_i^{\mathrm{faces}} = \sum_{j\in i} A_j\kappa(\mathbf n_j).
\]

The earlier evolving adapter instead computes

\[
W_i^{\mathrm{adapter}}
= A_i^{\mathrm{approx}}\kappa(\overline{\mathbf n}_i),
\]

where the normal is averaged before evaluating the nonlinear cubic function.
In general,

\[
\sum_j A_j\kappa(\mathbf n_j)
\ne
\left(\sum_j A_j\right)
\kappa\left(\frac{\sum_j A_j\mathbf n_j}
{\left|\sum_j A_j\mathbf n_j\right|}\right).
\]

### 2.1 Analytical perpendicular-face counterexample

For two equal-area perpendicular faces with normals `+x` and `+y`, the facewise
weight per unit face area is

\[
W^{\mathrm{faces}}/A = 2\exp(2\alpha/3).
\]

The bundled normal is `<110>`, so

\[
W^{\mathrm{adapter}}/A = 2\exp(\alpha/6),
\qquad
\frac{W^{\mathrm{adapter}}}{W^{\mathrm{faces}}}
= \exp(-\alpha/2).
\]

At `alpha = +(3/2) ln 2`, the adapter/facewise ratio is 0.5946, a 40.54%
underestimate. At negative alpha it is 1.6818, a 68.18% overestimate. The
implemented analytical check agrees within `2.22e-16` relative error.

### 2.2 Reconstructed-geometry audit

The 10 micrometre audit used identical plane, sphere, convex edge, concave
corner, and rotated-cube masks. Every internal salt/empty face was counted once,
with the Phase 2A.3b-3 L1 area correction and a 40 micrometre physical
plane-fit radius.

At nonzero alpha:

| geometry | maximum candidate discrepancy | total-weight discrepancy |
| --- | ---: | ---: |
| plane | floating noise only | 0% |
| convex edge | floating noise only | 0% |
| sphere | 1.48% | about 0.4% |
| concave corner | 4.46% | 0.35–0.77% |
| rotated cube | 3.49% | 0.23–0.27% |

The planar and convex controls agree because their candidate-owned faces have
equivalent reconstructed normals. This does not rescue the approximation; the
analytical counterexample and curved/concave results prove that normal averaging
and nonlinear kinetic evaluation do not commute.

## 3. Persistent local mass ownership

The implemented owner is `(front ID, candidate voxel coordinate)`, not a
transient face or triangle number. Distinct interface faces first contribute to
one candidate-integrated weight,

\[
w_i = \sum_{j\in i}A_j\kappa(\mathbf n_j),
\qquad
p_i = \frac{w_i}{\sum_k w_k}.
\]

Let `r_i` be deposited fractional solid volume stored at candidate `i`, and
`V_vox = Delta x^3`. The remaining capacity is

\[
g_i = V_{\mathrm{vox}} - r_i.
\]

For the current interface, the total additional volume that would first fill
any candidate is

\[
\Delta V_* = \min_i \frac{g_i}{p_i}.
\]

Each transaction allocates the lesser of the remaining input and `Delta V_*`.
All candidates reaching `V_vox` in that event batch become solid exactly once.
The interface is then updated and the weights are recomputed before more input
is consumed. This event substepping prevents a large input transaction from
overfilling a candidate and silently transferring the excess elsewhere.

The global identity after every update is

\[
V_{\mathrm{in}}
= N_{\mathrm{realized}}V_{\mathrm{vox}}
+ \sum_i r_i
+ V_{\mathrm{deferred}}.
\]

Fractional volume in `r_i` is physically deposited but geometrically unresolved
solid. A reservoir whose owner is temporarily inactive is reported separately
as inactive deposited fractional volume and is not redistributed. In contrast,
`V_deferred` is incoming precipitation that could not be assigned because no
admissible interface element existed. These terms are never double counted.

## 4. Interface evolution and local reconstruction

Every accepted event cycle performs:

1. enumerate unique internal salt/empty faces;
2. reconstruct a physical normal at every face;
3. correct its digital area by the local L1 projection factor;
4. evaluate `kappa` on each face normal;
5. integrate face weights by persistent candidate owner;
6. allocate the event-limited input to local reservoirs;
7. create all newly full voxels once;
8. update interface topology and locally affected normals;
9. retain all surviving and inactive owner inventories;
10. verify local bounds and global volume closure.

The local interface cache keys each face by axis and lower-cell coordinate.
Inserting a voxel changes only its six topological faces. All face normals
within `fitting radius + sqrt(3) Delta x` are conservatively invalidated and
recomputed with the same Gaussian-weighted PCA rule. Locked local caches had
exact face keys, validity flags, areas, and normals relative to a fresh complete
reconstruction within the registered tolerances.

Patch creation and disappearance were tested explicitly: dynamic insertion
created and removed face keys and activated new candidate owners, while the
cumulative local-allocation record retained all original owner identities.
Concave multi-face candidates integrate all owned faces before one reservoir is
updated. No exterior-domain candidate is synthesized.

## 5. Transaction safety

Before each update the implementation snapshots:

- all material masks and front versions;
- active and inactive local reservoirs;
- cumulative local allocations;
- global input, allocated, deferred, and realized inventories;
- per-front allocations and voxel counts;
- interface caches and reconstruction counters;
- transaction history and RNG state.

An injected failure after the first realization batch restored geometry,
reservoirs, every ledger, cache keys, and RNG state exactly.

## 6. Method development and locked protocol

### 6.1 Preserved development v1 failure

The original development configuration required at least 10% discrepancy in
the selected concave geometry. It observed 4.46%, so its stored decision is
`STOP_NO_GO`, while conservation, timestep, rollback, and cache gates passed.

This was an invalid success condition: the task requires detecting and
quantifying non-equivalence, not forcing one arbitrary geometry to exceed a 10%
effect. The v1 files and failed output were preserved unchanged.

### 6.2 Development v2 correction

Before locked outcomes, v2 changed only the discrepancy-detection threshold to
`1e-6`, safely above floating noise. The analytical counterexample, all
conservation thresholds, partition checks, rollback requirements, cache
tolerances, geometries, and dynamic smoke case were unchanged. V2 passed and
authorized the locked run.

The v2 configuration SHA-256 is
`25bade212b54761e70d32658de4185655cd52f46166bd5e66c3f9047189e3049`.

### 6.3 Locked protocol

The locked configuration SHA-256 is
`29ac6a5a2af919f6fb60d56a5bd327cd2d01723f8ddc8f4fbf4689dca45241c0`.
It fixed:

- `<100>`, `<110>`, and `<111>` circular front windows of 40 micrometre radius;
- a 200 micrometre cubic domain and 40 micrometre fitting radius;
- 20, 10, and 5 micrometre grids;
- one identical total input volume corresponding to 20 micrometres of
  alpha-zero advance on three exact circular patches;
- `alpha = 0` and `alpha = ±(3/2) ln 2`;
- actual area-weighted reconstructed-front position, partial-volume equivalent
  advance, realized-volume equivalent advance, surface area, inactive fraction,
  conservation, runtime, and memory metrics;
- 2, 5, and 10 equal input partitions at 10 micrometres;
- 10% finest alpha-zero actual orientation-spread limit;
- 5% finest alpha-zero partial-volume spread limit;
- 15% finest actual advance-error limit;
- 5% finest partial-volume advance-error limit;
- improvement relative to the 20 micrometre actual error;
- unchanged mass, cache, rollback, and ownership gates.

The actual calculation always used reconstructed areas and normals. Analytical
areas were used only to define the common physical input and reference advance.

## 7. Locked dynamic planar results

### 7.1 Alpha-zero refinement

| voxel size | actual orientation spread | partial-volume spread | maximum actual error | maximum partial error |
| ---: | ---: | ---: | ---: | ---: |
| 20 micrometres | 167.59% | 8.14% | 100.00% | 9.54% |
| 10 micrometres | 36.94% | 4.23% | 53.14% | 5.33% |
| 5 micrometres | 20.97% | 12.16% | 25.43% | 7.23% |

Actual front error improves strongly, so there is a refinement signal. The
finest grid nevertheless fails both actual gates. Partial-volume orientation
spread is non-monotone and worsens from 10 to 5 micrometres, indicating that
evolving reconstructed area and normal changes still introduce orientation
bias even before integer-front position is evaluated.

At 5 micrometres, alpha-zero actual `<100>/<110>/<111>` advances were 14.91,
17.10, and 18.44 micrometres for a common 20 micrometre reference. Their
partial-volume equivalent advances were 18.55, 20.23, and 20.98 micrometres.

### 7.2 Nonzero-alpha cases

At 5 micrometres:

| alpha | maximum actual relative error | maximum partial-volume error |
| ---: | ---: | ---: |
| `-(3/2) ln 2` | 22.65% | 0.19% |
| `0` | 25.43% | 7.23% |
| `+(3/2) ln 2` | 23.06% | 13.14% |

The negative-alpha partial allocation is close to its ideal ratios, but the
actual voxel front is not. Positive alpha also shows substantial error in the
partial-volume result, so integer quantization is not the only limitation.

### 7.3 Conservation and timestep partition

The maximum locked global volume residual was `1.51e-28 m3`, well below the
`1e-24 m3` limit. Inactive fractional volume never exceeded `2.53e-29 m3` and
was retained rather than redistributed.

At 10 micrometres, one cumulative input and 2, 5, or 10 equal increments gave:

- exactly identical final masks;
- identical realized voxel counts and event-batch counts;
- maximum local-reservoir difference below `1e-30 m3`;
- allocated-volume differences no larger than `5.05e-29 m3`;
- closed global ledgers.

Thus fractional ownership removes artificial timestep partition dependence for
the tested physics-free input sequence. This success does not compensate for
the failed spatial front-accuracy gates.

## 8. Computational cost

The host audit found 22 visible logical CPUs, no usable local GPU, unknown RAM
capacity, and about 66.9 GB free disk. Locked cases used at most three worker
processes and one native thread per process.

At 5 micrometres, each three-front case required 470–533 event batches and
279–325 seconds of worker time. The nine locked cases summed to about 930
seconds of per-case runtime. Python-tracked peak allocations were 9.5–10.6 MB;
this excludes NumPy native buffers and must not be interpreted as total process
RSS. Observed worker RSS during execution was about 86–87 MB.

The local patch update made the calculation feasible relative to full
reconstruction after every event, but the method remains expensive for a small
three-front verification. This is an implementation limitation, not evidence
about physical kinetics.

## 9. Reproducibility

From the repository root:

```powershell
$env:PYTHONPATH=(Resolve-Path '.\phase1_mass_budget\src').Path + ';' + `
  (Resolve-Path '.\phase2_ca\src').Path
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'

# Preserved failed development v1
python .\phase2_ca\run_phase2a3b4_validation.py `
  --config .\phase2_ca\configs\phase2a3b4_method_development.json `
  --output .\phase2_ca\outputs\phase2a3b4_development

# Corrected pre-lock development v2
python .\phase2_ca\run_phase2a3b4_validation.py `
  --config .\phase2_ca\configs\phase2a3b4_method_development_v2.json `
  --output .\phase2_ca\outputs\phase2a3b4_development_v2

# Locked verification
python .\phase2_ca\run_phase2a3b4_validation.py `
  --config .\phase2_ca\configs\phase2a3b4_locked_verification.json `
  --output .\phase2_ca\outputs\phase2a3b4_locked
```

All three output directories contain verified `SHA256SUMS.txt` inventories.
The combined Phase 1, Phase 2, and M4 suite passed 150 tests. Static unit and
uncertainty audits reported zero findings in the new model and runner.

## 10. Scientific interpretation and stop decision

The phase establishes that a dynamic, coordinate-owned fractional-volume ledger
can be conservative, transaction-safe, and invariant to input partitioning. It
also establishes that the previous candidate-area adapter is not generally
equivalent to facewise interface-flux integration.

It does **not** establish an acceptable dynamic crystallographic front law.
Even with exact mass closure, the reconstructed-area evolution and integer
geometry retain excessive orientation-dependent error. Because the end-to-end
planar gate failed, instruction/15 requires stopping here.

No single-grain comparison, multigrain growth, PNM, permeability, film
hydraulics, or vapour transport was run. A future phase would require a new,
explicit authorization and should first decide whether subvoxel/level-set front
geometry is necessary; the present result does not authorize that replacement.
