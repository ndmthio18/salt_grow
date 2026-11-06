# Phase 2A.2 morphology feasibility and numerical robustness audit

## Decision

**The sequential CA remains useful as a conservative numerical test bed, but
Phase 2B transport coupling is a NO-GO.** The transaction defect is fixed and
the tested mass ledger closes exactly. Nevertheless, the generated geometry is
not grid converged, an intrinsic pore volume is not uniquely identifiable, and
all apparent spanning voids are dominated by exterior-air-accessible space.

This result exposes an evidence–conclusion mismatch for a transport-ready
porous salt: mass conservation establishes inventory closure, not morphology
convergence; computational-box percolation establishes an open air route, not
an internal liquid pathway; and connected porosity does not establish
permeability.

## Scope and evidence boundary

This audit starts from commit 098f5bb8cec41ab28502fa3d6b1126057a14901a.
It does not implement a pore network, permeability calculation,
vapor-diffusion solver, or Phase 2B coupling. Simulation seeds quantify
algorithmic variability only and are not experimental replicates; tables
therefore use raw values, means, sample standard deviations, and ranges without
inferential p-values.

The figures and tables are general research outputs, not formatted against a
specific journal's current production requirements.

## Exception-safe sequential transactions

The prior grow rollback depended on the coordinate list returned by
_sequential_selection. If selection raised after inserting one or more voxels
but before returning, the caller's selected list remained empty and the
tentative geometry survived.

The corrected step snapshots:

- the full material geometry;
- every mutable allocator field, including pending proposal and transaction ID;
- the NumPy bit-generator state.

Any exception during proposal preparation, candidate generation, stochastic
selection, tentative insertion, or commit restores all three snapshots and
then rechecks the geometry/ledger invariants. This also protects against a
future commit implementation that mutates part of the ledger before raising.

Failure-injection tests raise immediately after the first and fourth sequential
insertions, during stochastic selection after prior insertions, during
candidate regeneration after prior insertions, and during synchronous or
sequential commit. In every case geometry, cumulative input, committed mass,
realized voxel count, pending state, transaction counter, and RNG state equal
their pre-step values.

## Void and connectivity definitions

Six-face connectivity is used throughout. Void means a computational voxel
whose material is EMPTY; it is not automatically an intrinsic salt pore.

| Analysis mask | Exact definition | Interpretation |
|---|---|---|
| Computational box | Entire domain above the substrate | Detects domain-scale air routes; mostly exterior space |
| Salt bounding box | Tight axis-aligned box containing all salt | Loose envelope that includes gaps between branches or grains |
| Column envelope | Each occupied x-y column from substrate surface to its local salt top | Surface-supported envelope; excludes unsupported lateral gaps |
| Face shell | Salt plus one six-face-neighbor voxel layer, clipped above substrate | Very local, explicitly grid-scaled envelope |

For each mask, the audit labels void components and records exterior-air
access, enclosed voids, bottom-connected volume, bottom-to-top spanning, and
bottom-to-lateral spanning. Bottom, top, and lateral are the minimum and
maximum coordinate planes of that analysis mask. Accessible pore volume is the
empty volume inside the named mask that belongs to the domain-boundary
external-air component.

This is deliberately conservative: if the tested masks differ, the analysis
flags intrinsic porosity as non-unique rather than selecting a preferred
number.

## Grid-sensitivity design

The grid campaign uses the stochastic sequential algorithm and two Phase 2A.1
weight sets:

- uniform: (1, 1, 1, 1);
- strong high coordination: (1, 4, 16, 64).

The physical domain is 200 x 200 x 180 µm, the substrate is 20 µm thick,
and the initial crystal is the same 20 µm cube at x=80–100 µm,
y=80–100 µm, z=20–40 µm on every grid. Thus the seed is represented by
4³, 2³, and 1³ voxels at 5, 10, and 20 µm, respectively; its physical size and
position do not drift with resolution.

All 18 cases use the same 600 s duration, 60 s step, effective evaporation area
3.0789180511e-8 m², and total supplied NaCl mass 3.4560000000e-10 kg
(1.6000000000e-13 m³). This realizes exactly 1280, 160, and 20 new voxels at
5, 10, and 20 µm. Pinned RNG seeds are 101, 202, 303.

## Numerical closure and cost

- maximum relative supplied-mass mismatch across grids: 1.50e-16;
- maximum mass-closure residual: 0 kg;
- every final salt body is one connected solid cluster;
- summed measured simulation-loop runtime in the final recorded run: 58.20 s;
- the slowest case took 10.21 s.

Mean runtime fell from about 8.7–9.2 s at 5 µm to 0.08–0.09 s at 20 µm.
These are local wall-clock measurements for this implementation, not a
hardware-independent scaling law.

## Grid sensitivity

Values are mean ± sample SD across three algorithmic seeds.

| Rule | Voxel [µm] | Height [µm] | Footprint [µm²] | Surface/volume [mm⁻¹] | Characteristic size 6V/A [µm] | Bounding-box void fraction | Column void fraction |
|---|---:|---:|---:|---:|---:|---:|---:|
| Uniform | 5 | 75.0 ± 0.0 | 6667 ± 521 | 363.6 ± 31.3 | 16.58 ± 1.38 | 0.806 ± 0.009 | 0.307 ± 0.052 |
| Uniform | 10 | 70.0 ± 10.0 | 6933 ± 58 | 250.6 ± 17.1 | 24.01 ± 1.58 | 0.809 ± 0.029 | 0.278 ± 0.063 |
| Uniform | 20 | 86.7 ± 11.5 | 4800 ± 400 | 154.0 ± 17.6 | 39.29 ± 4.26 | 0.733 ± 0.059 | 0.279 ± 0.075 |
| Strong high coordination | 5 | 55.0 ± 0.0 | 5025 ± 109 | 129.2 ± 6.1 | 46.52 ± 2.26 | 0.606 ± 0.014 | 0.087 ± 0.030 |
| Strong high coordination | 10 | 56.7 ± 11.5 | 5267 ± 643 | 126.0 ± 3.5 | 47.65 ± 1.34 | 0.606 ± 0.065 | 0.075 ± 0.012 |
| Strong high coordination | 20 | 66.7 ± 11.5 | 4133 ± 833 | 123.8 ± 4.8 | 48.51 ± 1.87 | 0.582 ± 0.061 | 0.147 ± 0.039 |

Relative to 5 µm, the uniform rule at 20 µm changes footprint by -28.0%,
surface/volume by -57.7%, and 6V/A by +137.0%. The strong rule is more stable
in surface/volume (-4.1%) and 6V/A (+4.3%) but still changes height by +21.2%,
footprint by -17.7%, and column void fraction by +68.8%.

The fixed 160 µm window has the same solid fraction (0.041015625) for every
case because it contains the full matched inventory. That equality is purely a
mass statement. In the 80 µm window, values vary with grid and seed
(0.245–0.280 mean for uniform and 0.271–0.294 for strong), demonstrating that
spatial placement is not fixed by inventory.

**Conclusion: neither rule has demonstrated morphology convergence.** The
strong rule has a relatively stable face-area metric over this limited range,
but the remaining observables do not jointly approach a grid-independent
value. Three grids and three simulation seeds are a sensitivity screen, not a
formal asymptotic convergence study.

## Reanalysis of Phase 2A.1 open and compact structures

Ten saved 10 µm masks were reused after verifying each file against the
Phase 2A.1 SHA-256 inventory. Every case has the same 211 salt voxels (four
seed voxels plus 207 grown voxels).

| Rule (n=5) | Bounding-box void fraction | Column void fraction | Face-shell void fraction | Column void externally accessible | Bounding-box bottom-to-top | Column bottom-to-lateral | Enclosed voids |
|---|---:|---:|---:|---:|---:|---:|---:|
| Uniform/open | 0.854 mean | 0.256 mean | 0.623 mean | 0.977 mean | 5/5 | 4/5 | 0–4 voxels |
| Strong/compact | 0.816 mean | 0.091 mean | 0.535 mean | 1.000 mean | 5/5 | 5/5 | 0 voxels |

All ten computational boxes and bounding boxes have a bottom-to-top void
component. That route is exterior air around or between salt branches. No
column envelope has a bottom-to-top path. Bottom-to-lateral column paths occur
in 9/10 cases, but 95.4–100% of column-mask void is connected to exterior air.
The compact rule even reports column paths in all five cases while having no
enclosed void at all; this is a direct demonstration that a binary percolation
flag can describe an exposed surface gap rather than an internal pore network.

The three salt-envelope definitions disagree on pathway signature for every
case and span mean void-fraction ranges of 0.599 (uniform) and 0.725 (strong).
Intrinsic porosity is therefore flagged as non-identifiable in 10/10 cases.

**No case demonstrates a continuous, internally preserved pore pathway
suitable for a liquid-supply model.** Some cases satisfy a necessary geometric
connectivity test under one envelope, but none establish an internal pore
boundary, hydraulic aperture, wetting state, capillary continuity, or
permeability.

## Physical interpretation of the CA rule

The six-face weights can be interpreted only as coarse-grained numerical
preferences for filling sites with different local occupancy. They are useful
for spanning compact-to-branched algorithmic morphologies at fixed inventory.
They cannot be identified with atomistic NaCl kink, step, or terrace kinetics:
a 5–20 µm voxel contains an enormous number of lattice sites, has no crystal
orientation, and the rule has no detachment, surface-energy, or solution
supersaturation dependence.

The minimum physically motivated extension should add a small number of
structured ingredients, in this order:

1. **Persistent grain identity and orientation.** Each seed/grain needs a
   crystallographic frame so faceted growth and grain impingement are distinct
   from isotropic occupancy coordination.
2. **Multiple interacting crystals with an explicit collision rule.** Grain
   fronts should meet without automatically filling every intercrystalline
   gap; this is the mechanism needed to preserve candidate voids.
3. **A defined surface-film support domain.** Replace global solid
   connectivity with an explicit film-access mask or thickness field tied to
   the substrate and exposed crystal surfaces.
4. **Local supply accessibility.** Allocate the already-conservative mass
   according to local film/supersaturation access rather than only solid
   coordination.

These ingredients should first be calibrated or rejected against independent
3D morphology evidence. Adding several free coordination weights, a PNM, or a
permeability closure before that comparison would increase parameter
non-identifiability.

## Figures and machine-readable evidence

- grid_sensitivity.png: all grid/seed observations with mean ± sample SD.
- grid_morphologies.png: 3D and orthogonal views for the paired uniform
  seed-101 grids.
- phase2a1_open_compact_pores.png: matched-inventory open and compact
  structures.
- envelope_sensitivity.png: raw void fractions and pathway classifications.

In categorical sections, salt is blue, initial seed is yellow, substrate is
gray, the explicitly defined column-mask candidate pore is vermillion, and
exterior air is pale blue. Exterior air is transparent in 3D to avoid
occlusion. Slice indices are recorded in figure_slice_indices.csv; raw masks,
configuration, tables, source-mask hashes, and a complete SHA-256 inventory
are retained under phase2_ca/outputs/phase2a2/.

The four PNGs are opaque RGB at approximately 220 dpi. Color is redundant with
panel position and category boundaries. Automated screening flags the yellow
seed and pale-cyan exterior against white below the 3:1 graphical-object
contrast heuristic, and flags several blue/gray/vermillion pairs for grayscale
review. The rendered sections retain hard category boundaries and direct
legends, but the screen is a review aid, not a claim of universal accessibility
or journal compliance.

## Tests

    python -m pytest phase2_ca/tests -q
    38 passed

    PYTHONPATH=prat_dynamic_growth/src;phase1_mass_budget/src;phase2_ca/src
    python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests phase2_ca/tests -q
    107 passed

Both results are from the final regenerated-output worktree.

## Recommendation

Stop at Phase 2A.2. The next model-development step should be a **small
grain-resolved morphology pilot**, not Phase 2B:

1. retain the conservative sequential transaction and fixed physical-grid
   comparison;
2. add grain labels/orientations and a single explicit intergrain-gap
   preservation rule;
3. add a binary, auditable surface-film access field without flow;
4. compare the resulting 3D envelope and gap statistics against independent
   experimental morphology;
5. reconsider transport only if pore identity and key physical dimensions are
   robust to grid, seed, and envelope choice.
