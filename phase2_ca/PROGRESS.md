# Progress

## 2026-11-06 — Phase 2A minimal CA

- Added separate 3D material, liquid-accessibility, and seed arrays.
- Added a prescribed bottom substrate and validated multiple surface seeds.
- Implemented six-face salt-network connectivity from wetted seed anchors.
- Allowed empty air targets adjacent to supply-connected salt; target liquid
  occupancy is not required.
- Added an explicit control that bypasses the wetting filter.
- Coupled geometry realization transactionally to the verified Phase 1.1 mass
  allocator and kept initial seed mass separate.
- Added morphology and mass-closure diagnostics plus a deterministic
  evaporation-driven demonstration.
- Passed 18 Phase 2A tests and all 87 combined repository tests.

Status: numerical prototype only. No film-flow physics or physical validation
is claimed.

## 2026-11-06 — Phase 2A.1 growth-rule validation

- Preserved the legacy deterministic synchronous control and added weighted
  stochastic synchronous/sequential modes with persistent NumPy RNG state.
- Added transaction abort/rollback and injected-failure tests.
- Ran 44 CPU cases across three weight sets, five rule seeds, three timesteps,
  two update modes, seed topology, wetting, and capacity controls.
- Verified exact sequential 30/60/120 s event-sequence invariance and strong
  synchronous timestep grouping sensitivity at matched salt inventory.
- Saved 3D masks, orthogonal sections, connectivity/void tables, and opaque RGB
  figures with a SHA-256 inventory.
- Passed 30 Phase 2 tests and all 99 combined repository tests.

Status: GO as a stochastic numerical test bed; NO-GO for Phase 2B transport
coupling because morphology and inferred porosity remain rule/envelope sensitive.

## 2026-11-06 — Phase 2A.2 morphology feasibility audit

- Fixed the sequential rollback hole by restoring complete geometry, allocator,
  and RNG snapshots after any failed tentative update.
- Added failure injection after the first and fourth sequential insertions and
  during candidate regeneration.
- Added computational-box, bounding-box, column, and face-shell void masks with
  explicit external-air, enclosed-void, and bottom/top/lateral connectivity.
- Ran 18 paired 5/10/20 µm stochastic sequential cases at fixed physical
  domain, seed, duration, evaporation area, and supplied salt mass.
- Reanalyzed ten hash-verified Phase 2A.1 open/compact masks.
- Found no joint grid convergence and no uniquely identifiable internal pore
  network; apparent spanning paths are dominated by exterior air.
- Saved raw tables, masks, categorical 3D/cross-section figures, provenance,
  and a SHA-256 inventory under outputs/phase2a2/.
- Passed 38 Phase 2 tests and all 107 combined repository tests.

Status: numerical transaction and diagnostics accepted; NO-GO for Phase 2B.
Next gate is a grain/orientation and film-access morphology pilot validated
against independent 3D evidence.

## 2026-11-06 — Phase 2A.2R corrected pore diagnostics

- Proved and tested that the column envelope's implicit upper plane contains no
  empty cells by construction.
- Replaced envelope-derived pathway boundaries with explicit physical inlet,
  top-outlet, and lateral-outlet masks.
- Added five analytical geometry benchmarks separating exterior bypass,
  through-channel, intergranular gap, closed cavity, and blind pore.
- Replaced exact-envelope-equality porosity logic with explicit fixed-window,
  envelope, intergranular, enclosed, and ambiguous-exterior quantities.
- Added fixed 20 µm shell comparisons across 5/10/20 µm grids.
- Hash-verified and reanalyzed six selected Phase 2A.2 masks.
- Designed, but did not implement, a minimal grain-ID/orientation/film-access
  growth architecture.
- Passed 44 Phase 2 tests and all 113 combined repository tests.

Status: conditional GO for a bounded grain-aware morphology pilot; NO-GO for
physical calibration claims, permeability, liquid transport, or Phase 2B.

## 2026-11-06 — Phase 2A.3a normals and crystallographic orientation

- Added signed-distance-gradient and physical-radius weighted-plane-fit normal
  estimators with an outward signed-distance convention.
- Separated development and locked sets; selected the 120 µm plane-fit radius
  on four development geometries before running the independent locked set.
- Verified a 33° plane, shifted sphere, and two non-grid-aligned rotated cubes
  at 40/20/10 µm without changing the locked parameters or thresholds.
- Passed analytical `<100>/<110>/<111>`, alpha-zero, sign, 24-fold cubic
  symmetry, and joint-rotation invariants.
- Demonstrated refinement of inverse-rotation normal and anisotropy fields;
  finest-grid joint-weight RMSE was 0.000749.
- Retained the misleading aligned one-voxel volume result and added a labelled
  post-lock translation/rotation audit; one- and two-voxel seeds remain NO-GO.
- Set bounded empirical minima of 4 voxels for robust volume, 8 for normals, 16
  for facet orientation, and 5 across a gap for geometric width resolution.

Status: conditional GO for a resolved single-grain CA orientation pilot only.
No growth physics, multi-grain interaction, film flow, permeability, vapour
diffusion, PNM coupling, or physical calibration was implemented.
