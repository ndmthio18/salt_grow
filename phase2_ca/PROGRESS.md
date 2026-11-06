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

## 2026-11-06 — Phase 2A.3b-0 orientation-to-growth interface

- Mapped every eligible empty growth voxel to reconstructed normals at all
  adjacent salt-interface faces; no Cartesian face direction was substituted
  for the reconstructed normal.
- Added smooth-average and isotropic-ambiguous edge/corner policies with
  candidate-level quality, class, eligibility, and weight records.
- Tested dimensionless fitting-radius transfer on 320/480 µm continuous cubes,
  16/24/32 cells across, two rotations, and two radius ratios.
- Preserved the original Phase 2A.3a locked data and added corrected resolution
  metadata; four and eight cells are excluded from facet evidence.
- Integrated cubic weights into the transactional sequential CA, including
  rollback of normal caches and diagnostics.
- Ran 42 paired CPU growth cases across two rotations, three alpha values, two
  edge policies, and three numerical RNG seeds; all mass closures were exact.
- Selected conservative local invalidation after exact agreement with full
  recomputation; rejected the faster stale adjacency-only cache.
- Passed 128 combined repository tests.

Status: conditional GO for a bounded full Phase 2A.3b single-grain numerical
pilot with both edge policies retained. NO-GO for a uniquely selected physical
edge rule, calibrated NaCl kinetics, multi-grain growth, or transport coupling.

## 2026-11-06 — Phase 2A.3b-1 anisotropic mechanism and morphology signal

- Separated method development from a SHA-256-locked verification protocol.
- Verified eight 100,000-draw fixed-front distributions on a sphere and rotated
  cube; maximum absolute and standardized family-probability errors were 0.002227
  and 2.073.
- Quantified face/edge/corner contributions at 16/24/32 cells across the grain;
  ambiguous-site influence did not diminish over the tested refinement range.
- Ran 36 paired single-grain trajectories to approximately 1%, 5%, and 10%
  added volume using two rotations, two edge policies, three alpha values, and
  three random seeds.
- Verified exact 256-event full/local cache sequences and exact periodic full
  checks across all 108 growth stages; median cache speedup was 7.76-fold.
- Confirmed event-family response to alpha, exact mass closure within 1.03e-25
  kg, and bounded joint-rotation error.
- Added a separate substrate-contact control with one prescribed supply anchor
  and Phase 1 RH/duration mass budget; no film hydraulics are claimed.
- The locked morphology signal/baseline ratios were only 0.893-0.958 versus the
  predeclared 1.25 minimum.

Status: numerical selector and cache accepted; **NO-GO for a complete
single-grain anisotropic morphology model**. Do not proceed automatically to
multigrain growth or transport coupling.

## 2026-11-06 — Phase 2A.3b-2 weak-anisotropy mechanism diagnosis

- Preserved the Phase 2A.3b-1 failed 1.25 morphology-signal gate unchanged.
- Calculated that 10% uniform cubic growth advances a 160 micrometre facet only
  2.582 micrometres or 0.258 of a 10 micrometre voxel; 42.3828% is required for
  one mean voxel layer.
- Reanalyzed 14,760 existing event records and separated alpha, seed,
  rotation/grid, and edge-policy distributions without rerunning them.
- Proved with frozen `<100>/<110>/<111>` planes that site-weighted displacement
  equals weight multiplied by orientation-dependent candidate density, not
  weight alone.
- Found that the isotropic ambiguous policy suppresses crystallographic weights
  on digital `<110>/<111>` staircase planes.
- Passed a 1735-event feasibility gate, then ran a frozen 15-case, three-seed,
  five-depth controlled matrix in a boundary-clear 48 cubed domain.
- Verified exact periodic full/local cache diagnostics and mass closure within
  1.03e-25 kg across 75 controlled stage states.
- Found signed area-normalized facet response, strongest for negative alpha,
  but positive-alpha response remained seed-sensitive and policy-dependent.
- All final Chamfer signal/baseline ratios remained 0.886-1.002, below the
  unchanged 1.25 threshold; Chamfer also decreased with point-cloud density.

Status: **NO-GO for the current site-weighted law and porous-salt continuation.**
A separately specified area/flux-based kinetic closure is required before any
further morphology physics. Multigrain and transport work remain out of scope.
