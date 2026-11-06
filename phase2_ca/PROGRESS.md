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
