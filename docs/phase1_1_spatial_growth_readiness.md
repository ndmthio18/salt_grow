# Phase 1.1 spatial-growth readiness

## Decision

**Conditional GO for an independent, CPU-based Phase 2 CA prototype; NO-GO
for treating the original 64-cubed, 24 h reference case as a feasible spatial
domain.** The accounting and transaction boundaries are now explicit and
tested. However, the reference salt inventory exceeds even a fully solid
64-cubed domain. Phase 2 must either stop before capacity, use a larger domain,
or change a physically justified input. It must never discard excess salt or
create solid outside the grid.

This is numerical readiness for a bounded prototype, not validation of salt
morphology, wetting, transport, or the illustrative evaporation coefficient.
No 3D growth algorithm is implemented in Phase 1.1.

## Reference-case capacity audit

The audit independently combines the Phase 1 salt rate with the prescribed
domain geometry. Inputs are a `64 x 64 x 64` grid, `10 um` cubic voxels,
`T=298.15 K`, `RH=0.30`, and 24 h duration. The Phase 1 reference has no
initial seed voxels; future seeds must be counted as pre-existing solid volume.

| Quantity | Result |
|---|---:|
| Domain side length | `6.4e-4 m` (`0.64 mm`) |
| Domain volume | `2.62144e-10 m3` (`0.262144 mm3`) |
| Initial seed solid volume | `0 m3` |
| Cumulative precipitated salt mass | `6.6206105853e-7 kg` |
| Cumulative precipitated salt volume | `3.0650974932e-10 m3` (`0.306510 mm3`) |
| Solid-volume fraction of the domain | `1.1692419` (`116.924%`) |
| Volume beyond a fully solid domain | `4.4365749318e-11 m3` (`16.924%` of domain volume) |
| Time to theoretical fully solid capacity | `73894.03 s` (`20.526 h`) |

The capacity time assumes the constant Phase 1 RH=0.30 precipitation rate and
zero seed volume. It is an inventory extrapolation, not a claim that a
space-filling crystal geometry can actually form.

For the same `0.64 mm x 0.64 mm` footprint and 24 h solid inventory, the
continuous domain heights required at representative bulk solid fractions are:

| Permitted solid fraction | Required height | Minimum 10-um voxel layers |
|---:|---:|---:|
| 0.3 | `2.49438 mm` | 250 |
| 0.5 | `1.49663 mm` | 150 |
| 0.7 | `1.06902 mm` | 107 |

The new `diagnose_domain_capacity` function reports domain, seed, precipitated,
total, remaining, and overflow volumes; an explicit status; and an optional
constant-rate capacity time. Its configurable maximum solid fraction reserves
porosity. It is deliberately only a scalar feasibility screen: a passing
result does not establish that eligible sites are connected, wet, or
geometrically accessible.

## Transactional mass allocation for Phase 2

The existing `MassToVoxelAllocator` remains available for the already verified
Phase 1 scalar ledger. Phase 2 should instead use the new two-stage
`TransactionalMassToVoxelAllocator`:

1. `prepare(available_precipitation_mass_kg, ...)` returns a non-mutating
   proposal containing available mass, requested whole voxels, requested solid
   mass, and sub-voxel mass.
2. The geometry layer determines eligible sites and attempts to create at most
   the requested number on a tentative geometry.
3. `commit(proposal, realized_solid_voxels, unallocated_state=...)` records only
   successfully created voxels as solid.
4. All uncommitted mass is retained explicitly either as deferred
   precipitation or as mass remaining in solution.

The accounting identity after every committed transaction is

`committed solid mass + deferred mass + solution mass = cumulative supplied mass`.

Deferred mass is reconsidered in the next transaction. Solution mass is not
reconsidered until the caller explicitly releases it. This makes the physical
interpretation a model decision and prevents insufficient spatial capacity
from being silently reported as precipitation. Explicit outcomes include
`zero_growth_insufficient_mass`, `no_spatial_capacity`,
`partial_spatial_capacity`, and `committed`.

The proposal is intentionally mass-based and does not accept an eligible-site
count. Capacity and geometry remain independent checks: the global domain may
have scalar room while the local growth front has zero eligible sites.

## Proposed Phase 2 physical state

A single state array is not recommended because material identity and fluid
phase are different axes of information. A single enumeration would create
ambiguous combinations and make wetting changes unnecessarily invasive.

Use the following minimal representation:

- `material: uint8[nx, ny, nz]`: `0` pore, `1` original PNM substrate, `2`
  precipitated NaCl crystal.
- `pore_phase: uint8[nx, ny, nz]`: liquid or air, meaningful only where
  `material == 0`.
- `eligible: bool[nx, ny, nz]`: a derived, recomputed mask, never an independent
  source of truth.
- Optional audit arrays: integer `seed_id`, crystallization step/time, and
  parent site. These are diagnostic metadata, not extra material states.

An eligible growth site must simultaneously be pore space, liquid-accessible,
face-adjacent to existing NaCl, inside the domain, and not already reserved by
the current tentative transaction. Original PNM solid is immutable and cannot
be converted to salt.

### Multiple seed initialization

Seed coordinates should be supplied explicitly or generated reproducibly from
a pinned random seed. Candidate seed cells are pore voxels that are
face-adjacent to the original PNM substrate and liquid-accessible. Reject
duplicates, out-of-domain coordinates, substrate cells, and air-only cells.
For automatically placed multiple seeds, select without replacement and use a
documented spacing rule (for example farthest-point selection) to avoid an
unstated clustering bias. Seed voxels are pre-existing morphology and their
volume enters the capacity diagnostic, but not the newly precipitated mass
ledger.

### Candidate selection and mass coupling

For each physical time increment:

1. Compute the newly available salt mass and create a non-mutating allocation
   proposal.
2. Recompute liquid accessibility and eligible sites.
3. Score eligible sites by local crystalline coordination, initially the
   number of face-adjacent NaCl voxels. Pin the tie-breaking rule and random
   seed. Coordination changes selection priority, not the global mass budget.
4. Tentatively select at most
   `min(requested_solid_voxels, eligible_site_count)` sites without replacement.
5. Apply the tentative material update atomically. Only after it succeeds,
   commit the actual count to the mass ledger. If the tentative update is
   rolled back, commit zero realized voxels so the entire mass remains explicit.
6. Retain every unallocated kilogram in the declared deferred or solution
   state and emit the explicit transaction status.

The initial prototype may assume continuous wetting: every pore voxel connected
to the liquid reservoir is labeled liquid-accessible at every step, with no
rate limitation. This must be an explicit configuration flag. A test must show
that turning liquid eligibility off for otherwise identical candidate sites
produces zero growth and retained mass. This assumption does not establish ion
supply, capillary continuity, or constant evaporation during real growth.

## Reproducibility changes

- Time-series points now derive from configured total duration and physical
  timestep, include the final time, and no longer depend on a hard-coded
  1441-point grid.
- Matplotlib SVG IDs use a fixed hash salt, SVG date metadata is removed, and
  creator metadata is fixed. PNG software metadata is also fixed.
- Two independent output directories are generated byte-for-byte identically
  in an automated test under the same environment.
- The default 24 h, 60 s reference still contains 1441 points and preserves the
  previous numerical RH sweep and cumulative-mass values.

The documented local environment is Python 3.12.4, Matplotlib 3.9.2, and
pytest 9.1.1. Reproducible hashes are asserted only within the same documented
software environment; rendering-library or font changes can legitimately
change image bytes.

## Verification executed

The Phase 1.1 suite covers capacity overflow, multiple heights, configurable
solid fractions, seed volume, excess requests relative to spatial capacity,
zero growth, variable-RH conservation, repeated sub-voxel allocations,
transaction consistency, configured time grids, and byte-reproducible outputs.

```text
python -m pytest phase1_mass_budget/tests -q
53 passed

python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests -q
69 passed
```

## NVIDIA resource check

Among the requested remote systems, `frodo` was reachable through its existing
verified SSH configuration and exposed an NVIDIA GeForce RTX 2060 with 6 GB
VRAM. At inspection it was lightly loaded. `sal` exposed `nvidia-smi` but the
current session lacked NVML permission. `cxhpc` was not used because its name
had no locally trusted SSH host key and a verified usable GPU had already been
found.

The `frodo` environment had a compatible NVIDIA driver but no CUDA toolkit and
no PyCUDA package required by the audited external reference CA. Installing or
changing that environment is outside Phase 1.1. Therefore no original-CA GPU
baseline is claimed. The present deliverables are CPU accounting and interface
verification only.

## Remaining limitations

- Capacity passing is necessary but not sufficient for geometric feasibility.
- No CA growth, PNM transport, permeability, liquid redistribution, vapor
  diffusion, morphology feedback, dissolution, or deliquescence is modeled.
- The continuous-wetting rule is only a proposed, testable prototype
  assumption.
- The illustrative transfer coefficient, fixed area, and continuous saturated
  brine supply remain uncalibrated and physically unvalidated.
- Coordination-based selection is a future numerical rule, not an established
  NaCl growth mechanism.
- The remote reference CUDA baseline remains blocked by the exact software
  environment, not by lack of detected NVIDIA hardware.

The appropriate next step is a small, deterministic CPU CA prototype with a
short-duration or enlarged-domain test case, transaction rollback tests, and
full mass closure. It should not begin automatically from this Phase 1.1 task.
