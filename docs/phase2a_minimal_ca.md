# Phase 2A minimal 3D evaporation-driven salt-growth CA

## Scope and evidence boundary

Phase 2A implements the first independent CPU 3D geometry prototype using the
verified Phase 1 evaporation supply and Phase 1.1 transactional mass ledger.
It demonstrates a reproducible spatial update, strict numerical mass closure,
and supply-connectivity diagnostics.

It does **not** establish that a continuous liquid film exists, predict its
flow rate, or validate NaCl morphology. No PNM flow, vapor diffusion,
permeability, capillary disconnection, dissolution, or morphology-dependent
evaporation model is present.

## Physical state

The state uses three separate same-shaped arrays:

| Array | Meaning |
|---|---|
| `material: uint8` | `0` empty, `1` original substrate, `2` precipitated NaCl |
| `liquid_accessible: bool` | Direct liquid-film supply anchor; it is not a material occupancy flag |
| `seed_mask: bool` | Initial NaCl seeds, retained separately from new precipitation |

The bottom `substrate_layers` planes form an immutable prescribed substrate.
Seeds must be unique and lie immediately above the substrate interface. Initial
seed mass is `N_seed rho_s dx^3`; it is never charged to the cumulative
precipitation budget.

An empty target may represent air or liquid. Target liquid occupancy is not an
eligibility requirement. The model deliberately separates the location where
solid may appear from the conceptual path by which dissolved salt reaches the
adjacent crystal surface.

## Continuous-wetting prototype rule

A seed is a film-supply anchor when its salt voxel is marked
`liquid_accessible`. The algorithm recomputes the six-face-connected component
of salt attached to any such anchor. An empty voxel is eligible when it shares
a face with this supply-connected salt network.

With `enforce_wetting=False`, every existing salt voxel can parent growth. This
is an explicitly nonphysical control for isolating the consequence of the
wetting filter; it is not an alternative film model.

Connectivity is non-periodic. Edge voxels do not wrap to the opposite face.

## Geometry and transactional mass update

At each step:

1. The Phase 1 evaporation model supplies newly available NaCl mass.
2. `TransactionalMassToVoxelAllocator.prepare` converts available plus deferred
   mass to a requested whole-voxel count without committing solid.
3. The CA recomputes supply-connected salt and eligible empty neighbors.
4. Candidates are ranked by descending local six-face salt coordination. Ties
   use a pinned mixed-coordinate priority that avoids a systematic low-`z`
   ordering bias and is reproducible across runs.
5. At most `min(requested, eligible)` candidates are changed to salt.
6. Only the realized count is committed. All other mass remains explicitly
   deferred or in solution.

The tie-break and coordination rules are numerical prototype choices. They are
not derived from NaCl crystallography, interfacial energy, or growth kinetics.

The enforced identities are

`N_salt = N_seed + N_realized,cumulative`

and

`M_new solid + M_deferred + M_solution = M_supplied,cumulative`.

Any external salt voxel not represented in the mass ledger is detected as an
invariant violation.

## Deterministic demonstration

The checked-in demonstration uses:

- domain `16 x 16 x 12`, `10 um` voxels, one substrate layer;
- four wetted surface seeds;
- `T=298.15 K`, `RH=0.30`, `a_w=0.753`, `k_v=0.005 m/s`;
- a `16 x 16` voxel evaporation footprint;
- 12 steps of 60 s;
- Phase 1 NaCl supply `2.8735288999e-11 kg` per step;
- fixed tie-break seed `20261106`.

Final numerical result:

| Diagnostic | Result |
|---|---:|
| Initial seed voxels | 4 |
| Newly committed voxels | 159 |
| Total NaCl voxels | 163 |
| Supply-connected NaCl voxels | 163 |
| Occupied salt layers | 3 |
| Maximum height above substrate | 3 voxels (`30 um`) |
| Footprint | 62 voxels |
| Initial seed mass | `8.64e-12 kg` |
| Newly committed salt mass | `3.4344e-10 kg` |
| Retained sub-voxel/deferred mass | `1.3834679823e-12 kg` |
| Mass closure error | `0 kg` in the recorded floating-point evaluation |

The result proves only that the specified algorithm grows a 3D, connected,
mass-conservative lattice under its imposed rules. It is not evidence that the
shape, height, or growth rate is experimentally representative.

## Verification

The Phase 2A tests cover:

- independent material, liquid-accessibility, and seed arrays;
- bottom substrate and seed validation;
- wetted versus disconnected salt components;
- air-target eligibility without target liquid;
- non-periodic six-face connectivity;
- wetting-filter bypass control;
- transactional realization and seed/new-mass separation;
- insufficient geometric capacity and explicit retained mass;
- zero growth when no seed is wetted;
- connectivity recomputation after growth;
- deterministic candidate selection;
- detection of salt not represented in the mass ledger;
- byte-reproducible demonstration output and nonzero 3D height.

The demonstration writes a SHA-256 manifest beside its JSON output.

```text
python -m pytest phase2_ca/tests -q
18 passed

python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests phase2_ca/tests -q
87 passed
```

## Remaining limitations and next gate

- The substrate is a prescribed voxel block, not an imported PNM geometry.
- Liquid accessibility is an imposed seed-anchor mask; film thickness, flux,
  pressure, ion transport, and rupture are absent.
- Coordination ranking and tie-breaking can materially affect morphology and
  require sensitivity analysis before interpretation.
- Evaporation is spatially uniform and unaffected by growing salt.
- The demonstration is only 12 minutes and is not compared with experiment or
  the external CUDA CA.

Phase 2A is a **GO as a numerical interface prototype** and a **NO-GO for
physical morphology claims**. A later task should first test sensitivity to
candidate-selection rules and wetting-anchor topology before adding transport
physics.
