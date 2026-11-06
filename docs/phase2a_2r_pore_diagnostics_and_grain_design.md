# Phase 2A.2R corrected pore diagnostics and grain-aware CA design

## Decision

**GO for a bounded grain-aware morphology pilot; NO-GO for permeability,
liquid transport, PNM coupling, or a physically calibrated porous-salt claim.**

The corrected diagnostics pass five independent geometry benchmarks and remove
the invalid column-envelope outlet. The Phase 2A.2 transport conclusion remains
conservative, but its evidential basis changes: the present structures do not
demonstrate a physically identified liquid-supply pathway because neither the
salt-body boundary nor the liquid-film domain is known. The previous
column-envelope bottom-to-top negatives were not valid evidence that such a
path was absent.

The proposed pilot is a falsifiable numerical-development step. It does not
authorize Phase 2B and does not assume that suitable calibration data already
exist.

## Scope and claim boundary

This audit starts from commit
3f6fc075e88fed9a175960554b8325c13ccb2045. It preserves the evaporation-derived
salt budget, sequential event logic, transaction rollback, RNG reproducibility,
and existing verification infrastructure.

It does not implement grain growth, permeability, vapor diffusion, liquid
transport, capillary flow, or PNM coupling. Geometric void connectivity,
exterior-air access, and a physically justified liquid pathway are reported as
different constructs.

## Critical correction: the column envelope has no usable void outlet

For each occupied horizontal coordinate (x,y), the column envelope contains
cells from the substrate surface through that column's highest salt cell.
Let h(x,y) be that highest salt index and H be the maximum h over all occupied
columns.

At z=H, an envelope cell exists only in columns with h(x,y)=H. By definition,
the cell at z=h(x,y) is salt. Therefore

    column envelope ∩ empty phase ∩ implicit global upper plane = empty set.

The former upper boundary was therefore unavailable to the void phase by
construction. A false bottom-to-top result under that boundary is a diagnostic
artifact, not evidence against a physical pore pathway.

The automated test
test_column_envelope_implicit_global_top_has_no_empty_voxels_by_construction
demonstrates this property for columns of unequal height.

### Corrected boundary treatment

Analysis masks and boundary surfaces are now separate inputs:

- the analysis mask states where a path is allowed;
- the inlet and outlet masks state which physical surfaces the path must touch;
- external-air connectivity is labeled independently in the full domain.

The corrected API accepts explicit inlet, top-outlet, and lateral-outlet
surfaces. Fixed rectangular windows use coordinates specified independently of
the salt or void realization. The column envelope remains useful for an
envelope-dependent volume fraction, but its implicit top is no longer treated
as a physical pore outlet.

This change prevents a morphology-dependent boundary from predetermining the
connectivity answer.

## Independent geometry benchmarks

All benchmarks use six-face connectivity and explicit fixed-window surfaces.
The whole-window path is evaluated separately from the path inside the
intended material mask.

| Geometry | Whole-window bottom-to-top | Material bottom-to-top | External material void [vox] | Enclosed void [vox] | Accessible face-adjacent intergranular void [vox] | Reference result |
|---|---:|---:|---:|---:|---:|---|
| Solid pillar in open air | yes | no | 0 | 0 | not applicable | exterior bypass only |
| Porous slab with through-channel | yes | yes | 7 | 0 | not applicable | true channel through material |
| Two grains with connected gap | yes | yes | 35 | 0 | 35 | connected intergranular pathway |
| Closed internal cavity | yes | no | 0 | 1 | not applicable | enclosed cavity |
| Externally connected blind pore | yes | no | 3 | 0 | not applicable | external access without through-path |

Every computed answer matches the analytical voxel count and connectivity.
Most importantly, the solid-pillar case proves that a bottom-to-top path in the
computational window can be an exterior bypass around a nonporous object.

The intergranular diagnostic requires persistent positive grain IDs on all
salt cells and identifies an empty cell only when it face-touches at least two
distinct grains. This local definition is exact for the one-cell benchmark gap;
it is not yet a general segmentation of wide intergranular regions.

## Corrected porosity quantities

The implementation now keeps the following quantities separate:

1. **Fixed-window solid volume fraction.** Salt volume divided by a
   predeclared physical observation window. It is reproducible across grids
   when the same coordinates are used, but is not intrinsic porosity.
2. **Envelope-dependent void fraction.** Void fraction inside bounding-box,
   column, or physical-shell masks. The range across masks quantifies boundary
   sensitivity.
3. **Accessible intergranular void.** Externally accessible empty cells that
   also satisfy an explicit grain-ID-based intergranular definition.
4. **Enclosed void volume.** Empty components not connected to the exterior
   domain boundary.
5. **Ambiguous exterior space.** Exterior-connected empty volume inside a
   provisional salt envelope when no physical material boundary is supplied.

Intrinsic porosity is now reported only when the caller supplies an explicit
material-boundary mask and names its source. Exact equality of several
algorithmic envelope masks is no longer an identifiability criterion. With no
physical boundary, intrinsic porosity is undefined rather than zero or
estimated.

For cross-grid shell comparisons the shell thickness is fixed at 20 µm. This
corresponds to four, two, and one face-dilation layers on the 5, 10, and 20 µm
grids. A test verifies that the same 20 µm crystal produces the same 60 µm
horizontal shell extent on all three grids.

## Reanalysis of selected Phase 2A.2 masks

Six source masks were verified against the Phase 2A.2 SHA-256 inventory:
uniform and strong-high-coordination cases at 5, 10, and 20 µm, all with RNG
seed 101.

The fixed observation window is x=60–140 µm, y=60–140 µm, and z=20–100 µm.
The physical inlet is the z=20 µm layer and the top outlet is the final cell
layer below z=100 µm. The physical shell thickness is 20 µm.

| Rule | Grid [µm] | Fixed-window solid fraction | Fixed-window bottom-to-top | Exterior fraction of window void | Bounding-box void fraction | Column void fraction | 20 µm-shell void fraction | 20 µm-shell bottom-to-top | Enclosed bounding-box void [vox] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Uniform | 5 | 0.289 | yes | 0.985 | 0.797 | 0.359 | 0.763 | yes | 44 |
| Uniform | 10 | 0.293 | yes | 0.994 | 0.802 | 0.266 | 0.730 | yes | 2 |
| Uniform | 20 | 0.297 | yes | 1.000 | 0.738 | 0.222 | 0.677 | yes | 0 |
| Strong high coordination | 5 | 0.309 | yes | 1.000 | 0.623 | 0.114 | 0.661 | no | 0 |
| Strong high coordination | 10 | 0.303 | yes | 1.000 | 0.627 | 0.082 | 0.650 | no | 0 |
| Strong high coordination | 20 | 0.281 | yes | 1.000 | 0.650 | 0.192 | 0.661 | no | 0 |

The fixed-window path exists in all six cases, but 98.5–100% of its void is
connected to exterior air. This is an exterior-space result, not a liquid-film
result. The 20 µm shell distinguishes the open and compact rules in this
selected realization, but the shell is still a geometric neighborhood around
salt and cannot identify fluid occupancy.

The range between tested envelope void fractions is 0.438–0.568. Enclosed
bounding-box void in the uniform case falls from 44 voxels at 5 µm to 2 at
10 µm and zero at 20 µm. This is strong grid dependence, not evidence of a
resolved pore population.

Every column envelope again has zero empty cells on its implicit upper plane,
as predicted analytically. Column void fractions remain reported only as
volume-sensitivity quantities.

No selected mask has persistent grain IDs, a measured material boundary, or a
physical liquid-film mask. Consequently:

- accessible intergranular void is unresolved;
- intrinsic porosity is undefined;
- a liquid-supply pathway is not evaluable;
- no permeability inference is permitted.

## Change to the Phase 2A.2 scientific conclusion

### What remains supported

- The mass ledger is conservative.
- Morphology remains grid and rule dependent.
- Envelope choice changes inferred void fraction substantially.
- Current outputs do not justify Phase 2B transport or permeability.

### What is withdrawn or narrowed

The former observation that no column envelope had a bottom-to-top path cannot
be used as evidence that an internal path was absent. The outlet was blocked by
the envelope definition itself.

The corrected statement is:

**The current data do not identify whether a physically meaningful liquid
pathway exists.** Fixed-window paths are dominated by exterior bypass; shell
paths are geometric surface neighborhoods; column paths are not applicable;
and the required grain and film states are absent.

This is an insufficiency-of-identification conclusion, not proof of
nonexistence.

## Minimal grain-aware CA architecture

### State

The smallest useful extension adds:

- material state: existing EMPTY, SUBSTRATE, and SALT array;
- grain_id: unsigned integer, zero outside salt and persistent positive ID on
  every salt cell;
- orientation: one normalized quaternion or rotation matrix per grain,
  canonicalized under cubic symmetry;
- liquid_film_accessible: a separate Boolean or bounded accessibility field,
  never inferred from solid connectivity;
- grain-contact graph: grain-pair contact area and first-contact time;
- existing mass allocator, RNG state, and transactional ledger unchanged.

Grain ID and orientation are state or initial conditions, not empirical growth
coefficients.

### Sequential event

1. Prepare the evaporation-derived salt-mass proposal.
2. Generate empty interface candidates adjacent to salt.
3. Apply the independently defined liquid-film-access mask.
4. Determine adjacent parent grain IDs.
5. Estimate a physical interface normal and transform it into each candidate
   grain's crystal frame.
6. Compute one facet-anisotropy preference.
7. Apply an explicit grain-impingement/contact rule.
8. Sample or rank the remaining candidates, tentatively update material,
   grain ID, and contact state, and commit only realized voxel mass.
9. On any failure, roll back geometry, grain IDs, contact state, ledger, and
   RNG together.

If a contact or accessibility rule rejects growth, unallocated salt remains in
the existing deferred or solution pool. Rejection must never create or destroy
mass.

### Minimal facet parameterization

Let n_c be the interface normal in the grain crystal frame and define

    q4 = n_x^4 + n_y^4 + n_z^4.

A one-parameter cubic preference can use a normalized positive weight such as

    w_facet = exp(alpha × (q4 - 1/3)).

The dimensionless alpha controls anisotropy; orientation rotates the same law
for each grain. This is a coarse-grained morphology parameter, not an
atomistic kink rate. A reliable normal estimator is required; using only the
six lattice directions would confound physical anisotropy with grid
orientation.

Because the total deposition mass is prescribed externally, the pilot does not
need a separate absolute kinetic rate. Only relative site competition is
modeled.

### Impingement and pore preservation

- A candidate touching one grain inherits that grain ID.
- A candidate touching multiple grains is a contact event, not an automatic
  merger.
- Deposition across a gap is allowed only under a predeclared physical contact
  or liquid-bridge criterion; otherwise the candidate is rejected and its mass
  remains unallocated.
- Contact records preserve which grains met, where, and when.

Intercrystalline pores can then arise from measured seed spacing, orientation-
dependent faceting, growth-front impingement, restricted film reach, and local
supply shadowing. They must not be imposed through a target porosity or by
randomly protecting a desired fraction of void.

## Parameters and calibration status

| Quantity | Role | Dimension | Source needed |
|---|---|---|---|
| Seed centroid and size per grain | Initial condition | m | image measurement |
| Grain orientation | Initial condition | dimensionless rotation | facet geometry, EBSD, diffraction, or 3D orientation measurement |
| Facet anisotropy alpha | Relative growth preference | dimensionless | time-resolved facet displacement or morphology calibration |
| Film-access length or measured access mask | Growth eligibility | m, or direct field | surface-film observation |
| Contact/bridge length criterion | Grain-gap closure | m | observed gap/contact evolution |
| Supplied salt per step | Conserved budget | kg | existing evaporation model |
| Voxel size and timestep | Numerical controls | m, s | convergence study, not calibration |
| RNG seed | Reproducibility control | none | not physical |

The facet anisotropy and film/contact criteria are the only new growth-rule
parameters. Orientation, seed geometry, and spacing should be measured inputs
whenever possible rather than fitted free parameters.

## Same physical crystal across 5, 10, and 20 µm grids

Define every initial grain as a continuous oriented polyhedron with a physical
centroid, dimensions, and rotation. Voxelize that same object independently on
each grid using exact polyhedron-cell intersection or fixed supersampling, and
record volume, centroid, surface area, and facet-normal errors.

The orientation metadata remain identical across grids; only its geometric
discretization changes. Seed mass remains separate from newly precipitated
mass, as in the current ledger.

A 20 µm crystal represented by one 20 µm voxel can match volume but cannot
resolve facet orientation. Such a grid must fail an orientation-resolution
gate even if its mass is exact. The pilot should predeclare acceptable
voxelization errors and exclude resolutions that cannot represent the smallest
measured grain or gap. Identical mass is not a grid-resolution criterion.

## Minimum experimental evidence

No project dataset was verified in this audit. The following list is a
measurement requirement, not a claim that the data are available.

### Quantities obtainable conditionally from 2D SEM or optical images

- projected crystal/grain size and aspect-ratio distributions;
- areal seed density and projected nearest-neighbor spacing;
- visible facet-edge angles and in-plane shape orientation;
- projected grain-contact length and contact frequency on exposed surfaces;
- time-resolved projected growth displacement when the same field is tracked;
- surface-film extent only if a validated optical contrast or tracer identifies
  liquid.

These are section or projection statistics. Converting them to 3D size,
orientation, contact, or gap distributions requires stereological assumptions
that must be stated and tested. Ordinary SEM shape alone does not give full
crystallographic orientation; surface EBSD could add local orientation where
applicable.

### Quantities that genuinely require 3D or equivalent volumetric evidence

- grain volume and 3D size distribution;
- full 3D grain orientation when shape does not uniquely determine it;
- grain-contact surface area and coordination number;
- intergranular gap and throat-size distributions;
- distinction between through-connected, blind, and enclosed voids;
- pore-network topology and tortuous path geometry;
- internal structural change during growth;
- volumetric liquid-film continuity, if this is a model target.

At least two growth states from the same specimen or statistically matched
specimens are required to constrain structural change. A single final image can
constrain shape but cannot identify growth kinetics or the sequence of contact
formation.

### Minimum pilot calibration set

Before a physical calibration claim, obtain:

1. a grain/seed size and spacing distribution;
2. facet shape or orientation information;
3. at least one contact/gap distribution;
4. initial and later growth states;
5. evidence defining where a surface liquid film is present, or an explicit
   statement that film accessibility remains hypothetical.

Without 3D gap/connectivity evidence, the pilot may be calibrated only to 2D
projected morphology and must not claim internal porosity prediction.

## Verification and outputs

Machine-readable outputs under phase2_ca/outputs/phase2a2r include:

- geometry_benchmarks.csv;
- selected_mask_reanalysis.csv;
- column_envelope_analytical_property.json;
- manifest.json with six verified source-mask hashes;
- SHA256SUMS.txt.

The corrected Phase 2 test suite contains 44 tests and the combined repository
suite contains 113 tests:

    python -m pytest phase2_ca/tests -q
    44 passed

    python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests phase2_ca/tests -q
    113 passed

## GO/NO-GO gate for the grain-growth pilot

**Conditional GO** for a minimal numerical pilot containing grain IDs,
orientation, one cubic anisotropy parameter, an independent film-access field,
and a contact rule. The pilot must:

- pass the five geometry benchmarks;
- preserve whole-step rollback across all added state arrays;
- conserve the existing salt budget exactly;
- represent the same continuous seeds across grids and report discretization
  error;
- retain rejected contact mass explicitly;
- report fixed-window, envelope, intergranular, enclosed, and ambiguous
  exterior volumes separately;
- avoid target porosity as an input.

**NO-GO** for claiming a physically predictive grain model until the required
2D/3D calibration data are verified, and NO-GO for Phase 2B until liquid-film
access and pore identity are independently constrained.
