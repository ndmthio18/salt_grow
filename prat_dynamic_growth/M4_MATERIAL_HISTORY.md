# M4 fixed-material histories and travelling densification profile

Date: 2026-11-05

Reference: M3 commit `a37759f90f0aa92b3fb37df7b4106cf711da2d02`

Decision: **A — the calculations support an approximately stable moving
reaction–densification profile, with explicit numerical qualifications.**

This is a statement about the present continuum/ALE model. It is not physical
validation of pore-scale salt filling, the constant birth porosity, or the
assumed unlimited liquid supply.

## What was tracked

The solid skeleton has no translation or deformation velocity in the current
model. A material marker is therefore a fixed laboratory coordinate `z0`.
For material born after the initial state,

```text
h(tb) = z0
age = t - tb
x(t) = z0/h(t)
d(t) = h(t) - z0.
```

Three coordinates answer different questions:

- fixed `z=z0` follows one material location;
- fixed `x=z/h(t)` follows the ALE mesh and changes physical material;
- fixed `d=h(t)-z` follows the front-frame profile and changes material.

No liquid Darcy velocity was used to move a solid marker. The synthetic test
in `tests/test_material_tracking.py` verifies the fixed-`z` trajectory, birth
time and exact reaction integral, and demonstrates why a fixed grid index is
not a material trajectory.

The primary reconstruction is bounded linear interpolation of finite-volume
cell averages in the actual coordinates `z=h(t)x`, including the prescribed
birth state and solved top concentration. A piecewise-constant reconstruction
is retained as a deliberately lower-order sensitivity. No unconstrained
high-order interpolation is used.

## Time resolution and restart checks

At 6 h,

| Scale | Value |
|---|---:|
| Reference reaction thickness `delta_ref` | 0.513631 micrometres |
| `tau_delta = delta_ref/V` | 17.0881 s |
| `tau_95 = L95/V` | 32.6611 s |

Newborn cohorts use birth times 5 min, 1 h, 3 h and 5 h for B1 and B2.
Each checkpoint window covers 300 s, over 17 `tau_delta`. Fine output spacing
is 0.170843 s and the comparison spacing is 0.341686 s. These are dense-output
sample intervals, not accepted BDF steps. Regular windows use a 5 s maximum
step; the B2–640 1 h cohort was repeated with a 1 s maximum step.

The 5 s/1 s comparison changes final porosity by `1.1e-15`, cumulative local
reaction by `6.23e-12` relative, and the 95% densification time by zero at the
reported sampling resolution. Halving output spacing changes cumulative
reaction by at most `2.87e-5` relative and the 95% time by at most 0.524%.
All window warning lists are empty. Restart endpoints agree with the archived
long runs; exact values are recorded in `manifest.json`.

Initial material at `h0-delta_ref`, `h0-3 delta_ref` and `h0/2` is a secondary
control. It uses the M3 T1/T3 matched quasi-steady concentration initialization
to avoid repeating the unrelated millisecond concentration startup. M3 showed
that equilibrium and quasi-steady initializations become structurally
indistinguishable after startup. These initial markers are not included in the
newborn-cohort convergence claim.

## H1: the global minimum is stable, but its material identity is not

For B2–640, `epsilon_min` is 0.23609365 at 5 min and 0.23612005 at 6 h, an
increase of only `2.64e-5`. B1 changes from 0.23777016 to 0.23779349.
The physical position of the minimum moves substantially:

| B2 time | Distance behind front | Approximate birth time of minimum location |
|---|---:|---:|
| 5 min | 4.92 micrometres | 136 s |
| 1 h | 75.60 micrometres | 1085 s |
| 3 h | 244.70 micrometres | 2658 s |
| 5 h | 424.72 micrometres | 3869 s |
| 6 h | 516.06 micrometres | 4430 s |

Thus H1 is true, but the minimum does not remain attached to one fixed marker.
The small rise is consistent with position migration, finite-volume sampling
and the weak secular rise of later cohort plateaus. Signed dissolution flux is
about `1e-13 kg m-2 s-1` or smaller and does not support dissolution as the
cause.

Volume-weighted regional averages use exact cell overlap at `h0`, rather than
cell-centre classification. For B2–640, the initial-region mean changes from
0.42987 at 5 min to 0.42740 at 6 h, the newborn-region mean from 0.25007 to
0.24015, and the whole-domain mean from 0.42826 to 0.35369. A stable minimum
therefore does not imply a stable mean structure.

## H2 and H3: fixed newborn material reaches a common finite-window platform

### B2

| Birth time | `z0` (mm) | Porosity at age 300 s | 95% densification time |
|---|---:|---:|---:|
| 5 min | 1.009010 | 0.2360946 | 32.802 s |
| 1 h | 1.108188 | 0.2361142 | 32.802 s |
| 3 h | 1.324584 | 0.2361536 | 32.802 s |
| 5 h | 1.540990 | 0.2361894 | 32.802 s |

### B1

| Birth time | `z0` (mm) | Porosity at age 300 s | 95% densification time |
|---|---:|---:|---:|
| 5 min | 1.009029 | 0.2377705 | 40.034 s |
| 1 h | 1.108425 | 0.2377881 | 40.034 s |
| 3 h | 1.325295 | 0.2378235 | 40.034 s |
| 5 h | 1.542175 | 0.2378556 | 39.864 s |

These values establish H2 for newborn material: fixed `z0` markers actually
move from the prescribed 0.43 birth porosity to the reported finite-window
platform. They also support H3: within each closure, the four 640-cell terminal
values span less than `1e-4`, and histories plotted against material age
visually overlap.

The word “terminal” here means age 300 s. It is not an infinite-time limit.
When local reaction falls below 0.1% of its cohort maximum, at about 75.7 s
for B2 and 93.5 s for B1, subsequent reconstructed porosity changes by only
about `1.0e-4` and `1.4e-4`, respectively. Because reaction is then negligible,
that small drift is treated as reconstruction/numerical mixing, not continuing
physical precipitation.

B2 densifies faster and farther than B1: approximately 32.8 versus 40.0 s,
and about 0.00167 lower porosity at matched material age. This is a model
closure difference arising from evolving `a_v` and `D_eff`, not experimental
evidence that either closure is physically correct.

Initial material does not share the full newborn trajectory. At age 300 s,
the positions initially at `h0-delta_ref` and `h0-3 delta_ref` retain
porosities 0.3468 and 0.3991 in B2 (0.3439 and 0.3967 in B1), while the
mid-height location remains 0.43. These locations began inside or behind the
initial reaction zone and never experienced the missing part of a 0.43-to-wake
birth path. They are not counterexamples to H3, whose population is newly
formed material.

## Material reaction balance and its numerical limit

For fixed material, the continuum relation is

```text
d epsilon(z0,t)/dt = -R(z0,t)/rho_crystal.
```

The comparison was independent in the relevant sense: `epsilon_field` comes
from the stored porosity field, while `epsilon_reaction` integrates the local
reaction reconstructed from concentration and `a_v`. Reaction was never
inferred from porosity differences. Positive precipitation and negative
reaction were integrated separately.

For newborn cohorts, the final 300 s discrepancy is at most 0.147% of the
total porosity loss. The maximum discrepancy during cell crossing is larger:
1.24–1.57% on 640 cells and up to 2.10% on 480 cells. It peaks near material
age 10 s, while the marker crosses the strongly varying, differently sized top
cells. Piecewise-constant reconstruction is worse (about 2.4–3.2%), whereas
the final linear-versus-piecewise porosity difference is below `2.1e-8`.

Therefore the proposed 1% material-balance target is met at the window endpoint
but **not at every material age**. Output integration and BDF step checks are
far smaller than this defect. Its decrease under grid refinement and its
association with cell crossings point to first-order ALE transport,
cell-average reconstruction and numerical mixing. It is not interpreted as a
new physical term.

## Travelling reaction and densification profiles

On a common physical distance axis, the reaction profiles at 5 min, 1 h, 3 h,
5 h and 6 h differ from the 6 h profile by at most 0.058% in B2 and 0.031% in
B1. Growth speed changes only 0.0146% (B2) and 0.0132% (B1) from 5 min to 6 h.
This supports a nearly stationary reaction profile and nearly constant speed.

Within the first 5 micrometres, the 5 min porosity profile differs from the 6 h
profile by at most 0.000495 (B2) and 0.000406 (B1). At 5 min the newborn layer
is only about 9 micrometres thick, so a 10-micrometre comparison crosses into
the initial structure and produces a large, non-comparable difference. A wake
platform is not declared at 5 min. From 1 h onward it is detected, and its
change from 1 h to 6 h is only `9.22e-5` in B2 and `8.27e-5` in B1.

The explicit front-frame time derivative is below `7e-6` of the reaction
source scale after 1 h. It is not assumed zero a priori. Nevertheless, a
pointwise finite-difference check of

```text
partial_t epsilon|d + V partial_d epsilon = -R/rho_crystal
```

has a 7–12% maximum normalized residual on the reconstructed profile, and the
integral form retains a 1.2–2.1% transient/grid defect. Both improve from 480
to 640 cells but do not meet a strict pointwise closure standard. This limits
the spatially local claim; it does not erase the directly observed overlap of
profiles and fixed-material histories.

The more robust integrated wake check is consistent after 1 h. For B2 at 6 h,

```text
G_internal = 1.26403e-5 kg m-2 s-1
chi_internal_inst = 0.253726
epsilon_wake measured = 0.23620585
epsilon_wake predicted = 0.23620557
difference = 2.85e-7.
```

B1 gives a 6 h difference of `1.04e-6`. At 5 min the B2 wake relation misses
by 0.00228 and no platform is detected, so the steady relation is correctly
classified as inapplicable then. The cumulative M3 internal fraction was not
inserted into this relation.

## Numerical acceptance

| Diagnostic | Result | Status |
|---|---:|---|
| 480/640 final cohort porosity difference | <= 0.000254 | pass (target 0.001) |
| Output-halving reaction difference | <= 2.87e-5 relative | pass (target 1%) |
| Output-halving 95% time difference | <= 0.524% | pass (target 1%) |
| Representative 5 s/1 s integration-step porosity difference | `1.1e-15` | pass |
| Final material-balance mismatch | <= 0.147% of loss | pass (target 1%) |
| Maximum all-age material-balance mismatch | <= 2.10% of loss | **does not pass** |
| Stable-profile pointwise differential residual | 7–12% maximum | qualification |
| Newborn history agreement across birth times | < 0.001 | pass |

The accepted conclusion is therefore not “all local diagnostics are converged.”
It is that the stable moving profile and common newborn history are supported
at the observable level, while pointwise material conservation during cell
crossing remains limited by the first-order finite-volume representation.

## Answers to the ten required questions

1. **Is 0.236 only a stable global minimum?** No. H1 is true, but fixed-marker
   evidence separately shows a material plateau. The minimum itself migrates.
2. **Does fixed material go from 0.43 to about 0.236?** Yes for B2 newborn
   material; B1 reaches about 0.238. Initial material follows partial histories.
3. **Do different birth cohorts reach one platform?** Approximately yes within
   each closure; 640-cell 300 s values span less than `1e-4`.
4. **How long is main densification?** About 32.8 s in B2 and 40.0 s in B1 for
   95% of the observed drop.
5. **Does porosity change after leaving the reaction zone?** Only by about
   `1.0e-4` (B2) or `1.4e-4` (B1) after the local reaction is below 0.1% of its
   maximum; this remainder is assigned to numerical reconstruction/mixing.
6. **Does local reaction explain porosity loss?** At the endpoint to better
   than 0.15% of the loss, but the transient cell-crossing error reaches 2.10%.
7. **Is there a stable moving profile?** Approximately yes. Reaction profiles
   and speed are stable early; a complete porosity wake is identifiable from
   about 1 h onward.
8. **Is the conclusion robust?** Final cohort values pass the two-grid,
   output-halving and representative time-step checks. Pointwise differential
   balance remains only partially converged and is explicitly limiting.
9. **How do B1 and B2 differ?** B2 densifies about 7 s faster and reaches a
   porosity about 0.0017 lower; both exhibit the same qualitative moving state.
10. **Is this physical validation?** No. It is an internally stable state of
    the present computational closure. Crystal-scale validity, birth morphology
    and hydraulic sustainability remain untested.

## Reproduction and artifacts

Run:

```powershell
$env:PYTHONPATH = (Join-Path (Resolve-Path '.').Path 'src')
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
python -m prat_growth.m4_analysis
python -m pytest -q
```

Primary evidence is in:

- `outputs/m4/material_histories.csv`
- `outputs/m4/cohort_summary.csv`
- `outputs/m4/spatial_statistics.csv`
- `outputs/m4/material_balance.csv`
- `outputs/m4/travelling_profile_checks.csv`
- `outputs/m4/travelling_profile_summary.csv`
- `outputs/m4/manifest.json`

The manifest records the exact reference commit, M3 input hashes, parameters,
solver statistics, accepted-step ranges, restart comparisons and diagnostic
thresholds. M3 data were read as immutable archives and were not overwritten.
