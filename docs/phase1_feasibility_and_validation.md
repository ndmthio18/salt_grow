# Phase 1 feasibility and validation

## Decision summary

**Phase 1 result: GO for the independent evaporation/mass-allocation module;
NO-GO for automatic progression to Phase 2 pending review.** The implemented
module closes the prescribed steady mass budget and conservatively carries
sub-voxel salt mass. This is numerical verification of a deliberately simple
accounting model, not validation of salt-growth physics or morphology.

The original CUDA CA was audited at reference commit
`995f197192988483140e6f7e773c496aa16cf76b` (retrieved 2026-11-06). No explicit
`LICENSE`, `COPYING`, or `NOTICE` file was present. The reference checkout was
therefore kept outside this repository, and none of its source was copied.

## Evidence boundaries

| Claim | Evidence level |
|---|---|
| Phase 1 equations and unit handling are implemented | Executed CPU tests |
| Mass-to-voxel rounding conserves the supplied mass budget | Executed CPU tests |
| Reference RH trends and tables are reproducible | Executed script and inspected plot |
| Original CA is intended to conserve occupied lattice sites after initialization | Source-code audit only |
| Original GPU example runs on this machine | **Not established**; CUDA baseline not run |
| The illustrative transfer coefficient or fixed area represents an experiment | **Not established** |
| The resulting voxel count predicts morphology | **Not established**; spatial placement is deferred |

## Audit of the reference CA

### Representation and update rules

- The lattice is a three-dimensional `int32` array: `0` is empty, `1` is a
  mobile particle, and `2` is a crystalline particle.
- Each CA iteration first applies a synchronous attachment/detachment kernel,
  then performs diffusion. A mobile site changes `1 -> 2` on attachment and a
  crystalline site changes `2 -> 1` on detachment; the site remains occupied.
- A mobile particle attempts one nearest-neighbor hop per diffusion sweep. The
  direction is selected from the six face directions, with optional independent
  x/y/z bias. A hop is accepted only if the target is empty. The 64 shuffled
  `(cycleX, cycleY, cycleZ)` offsets form a 4-by-4-by-4 coloring intended to
  avoid simultaneous conflicts during the in-place CUDA update.
- The rule index is the seven-trit base-3 encoding

  `729*self + 243*right + 81*down + 27*left + 9*up + 3*back + front`,

  which spans `0..2186`. `Ru[index]` is the transition probability. In the
  default ruleset, attachment depends on crystalline-neighbor count through
  `Pn`, `Ps`, `Pk`, and unity for more highly coordinated sites; detachment
  uses `Pne`, `Pse`, and `Pke` with the reverse transition.

### Conservation, boundaries, and initialization

- By inspection, both rule transitions preserve the number of occupied sites,
  and diffusion moves one occupied state from source to destination. There is
  no active source/supply kernel. The allocated `depositions_gpu`, `flux`, and
  `p0` variables do not implement a supply mechanism.
- This conservation statement is conditional on the intended collision-free
  CUDA schedule and was not empirically verified here because the GPU baseline
  could not run.
- Each domain face can independently be a blocking wall; otherwise that
  direction is periodic. Across a wall, neighbor reads return empty and hops
  are blocked.
- Initialization modes are centered single seed, multiple random seeds, flat
  bottom, stepped bottom, a pair of seeds, sinusoidal mounds, random occupancy
  only, and user-provided Python code.
- The code contains iteration counts and lattice coordinates, but no physical
  timestep, voxel length, particle mass, density, or calibration mapping. The
  probabilities are dimensionless algorithmic inputs.

### Material issues for the intended application

1. **License gate.** No explicit license was found. Redistribution or a derived
   codebase cannot be assumed permissible.
2. **No physical mass supply.** Closed particle redistribution cannot represent
   continuously supplied saturated brine without a separate mass ledger and a
   controlled creation mechanism.
3. **No physical scale.** CA steps and voxels have no mapping to seconds, metres,
   or kilograms.
4. **Grid divisibility is unchecked.** Diffusion launches with `size//4` while
   its indexing macros assume `4*gridDim`; dimensions not divisible by four can
   leave remainder planes outside the diffusion sweep.
5. **Partial stochastic reproducibility.** `rng_seed` seeds NumPy initialization,
   but the PyCUDA XORWOW generator is constructed without that seed. A complete
   stochastic run is therefore not pinned by the documented configuration.
6. **Misleading internal scalar.** `gpuarray.sum(atoms_gpu)` weights crystalline
   states twice because their label is `2`; it is not a particle count. Separate
   mobile/crystalline counts in saved snapshots are the appropriate conservation
   diagnostic.

Conceptually reusable elements are the discrete state convention, synchronous
local-rule evaluation, base-3 neighborhood lookup, boundary configuration, and
separate mobile/crystalline diagnostics. The empirical rule probabilities,
unscaled iteration time, and closed particle inventory are not yet physical
salt-growth closures.

## GPU baseline status

The original GPU example was **not executed**. The host exposes only an Intel
Arc graphics adapter; `nvidia-smi`, CUDA environment variables, PyCUDA, CuPy,
and Numba were absent. Consequently there are no defensible initial/final
particle counts, runtime, conservation error, or morphology output for the
original model on this host. A CUDA-to-CPU rewrite was intentionally not made.

A later Phase 1.1 remote-resource check found an NVIDIA RTX 2060 on `frodo`,
but that execution environment lacked both the CUDA toolkit and PyCUDA required
by the audited reference implementation. The hardware observation therefore
does not change the baseline status: the original CA still has not been run,
and no GPU morphology or conservation result is claimed.

## Phase 1 physical model

The implemented budget is

\[
c_{v,i}=a_w c_{v,sat}(T),\qquad
c_{v,\infty}=RH_\infty c_{v,sat}(T),
\]

\[
J_w=k_v\max(c_{v,i}-c_{v,\infty},0),\qquad
\Delta M_w=J_w A_{evap}\Delta t,
\]

\[
\Delta M_s=\frac{w_s}{1-w_s}\Delta M_w,\qquad
\Delta V_s=\frac{\Delta M_s}{\rho_s}.
\]

The saturation pressure uses the Buck (1981) liquid-water relation, followed by
the ideal-gas conversion `c = p M_w/(R T)`. At 298.15 K the implemented value is
`0.0230266116 kg/m3`, consistent with approximately `23 g/m3`. A validated
saturation concentration can instead be supplied explicitly.

The proof-of-concept assumptions are: isothermal conditions; continuously
supplied saturated NaCl brine; pre-existing seeds; stationary liquid inventory;
no dissolved-salt export; fixed evaporation area; no morphology feedback on
evaporation; sufficient wetting; and no primary nucleation, condensation,
deliquescence, or dissolution. These are imposed assumptions, not measured or
validated properties of the developing porous salt.

For `RH > a_w`, positive evaporation is zero and the result is flagged
`above_equilibrium_unmodeled`; the code does not reinterpret this as a valid
condensation calculation.

## Mass-to-voxel ledger

For cubic voxels,

\[
m_{voxel}=\rho_s\Delta x^3,\qquad
N_{new}=\left\lfloor\frac{\Delta M_s+M_{residual}}{m_{voxel}}\right\rfloor.
\]

Only `N_new*m_voxel` is allocated; the remainder is carried to the next call.
The stateful ledger verifies

\[
M_{allocated}+M_{residual}=M_{cumulative\ budget}
\]

within floating-point tolerance. Initial seed voxels are stored separately and
never included in newly precipitated mass. No spatial voxel creation or CUDA
kernel modification is present.

## Reference computation

Inputs were `T=298.15 K`, `a_w=0.753`, `w_s=0.264`,
`rho_s=2160 kg/m3`, `k_v=0.005 m/s`, a `64^3` grid with `10 um` voxels,
`A_evap=(64*10 um)^2=4.096e-7 m2`, and `dt=60 s`. Thus
`w_s/(1-w_s)=0.3586956522` and `m_voxel=2.16e-12 kg`.

| RH | Regime | `J_w` [kg m-2 s-1] | water/60 s [kg] | NaCl/60 s [kg] | new voxels | residual [kg] |
|---:|---|---:|---:|---:|---:|---:|
| 0.300 | evaporation | 5.21553e-5 | 1.28177e-9 | 4.59765e-10 | 212 | 1.84462e-12 |
| 0.500 | evaporation | 2.91287e-5 | 7.15866e-10 | 2.56778e-10 | 118 | 1.89804e-12 |
| 0.700 | evaporation | 6.10205e-6 | 1.49964e-10 | 5.37914e-11 | 24 | 1.95145e-12 |
| 0.753 | equilibrium | 0 | 0 | 0 | 0 | 0 |
| 0.800 | above-equilibrium, unmodeled | 0 | 0 | 0 | 0 | 0 |

At constant conditions, the 24 h cumulative NaCl masses are
`6.62061e-7`, `3.69760e-7`, and `7.74597e-8 kg` for RH 0.30, 0.50, and
0.70, respectively. These are model outputs under an illustrative, uncalibrated
`k_v`; they are not experimental predictions.

Generated evidence is in `phase1_mass_budget/outputs/reference/`:

- `rh_sweep.csv`
- `cumulative_salt_mass.csv`
- `reference_summary.json`
- `phase1_reference_results.png` and `.svg`
- `SHA256SUMS.txt`

The plot was visually inspected: axes, units, equilibrium marker, monotonic RH
trend, and zero high-RH curves are rendered correctly.

## Verification actually executed

```text
python -m pytest phase1_mass_budget/tests -q
34 passed in 0.10s

$env:PYTHONPATH=(Resolve-Path '.\\prat_dynamic_growth\\src').Path
python -m pytest prat_dynamic_growth/tests phase1_mass_budget/tests -q
50 passed in 12.97s

python phase1_mass_budget/run_reference_case.py
completed successfully and wrote all reference outputs
```

The tests cover positive/equilibrium/high-RH behavior, monotonicity, zero area
and zero transfer coefficient, analytical water/salt balances, timestep
partition independence, invalid units/ranges, multiple voxel sizes, sub-voxel
carry-forward, seed exclusion, and cumulative closure.

## Scientific limitations and Phase 2 gates

Phase 1 establishes only a scalar supply budget. It cannot identify where salt
appears, whether a site is wet, whether ion transport sustains the assumed
supply, or whether the fixed area and `k_v` remain valid as morphology evolves.
The equality between incoming dissolved salt and precipitation follows from the
stationary-inventory/no-export assumptions; it is not independently predicted.

Before Phase 2:

1. obtain an explicit license or continue with a clean, independent
   implementation of only documented concepts;
2. execute a pinned 32-cubed or 64-cubed CUDA baseline on compatible hardware,
   with initial/final occupied counts and morphology artifacts;
3. make CUDA stochastic seeding complete and enforce dimensions divisible by
   four before using the reference implementation quantitatively;
4. calibrate or bound `k_v`, `a_w`, `w_s`, effective area, and their temperature
   dependence with experimental provenance and uncertainty;
5. define a deterministic spatial-allocation policy, wetting eligibility, seed
   bookkeeping, and a coupled global conservation test before creating voxels;
6. retain a strict distinction between mass closure, numerical convergence,
   morphology plausibility, and experimental validation.

No permeability model, vapor-field solver, PNM coupling, morphology feedback,
or advanced growth rule has been implemented.

## Reference

Buck, A. L. (1981), *New Equations for Computing Vapor Pressure and Enhancement
Factor*, Journal of Applied Meteorology 20, 1527-1532,
`doi:10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2`.
