# M2 critical review: interface closure, attribution, and supply limitation

Historical review of baseline commit `eb8b352`. Line references below describe
that version, not the subsequently edited source. M3 identified a local ALE
liquid/salt mesh-flux inconsistency not detected by the global audits discussed
here. The baseline numerical values and the claim that the constitutive
effects are individually weak are superseded by [RESULTS.md](RESULTS.md).
The Q-control formulation is made explicit in [M3_METHODS.md](M3_METHODS.md).

Date: 2026-11-05
Scope: audit and model-design stage only. No solver equations or production
code are changed in M2 at this stage.

## 1. Executive assessment

The existing implementation is globally conservative by construction and is a
useful research prototype. Its main baseline result is internally coherent:
large densification and a small velocity change are compatible because the
model **self-replicates the front with fixed birth porosity** and computes front
speed from interfacial supersaturation, not from the minimum porosity left in
the wake.

The present evidence does **not** identify porosity evolution as the cause of
the approximately 3% Model A--Model B difference. The current campaign mixes
five changes relative to the analytical law: transient initialization,
advection/ALE transport, the exact conservative nonlinear reaction form,
porosity evolution, and constitutive updates. B1/B2 isolates only the last
item. A controlled frozen-snapshot decomposition is required before assigning
causality.

The most important physical limitation is the interface closure, not a failed
mass balance. `epsilon_birth = 0.43` fixes how much solid accompanies every
increment of height. Material then densifies behind that newly created front,
but the densified wake cannot directly change the solid fraction of later-born
layers. This is a self-replication assumption. It may be appropriate, but it is
not independently validated.

## 2. What the current model actually solves

### 2.1 State and structural closures

The dynamic state contains cell-average porosity, concentration, height, and
five cumulative inventories
([`dynamic.py:17-35`](src/prat_growth/dynamic.py#L17-L35),
[`dynamic.py:64-70`](src/prat_growth/dynamic.py#L64-L70)). B1 holds `a_v` and
`D_eff` at reference values, whereas B2 uses
`a_v = 36(1-epsilon)/d_sub` and `D_eff = sqrt(epsilon) D_s`
([`dynamic.py:54-62`](src/prat_growth/dynamic.py#L54-L62)).

The reference and birth porosities are independently named but both equal 0.43
in the baseline
([`parameters.py:7-19`](src/prat_growth/parameters.py#L7-L19),
[`baseline.json:11-22`](configs/baseline.json#L11-L22)).

### 2.2 Independent interfacial kinetic law

The top concentration is obtained from a quadratic form of the moving-interface
salt balance. The velocity is then computed independently as

\[
V = \dot h = \beta\frac{\rho}{\rho_{cr}}k_r(C_{top}-C_{eq}).
\]

This is implemented at
[`dynamic.py:72-93`](src/prat_growth/dynamic.py#L72-L93). Model A's analytical
velocity is not inserted into Model B.

The top quadratic deliberately mixes two porosities:

- `epsilon_birth` controls the front solid mass through
  `beta rho (1-epsilon_birth) k_r`;
- the adjacent evolving `epsilon_last` controls the diffusive conductance
  `rho epsilon_last D_eff / distance`.

This hybrid is visible at
[`dynamic.py:75-81`](src/prat_growth/dynamic.py#L75-L81). It is mathematically
consistent with the chosen birth-state boundary, but it is a constitutive
choice, not a consequence of conservation.

### 2.3 Internal precipitation and porosity evolution

The signed volumetric reaction is

\[
R=\rho a_v k_r(C-C_{eq}),
\]

at [`dynamic.py:95-101`](src/prat_growth/dynamic.py#L95-L101). It enters pore
volume with `-R/rho_cr` and dissolved salt with `-R`
([`dynamic.py:178-192`](src/prat_growth/dynamic.py#L178-L192)). Because solution
mass and salt mass use the same signed reaction, their combination recovers the
paper's nonlinear concentration sink proportional to
`(1-C)(C-C_eq)`. Dissolution is allowed; it is not clipped.

### 2.4 Prescribed evaporation and reconstructed liquid supply

Evaporation is imposed as `j = rho U_evap`
([`parameters.py:21-23`](src/prat_growth/parameters.py#L21-L23)). The top liquid
flux equals new-structure mass transport plus evaporation, and the liquid flux
is integrated downward through the reaction-induced divergence
([`dynamic.py:103-111`](src/prat_growth/dynamic.py#L103-L111)). Thus the bottom
reservoir flow is an output, not a separately prescribed incompatible flux.

At the accepted B2 endpoint, the reconstructed bottom Darcy velocity is
`1.65e-7 m/s`, or 1.42 times `U_evap`. The present model always supplies that
flow because it has no pressure or permeability constraint.

### 2.5 Conservative ALE creation of new volume

Mesh motion is `x V`; the top pore-volume mesh flux injects
`epsilon_birth V`, and the top salt flux removes dissolved salt at exactly
`rho_cr (1-epsilon_birth) V`
([`dynamic.py:113-117`](src/prat_growth/dynamic.py#L113-L117),
[`dynamic.py:149-156`](src/prat_growth/dynamic.py#L149-L156)). The cumulative
new-layer and internal-solid inventories use the same two source terms
([`dynamic.py:194-203`](src/prat_growth/dynamic.py#L194-L203)). Therefore

\[
\dot M_{front}=\rho_{cr}(1-\epsilon_{birth})V
\]

is hard-wired. Agreement of `M_front` with
`rho_cr(1-epsilon_birth)(h-h0)` verifies implementation consistency, not the
physical validity of constant birth porosity.

## 3. Why strong densification changes velocity only weakly

The accepted baseline values are recorded in
[`run_summary.csv:11-12`](outputs/run_summary.csv#L11-L12).

| Quantity at 6 h | B1 | B2 |
|---|---:|---:|
| `V` (m/s) | `3.04043e-8` | `3.03811e-8` |
| Model A `V` (m/s) | `2.94034e-8` | `2.94034e-8` |
| `epsilon_min` | 0.24389 | 0.24313 |
| top-cell `epsilon` | 0.42791 | 0.42790 |
| `C_top-C_eq` | `1.1952e-5` | `1.1943e-5` |
| cumulative front solid (kg/m2) | 0.81144 | 0.81067 |
| cumulative internal solid (kg/m2) | 0.26465 | 0.26542 |

The quantitative explanation has four parts.

1. **Velocity has no direct `epsilon_min` factor.** The explicit kinetic law
   depends on `C_top-C_eq`; porosity acts only indirectly through the transport
   problem and top balance.
2. **The most densified material is mostly a wake state.** In B2 the top cell is
   still 0.4279 because fresh volume continuously enters at 0.43. Porosity does
   fall within the active layer--approximately 0.368 at 0.12 micrometres,
   0.271 at 0.56 micrometres, and 0.248 at 1.08 micrometres below the top--but
   the later wake minimum of 0.243 is not the boundary value used to create the
   next layer.
3. **B1 already includes the strongest porosity coupling.** Although B1 freezes
   `a_v` and `D_eff`, evolving `epsilon` still changes pore storage, the
   `epsilon D_eff` diffusive coefficient, ALE pore-volume flux, and liquid
   divergence. B2 adds only constitutive updates. At the top, B2 changes `a_v`
   by about +0.37% and `D_eff` by about -0.24% relative to B1, so a below-0.1%
   B1/B2 velocity difference is unsurprising.
4. **The birth rule suppresses front-memory feedback.** Every new layer again
   starts at 0.43. Densification can alter the concentration gradient through
   the reaction layer, but it cannot lower the solid fraction assigned to the
   next increment of height. This is precisely the self-replication regime.

Relative to Model A, the terminal deviations are +3.40% for B1 and +3.33% for
B2; B2 is only -0.076% below B1. These numbers establish that the *additional*
B2 constitutive update is negligible. They do not establish that porosity is
responsible for all, or none, of the remaining 3%.

## 4. Why the source of the 3% is not yet identified

The analytical verification solves a quasi-steady, no-advection, linearized,
fixed-structure problem
([`analytic.py:40-68`](src/prat_growth/analytic.py#L40-L68),
[`analytic.py:85-125`](src/prat_growth/analytic.py#L85-L125)). The dynamic model
simultaneously changes transport time dependence, advection/ALE terms, reaction
nonlinearity, porosity, and constitutive coefficients. The current tests prove
the reduced analytical discretization and short-run global conservation
([`test_analytic.py:5-17`](tests/test_analytic.py#L5-L17),
[`test_conservation.py:8-16`](tests/test_conservation.py#L8-L16)); they do not
perform a one-factor attribution.

Plausibility estimates already constrain expectations:

- after startup, the boundary-layer diffusion time is far below the 5 min
  output interval, so a persistent 3% difference is unlikely to be solely the
  equilibrium initial condition;
- `(C_top-C_eq)/(1-C_eq)` is only about `1.6e-5`, so linearization alone should
  be a very small effect;
- `Pe` is about 0.53 while `sqrt(Da)` is about 3200 at 6 h, suggesting a small
  boundary-layer advection correction;
- the accepted power-2.5 versus power-3 grid changes terminal velocity by
  0.44%, so spatial discretization is nonzero but does not explain all 3%.

These are scale arguments, not a causal decomposition. The controlled matrix
below is required.

## 5. Proposed controlled comparison matrix

All comparisons must use the same `h`, `j`, densities, bottom `C=C_eq`, and top
mass balances. The primary outcome is `V`; secondary outcomes are `C_top`,
integrated reaction, `L95`, and salt/front residuals.

### 5.1 Frozen-snapshot attribution

At selected saved states (5 min, 1 h, 3 h, 6 h), freeze `h` and the structural
fields, then solve the steady boundary-value problem. This avoids the
solid-volume inconsistency that would arise from evolving a nonzero internal
precipitation rate while artificially forcing porosity to remain constant.

| ID | Steady operator | Advection | Reaction | Porosity field | `D_eff` | `a_v` | Attribution |
|---|---|---|---|---|---|---|---|
| Q0 | paper equation | off | linearized | uniform 0.43 | reference | reference | analytical/numerical identity |
| Q1 | conservative | off | exact nonlinear | uniform 0.43 | reference | reference | linearization effect |
| Q2 | conservative | on, including mesh term | exact | uniform 0.43 | reference | reference | advection effect |
| Q3 | conservative | on | exact | saved `epsilon(z)` | reference | reference | porosity/storage effect |
| Q4 | conservative | on | exact | saved `epsilon(z)` | evolving | reference | diffusivity-only effect |
| Q5 | conservative | on | exact | saved `epsilon(z)` | reference | evolving | surface-area-only effect |
| Q6 | conservative | on | exact | saved `epsilon(z)` | evolving | evolving | combined B2 and interaction |

Use ordered differences `Q1-Q0`, `Q2-Q1`, `Q3-Q2`, and `Q6-Q3`. Report the
non-additive interaction as `Q6-Q4-Q5+Q3`. Q0 must reproduce Model A before the
other contrasts are interpreted.

### 5.2 Transient attribution

| ID | Dynamic model | Initial concentration | Purpose |
|---|---|---|---|
| T0 | B1 | uniform `C_eq` | current baseline startup |
| T1 | B1 | matched quasi-steady Q3 profile | isolate initialization transient |
| T2 | B2 | uniform `C_eq` | current full model |
| T3 | B2 | matched quasi-steady Q6 profile | initialization sensitivity with evolving coefficients |

At each output time compare the transient velocity with the matching frozen
quasi-steady solution. The transient contribution is `T-Q`, not the full
`T-Model A` difference.

### 5.3 Interpretation rule

No single contrast may be called a porosity effect unless only the porosity
field changes. Fixed-porosity transient runs with nonzero net precipitation
must not be presented as physical solid-volume models; frozen-snapshot solves
are the acceptable diagnostic alternative.

## 6. Alternative front-growth closure

### 6.1 Limitation of the current closure

The current model specifies both a constant birth solid fraction and a kinetic
front speed. The top concentration adjusts so that the corresponding salt flux
is available. This is conservative, but it assumes the morphology of every new
layer is a replica of the reference structure.

Simply setting `epsilon_birth = epsilon_top_cell` would feed a grid-dependent
cell average back into a physical boundary state, double-count densification,
and make the closure resolution-dependent. It is rejected.

### 6.2 Proposed morphology-state closure

Introduce an interfacial salt deposition flux `J_f`, an exposed crystal-area
ratio `A_f` (crystal area per projected area), and a birth solid fraction
`phi_b = 1-epsilon_b` determined by a crystal-scale morphology state:

\[
J_f = \rho k_r A_f\,[C_{top}-C_{eq}],
\]

\[
V = \frac{J_f}{\rho_{cr}\phi_b},
\]

\[
(A_f,\phi_b)=\mathcal{M}(N_f,\ell_f,\theta_f,\text{orientation},
\text{supersaturation history}).
\]

Here `N_f` is active crystal/nucleus density, `ell_f` a crystal size state, and
`theta_f` a coverage/packing state. `M` must come from a specified cubic-crystal
population model, direct simulation, or measurements; it is not a tuning
function chosen to enlarge the deviation.

The conservative coupling is then:

- top dissolved-salt ALE flux: `G_s,top = J_f`;
- new-layer solid volume rate: `phi_b V = J_f/rho_cr`;
- new-layer pore volume rate: `(1-phi_b) V`;
- internal solid volume rate in existing cells: `integral(R/rho_cr dz)`;
- dissolved redistribution: obtained from the conservative salt-flux
  divergence, never from an assigned precipitation fraction.

Therefore

\[
\frac{d}{dt}\int_0^h (1-\epsilon)\,dz
=\frac{J_f}{\rho_{cr}}+\int_0^h\frac{R}{\rho_{cr}}\,dz
\]

and total salt remains equal to initial salt plus the full bottom salt influx.

### 6.3 Additional information required

The alternative closure requires independent values or laws for nucleation
density, active surface-area ratio, crystal size and orientation distribution,
coverage-to-packing relation, burial/occlusion of crystal surfaces, and the
relation between supersaturation history and front morphology. Contact angle,
capillary connectivity, and a defensible upscaling of `beta` are also needed.
Without these inputs, the current constant-birth closure is more transparent
than an underdetermined morphology model.

## 7. Minimal hydraulic supply limitation

### 7.1 Proposed feedback

Use an explicitly exploratory permeability law, for example

\[
K(\epsilon)=K_0
\left(\frac{\epsilon}{\epsilon_0}\right)^3
\left(\frac{1-\epsilon_0}{1-\epsilon}\right)^2,
\]

and define the hydraulic resistance

\[
\mathcal R_h=R_{substrate}+\int_0^h\frac{\mu}{K(\epsilon)}\,dz.
\]

The maximum bottom liquid supply is

\[
U_{cap}=\frac{\Delta p_{available}}{\mathcal R_h}.
\]

The required `U_bottom` must still be obtained from the same liquid mass
balance used by the current model. Atmospheric demand `j_atm` becomes an upper
bound, with the actual evaporation flux determined by

\[
0\le j\le j_{atm},\quad U_{req}(j)\le U_{cap},\quad
(j_{atm}-j)(U_{cap}-U_{req})=0.
\]

This complementarity condition preserves the prescribed-flux regime whenever
supply is adequate and reduces evaporation only when the pressure budget is
exhausted.

### 7.2 Baseline diagnostic

For the current B2 endpoint, `U_req = 1.65e-7 m/s`. With the above relative
Kozeny--Carman form, `K_min/K0` is approximately 0.103 and
`integral(dz/K_rel) = 6.94e-3 m`. For an illustrative viscosity of
`1.5e-3 Pa s`, the critical reference permeability is

\[
K_{0,crit}=\frac{1.71\times10^{-12}}{\Delta p_{available}[Pa]}\;m^2.
\]

This expression is a screening relation, not a conclusion that supply is or is
not adequate. `K0`, substrate resistance, viscosity, capillary pressure, and
their uncertainty are not yet evidenced. A porosity-stop event must therefore
remain a model-domain failure, not experimental validation of pore blocking.

## 8. Numerical verification required before M2 implementation is accepted

### 8.1 Sparse implicit formulation

The current BDF call supplies neither a Jacobian nor Jacobian sparsity
([`dynamic.py:235-257`](src/prat_growth/dynamic.py#L235-L257)). Eliminating all
face liquid fluxes makes numerical differentiation appear globally coupled and
caused the high-resolution runs to become expensive.

The preferred M2 formulation is a custom backward-Euler/BDF Newton solve with
face liquid fluxes retained as algebraic unknowns. Cell pore-volume, salt, and
liquid-flux residuals are then local; the Jacobian is block-banded plus a small
number of interface rows. Assemble it with `scipy.sparse` and solve with sparse
LU. An analytic or verified finite-difference Jacobian must be tested against
directional derivatives.

### 8.2 Independent residuals

For every accepted step, record independently recomputed residuals for:

1. each cell's liquid balance;
2. each cell's dissolved-salt balance;
3. the top liquid and salt moving-interface balances;
4. `J_f-rho_cr phi_b V`;
5. new-layer solid volume versus the integrated front source;
6. global salt and solution-plus-solid mass.

The present global residuals
([`dynamic.py:296-339`](src/prat_growth/dynamic.py#L296-L339)) reuse the same
flux and inventory definitions as the ODE. They are necessary but not an
independent local verification.

### 8.3 Acceptance tests

- maximum normalized local balance residual below `1e-8`, with both `L_inf`
  and volume-weighted `L1` reported;
- moving-interface and front solid-volume residuals below `1e-10` relative to
  their characteristic fluxes;
- global salt and bulk-mass residuals below `1e-8`;
- at least 8 cells across `delta_C` for the formal convergence run, with a
  second criterion based on `L95`;
- terminal height, velocity, `epsilon_min`, `L95`, and internal/front partition
  changing by less than 0.5% between the two finest grids;
- independent maximum-step refinement on the accepted finest grid;
- matched equilibrium and analytical/quasi-steady initial profiles converging
  to the same post-startup trajectory;
- physical root required for `C_top`; the fallback to the closest unphysical
  quadratic root at
  [`dynamic.py:87-91`](src/prat_growth/dynamic.py#L87-L91) must become a hard
  failure with diagnostics;
- cached runs may be reused only after a full configuration and code hash
  match. Current reuse checks only the case filename
  ([`campaign.py:154-166`](src/prat_growth/campaign.py#L154-L166)), which can
  silently confound comparisons after configuration changes.

## 9. Priority and decision gate

1. Implement the frozen-snapshot Q0--Q6 solver and the T0--T3 initial-condition
   controls before changing the front closure.
2. Add independent local/interface residuals and hash-checked provenance.
3. Replace the dense numerical BDF pathway with the sparse residual/Jacobian
   formulation and repeat convergence.
4. Add the hydraulic supply diagnostic with uncertainty ranges for `K0` and
   available capillary pressure.
5. Implement the morphology-state front only after at least one independent
   constraint on `A_f` or `phi_b` is available.

The M2 scientific decision is therefore **not yet whether internal
precipitation changes the regime**. The immediate decision is whether the 3%
difference survives controlled attribution and independent residual checks. A
fundamental regime change can be claimed only if the effect persists after
those controls and is tied to a separately evidenced front-morphology or
hydraulic-feedback mechanism.
