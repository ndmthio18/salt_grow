# M3 numerical attribution: equations and interpretation

This study implements `instruction/3.md` with unchanged physical parameters.
Morphology-dependent birth rules and hydraulic feedback are outside this study.
The baseline is commit `eb8b352`; its results remain in `outputs/` as historical
artifacts. Corrected calculations are stored separately in `outputs/m3/`.

## Frozen means instantaneous structure, not zero pore-volume rate

Let `E_i = epsilon_i h Delta_x_i`, `Q_i = rho C_i E_i`, and `R_i` be internal
precipitation per bulk volume. Let `H = -w epsilon_upwind` be the material
pore-volume mesh flux, with `H_top = -V epsilon_birth`. The equations are

```
dE_i/dt = H_left - H_right - R_i volume_i/rho_cr
F_left - F_right = (1-rho/rho_cr) R_i volume_i
L = F + rho H
rho dE_i/dt = L_left - L_right - R_i volume_i
dQ_i/dt = S_left - S_right - R_i volume_i
```

The physical liquid flux is `F`; the liquid flux relative to the mesh is `L`.
The top conditions are

```
J_f = rho_cr (1-epsilon_birth) V
V = beta rho/rho_cr k_r (C_top-C_eq)
F_top = j + [rho epsilon_birth + rho_cr (1-epsilon_birth)] V
L_top = j + J_f
S_top = J_f = C_top L_top - rho epsilon D_eff grad(C)_top.
```

In a frozen solve, `h` and `epsilon` are **coefficient snapshots**. The pore-rate
equation above is retained. Set only the instantaneous ALE concentration rate
to zero, `dC_i/dt = 0`, so that `dQ_i/dt = rho C_i dE_i/dt`. Subtraction gives

```
0 = D_left - D_right
    + L_left (C_face,left-C_i) - L_right (C_face,right-C_i)
    - (1-C_i) R_i volume_i,
D_face = -rho (epsilon D_eff)_face grad(C)_face.
```

This is not a permanent steady solid structure. Setting `dE/dt = 0` instead
would change liquid divergence to `F_left-F_right = R volume`, remove the
precipitation-induced pore storage, and answer a different question. The
frozen result is also coordinate-specific: `dC/dt=0` is imposed at fixed ALE
coordinate, not fixed laboratory distance. The transient remainder explicitly
measures this approximation against the full dynamic solution.

Q0 and Q1 omit the two **concentration-advection** terms from the reduced
equation. They retain the mass-balance-derived sink and interfacial closure;
they do not set the liquid supply to zero. Accordingly they are asymptotic
operator controls, not full stationary conservative flows with zero velocity.

## Controls and additive accounting

Q0 uses uniform porosity and reference coefficients, replaces `1-C` by
`1-C_eq`, and linearizes the top condition to
`D_top = (1-C_eq) J_f - C_eq j`. It must converge to the Prat analytical law.
Q1 restores both nonlinear concentration factors. Q2 adds concentration
advection including mesh motion. Q3 uses the saved B2 porosity. Q4 changes
only diffusivity, Q5 only specific area, and Q6 changes both. All share the
same height, grid, evaporation, densities, reservoir, and birth porosity.

Two extra controls prevent top-boundary nonlinearity being mislabeled as
reaction nonlinearity: Q0R restores only the nonlinear sink, while Q0B
restores only the exact top condition. The underlying kinetic `R(C)` is
linear; the "reaction nonlinearity" refers to the reduced sink `(1-C)R`.

Every contribution is divided by the **same analytical velocity**, yielding
percentage points of `(V-V_A)/V_A`. The telescoping sum is

```
T2-A = (Q0-A)
       + (Q0R-Q0) + (Q0B-Q0) + (Q1-Q0R-Q0B+Q0)
       + (Q2-Q1) + (Q3-Q2)
       + (Q4-Q3) + (Q5-Q3) + (Q6-Q4-Q5+Q3)
       + (T2-Q6).
```

Q0-A is a discretization error. Remaining contrasts are conditional,
ordered operator effects, not unique path-independent causal effects.
The diffusivity/area interaction is reported explicitly; interactions of
these updates with porosity are conditioned on the saved porosity field.
The field itself is the result of the full B2 history. Q3 is therefore not
the B1 dynamic trajectory. T0/T1 use B1 and T2/T3 use B2; their own matched
frozen snapshots are compared separately.

## Correction to the original discretization

The original salt ALE flux used the arithmetic mean of adjacent porosities
for `-rho epsilon w`, whereas the pore-volume equation used upwind porosity.
Consequently its salt advection did not use the liquid flux obtained from
its own local pore-volume balance. Global inventory closure alone did not
expose this mismatch. Both now use the same `H` (upwind for either sign of
mesh velocity). A rerun at the original baseline grid separates this code
correction from subsequent refinement. It also contains a numerically stable
top-root evaluation; these two changes are not silently credited to physics.

The top quadratic is solved for excess concentration with a cancellation-safe
formula. Absence of a root in `[0,1]` raises an exception with state diagnostics.
No fallback returns an unphysical concentration.

## Verification and provenance

`verification.py` independently rebuilds face diffusion, liquid supply,
upwind mesh transport, interface salt flux, and cell balances. It does not
call the model's `flux_state`. It compares these with derivatives of the
stored primitive variables and geometric volume. These are semidiscrete
operator checks, not an independent experiment or proof of time integration.
Time-step refinement and global inventories provide complementary checks.
An additional derivative audit integrates 0.5 s windows from the four B2
snapshots at tighter tolerance, then uses fourth-order central differences
of the dense trajectory at offset 0.25 s. Spacings are 0.001, 0.01 and 0.05 s.
Its strict gates are reported separately and are not all passed; height
differentiation error is strongly amplified in the top diffusion condition.
`RESULTS.md` documents this limit rather than treating every audit as passed.
Residuals use the larger of evaporation/front mass flux as a fixed physical
scale; weighted L1 uses cell volume weights. The top concentration difference
is evaluated in excess variables to avoid subtracting two nearly equal
absolute mass fractions. A deliberate derivative perturbation test verifies
that the local audit detects a damaged balance.

Grid sequences use the same power-3 mesh family and 160, 240, 320, 480, 640
cells. Temporal comparisons use maximum steps 120, 60, 30 seconds and identical
tolerances. L95 is interpolated through cumulative piecewise-constant cell
reaction, rather than snapped to the first crossing face. Both finest grids
repeat all Q contrasts, bounding error in attribution as well as total speed.
All accepted dynamic cases must reach 6 hours without a stop event.

Cache reuse requires the SHA-256 of the full dataclass configuration,
closure and initialization, the names and bytes of all package Python files,
and the saved NPZ bytes. Legacy caches lack these identities and are never
silently reused as current calculations. Archival reading remains possible.
Quasi-steady initialization uses its actual initial dissolved-salt inventory;
it does not subtract an equilibrium inventory.

These calculations test this continuum closure internally. They do not
validate crystal-scale homogenization, constant birth morphology, or unlimited
liquid supply. No parameter was tuned to increase the difference from Prat.
