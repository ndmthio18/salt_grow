# Results

## Executive finding

The analytical law is numerically reproduced under its own assumptions. In the
fully conservative evolving-porosity baseline, the post-startup interface
velocity is about 3.1-3.3% higher than Prat's analytical prediction, but
updating `a_v` and `D_eff` (B2 versus B1) changes terminal velocity by only
-0.08% and terminal height by -0.04%. The present closure therefore produces a
small modified growth coefficient, not a demonstrated mechanism change.

This conclusion is numerical, not physically validated: the predicted active
reaction layer is smaller than the estimated crystal size.

## Verification and numerical acceptance

- Reference groups at `h0 = 1 mm`: `Pe = 0.31575`,
  `Da = 3.7905e6`, `Da_top = 5264.59`, `Da* = 7.31193`, and
  `delta_C = 0.5136 micrometres`.
- The 1280-cell Stage A finite-volume solution reproduces the general
  `tanh(sqrt(Da))` velocity with relative error `6.83e-6`.
- At baseline `Da`, the general and large-`Da` analytical velocities agree to
  displayed precision: `2.94034e-8 m/s` (`2.54045 mm/day`).
- On the dynamic grid, `power=2.5` versus the more strongly clustered
  `power=3.0` changes terminal height by 0.15% and terminal velocity by 0.44%.
  The accepted mesh has seven cells across final `delta_C`.
- `max_step = 120, 60, 30 s` produces negligible differences relative to the
  1% observable gate.
- Accepted B2 normalized residuals at 6 h are `-1.21e-11` for total salt and
  `-1.49e-12` for solution-plus-solid mass.

## Baseline comparison at 6 h

| Quantity | Model A | B1 | B2 |
|---|---:|---:|---:|
| Height (mm) | 1.63511 | 1.65603 | 1.65541 |
| Velocity (mm/day) | 2.54045 | 2.62693 | 2.62493 |
| Minimum porosity | constant 0.43 | 0.24389 | 0.24313 |
| `L95` (micrometres) | not an inventory | 1.2427 | 1.2423 |

The equilibrium initial condition makes the Model B interface velocity nearly
zero at `t=0`. By 5 min, B2 reaches `2.6183 mm/day`, 3.06% above Model A. The
difference then persists and slowly grows to 3.33% relative to Model A (3.22%
under the specified error denominator) at 6 h. It is therefore a short startup
transient followed by a modest persistent coefficient difference.

## Location and motion of precipitation

After startup, B2's `delta_C` remains about 0.511-0.512 micrometres and `L95`
remains about 1.01-1.31 micrometres when expressed as distance below the moving
top. At 6 h, `L95/delta_C = 2.43`. This is consistent with a nearly self-similar
reaction zone that travels with the interface and leaves a densified wake.

The final newly formed region has mean porosity 0.251 and minimum porosity
0.243, a maximum local reduction of 0.187 from the birth value (43.5% relative).
In the original `z < h0` material, only the upper part is traversed: its minimum
porosity is 0.339 and its volume-weighted mean remains 0.424.

## Actual salt partition

At 6 h, B2 contains `0.8107 kg/m2` of cumulative new-layer solid and
`0.2654 kg/m2` of cumulative internal solid; the dissolved inventory is
`0.1868 kg/m2`. Internal precipitation is therefore 24.7% of the cumulative new
solid added after `t=0`. The corresponding analytical complement
`1 - sqrt(Da*)/(1+sqrt(Da*))` is about 26.9%. Their numerical proximity does not
make them identical: the former is a time-integrated conservative inventory,
whereas the latter is an analytical flux partition under fixed properties.

## Sensitivities

| Case | End/stop | Velocity (mm/day) | `epsilon_min` | `L95` (um) | Resolution |
|---|---:|---:|---:|---:|---:|
| `d_sub=10 um` | 6 h | 2.378 | 0.163 | 0.638 | 6 cells/delta |
| `d_sub=100 um` | 6 h | 2.963 | 0.330 | 2.681 | 9 cells/delta |
| `U_evap=0.1 cm/day` | 6 h | 0.262 | 0.242 | 1.040 | 8 cells/delta |
| `k_r=2.3e-6 m/s` | stop 220 s | 0.860 | 0.020 | 30.65 | 20 cells/delta |
| `k_r=2.3e-5 m/s` | stop 110 s | 1.709 | 0.020 | 3.168 | 12 cells/delta |
| `k_r=2.3e-4 m/s` | 6 h | 2.090 | 0.0457 | 1.484 | 6 cells/delta |
| `h0=0.5 mm` | 6 h | 2.649 | 0.247 | 1.132 | 4 cells/delta |
| `h0=2 mm` | 6 h | 2.633 | 0.246 | 1.064 | 6 cells/delta |

All completed baseline-family cases satisfy `Pe << Da`; the large-`Da`
concentration approximation is mathematically justified. The lowest-`k_r` case
still has `Da approximately 4.1e3`, but its `Da*` is only 0.008 and
`Da_top approximately 5.75`, which strongly reduces top growth relative to
internal residence and drives the chosen closure to its porosity limit.

## Answers to the ten research questions

1. **Reproduction:** yes for the reduced fixed-property equations under
   matching assumptions; this is verification, not validation.
2. **Location:** within roughly the upper 1.24 micrometres for 95% of baseline
   instantaneous precipitation.
3. **Motion:** the zone follows the moving interface and is nearly stationary
   in the distance coordinate after startup.
4. **Porosity reduction:** newly formed material can fall from 0.43 to about
   0.243 locally; the original column is affected only near its initial top.
5. **Velocity effect:** the fully dynamic model is about 3% faster than Model A,
   but B2 versus B1 differs by less than 0.1%.
6. **Duration:** the startup discrepancy is transient; the approximately 3%
   coefficient offset persists over the simulated 6 h.
7. **Self-similarity:** yes mathematically in the moving coordinate, within the
   tested interval and closures.
8. **Mechanism:** the baseline supports a modified coefficient, not a new
   physical growth mechanism. Low-`k_r` stops expose a missing hydraulic
   feedback rather than validate clogging.
9. **Continuum validity:** doubtful. Baseline `L95 = 1.24 um` and
   `delta_C = 0.511 um` are below the estimated `d_c = 3.5 um`.
10. **Needed evidence:** crystal-resolved morphology and size distributions,
    local interfacial area evolution, a defensible upscaled `beta`, independent
    kinetic measurements, permeability/capillary connectivity versus deposited
    mass, and spatially resolved concentration or precipitation data.
