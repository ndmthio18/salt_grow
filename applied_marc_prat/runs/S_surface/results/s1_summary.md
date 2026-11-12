# s1_run.py output

## 1. vapour benchmark: full-plane wet bottom at a_w(C_sat) against the 1-D flux
  j_1d 1.702089e-05  j_num 1.702089e-05  rel err 8.8e-12  imbalance 9.3e-12  (2014 mean flux J_MEAN 1.5854e-05: the legacy air constants; here moist_air at 21 C)

## 2. first saturation without growth: prototype depth vs the 50-layer reference (G2 head 22 mm)
  n_layers  demand      t_sat [h]   vs D-50     E0 [kg/s]    E(t_sat) [kg/s]  p_req G2 / G1 [Pa]  P_open  res salt / mass / water [kg]
    50      prescribed  1.917636   +0.0000%   7.7685e-12   7.7685e-12       257.0 /  408.6     556   -6.9e-23 / -1.1e-16 / -1.1e-16
    10      prescribed  1.895034   -1.1786%   7.7685e-12   7.7685e-12       257.4 /  409.3     556   +4.0e-21 / -4.1e-17 / -4.1e-17
    15      prescribed  1.917420   -0.0113%   7.7685e-12   7.7685e-12       257.1 /  408.8     556   +3.6e-21 / -1.1e-16 / -1.1e-16
    25      prescribed  1.917636   -0.0000%   7.7685e-12   7.7685e-12       257.0 /  408.6     556   -5.3e-23 / -1.1e-16 / -1.1e-16
    D-10 vs D-50 (prescribed): -1.1786%
    D-15 vs D-50 (prescribed): -0.0113%
    D-25 vs D-50 (prescribed): -0.0000%
    50      feedback   1.672215   +0.0000%   8.5615e-12   8.3065e-12       257.0 /  408.6     556   -2.4e-21 / -6.6e-17 / -6.6e-17
    10      feedback   1.660430   -0.7047%   8.5615e-12   8.3062e-12       257.3 /  409.2     556   +3.3e-21 / -4.6e-17 / -4.6e-17
    15      feedback   1.672150   -0.0039%   8.5615e-12   8.3065e-12       257.0 /  408.7     556   +9.4e-22 / -6.6e-17 / -6.6e-17
    25      feedback   1.672215   -0.0000%   8.5615e-12   8.3065e-12       257.0 /  408.6     556   -3.2e-22 / -6.6e-17 / -6.6e-17
    D-10 vs D-50 (feedback): -0.7047%
    D-15 vs D-50 (feedback): -0.0039%
    D-25 vs D-50 (feedback): -0.0000%
  Reading: the prescribed D-50 value is the kernel cross-validation anchor (1.917636 h with network_step).  The
  lumped prototype keeps the series resistance of the unresolved throats but drops their pore storage, so it
  enriches EARLIER (less downward capacity); the Stage-T shallow 10-layer column (reservoir at 7 mm, a perfect
  sink) was later by +0.25 %.  Convergence with the resolved depth is tabulated above.

## 3. growth: single-scale Prat structure after nucleation at C_sat
  case                  t_nuc [h]  m_seed [kg]  h(end) [um]   v_end [nm/s]  v/v_prat  X/E (end)  conv C_eq/(1-C_eq)  E(end)/E0  res salt / mass / water  wall [s]
  growth_D50_C025         1.6750   1.435e-11      531.56      7.6512    0.81328   0.32216   0.35870            1.0079   -1.1e-22 / -6.6e-17 / -6.6e-17  638
  growth_D25_C025         1.6750   1.435e-11      537.20      7.7810    0.82708   0.32763   0.35870            1.0079   -3.2e-22 / -6.6e-17 / -6.6e-17  619
  growth_D15_C025         1.6750   1.444e-11      562.44      7.9509    0.84272   0.33382   0.35870            1.0108   +3.4e-21 / -6.6e-17 / -6.6e-17  698
  growth_D10_C025         1.6667   1.951e-11      587.41      7.9190    0.83934   0.33248   0.35870            1.0108   +4.5e-21 / -4.6e-17 / -4.6e-17  592
  growth_D50_Csat         0.0083   1.706e-10      198.22      9.2233    1.00000   0.39613   0.35870            1.0111   -4.8e-21 / +2.8e-20 / +3.3e-20  157
  growth_D10_Csat         0.0083   1.706e-10      198.22      9.2233    1.00000   0.39613   0.35870            1.0111   +3.4e-22 / +3.4e-20 / +3.4e-20  182
  Reading: with the reservoir at C_sat the column delivers exactly C_eq/(1-C_eq) per unit of evaporated water and
  dh/dt equals prat2026.core.growth_rate (eq 29) to round-off; with C0 = 0.25 the column delivers less (X/E below
  the conversion) and the structure grows slower than eq 29 by that ratio -- a supply statement, not a kinetic one.
  E(end)/E0 > 1 is the far-field approach of a growing wet structure (plan/19 section 3.4), not screening.
