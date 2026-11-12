# solver_audit.py output

## A. layer 2 -> 3 transition and closure under discretisation changes (27.5 um, eps_new 0.43, consistent)
  setting                 N    F2        F3        tau_23 [dC/v]  tau_34   closure layer  h [um]    t [h]    max Picard  cum res salt/solid
  n_band=20               169  0.69461  0.63482     2.67          3.50         5          412.71   4.577    18        1.7e-11 / -5.4e-12
  n_band=40 (reference)   192  0.69438  0.63425     2.67          3.49         5          412.69   4.584    29        2.1e-11 / -5.8e-12
  n_band=80               216  0.69427  0.63396     2.67          2.99         5          412.67   4.587    30        8.0e-12 / -1.3e-13
  n_band=160              239  0.69421  0.63381     2.67          3.48         5          412.70   4.589    30        1.3e-10 / -1.5e-11
  dt x0.2                 192  0.69438  0.63425     1.75          2.09         5          412.67   4.584    30        3.9e-11 / -9.6e-12
  dt x4                   192  0.69438  0.63425     4.58           nan         4          193.16   2.000    30        3.5e-11 / -1.6e-11
  stretch 1.015           335  0.69421  0.63393     2.67          2.99         5          412.70   4.587    30        -4.6e-12 / -8.1e-12
  L = 120 delta_C         216  0.69438  0.63425     2.67          3.49         5          412.69   4.584    30        1.5e-11 / -2.2e-14
  harmonic face average   192  0.69439  0.63429     2.67          3.49         5          412.69   4.583    30        8.8e-12 / -5.5e-12
  hybrid variant          192  0.66819  0.58745     3.53           nan         4          193.06   2.123    30        1.9e-10 / 2.7e-07
  dt x4, converged check  192  0.69438  0.63425     4.58          5.98         5          413.46   4.593    30        5.1e-11 / -2.3e-14
  dt x8, converged check  192  0.69438  0.63425     9.15         11.96         5          414.28   4.602    30        9.5e-13 / -2.5e-11
  dt x8, iterate check    192  0.69438  0.63425     9.15           nan         4          194.02   2.009    28        2.0e-13 / -6.6e-13
  Reading: the plateau F of a layer is grid-, dt- and averaging-independent to the digits shown and the transient after
  a boundary lasts a few delta_C/v (seconds to minutes) against layer transit times of hours.  DEFECT of the original
  solver: closure is tested on every Picard iterate, and the first iterate uses the lagged sink, so a step of
  >= 2 delta_C/v overshoots and closes one layer EARLY (layer 4 instead of 5).  Testing the converged step (extended copy,
  closure_check='converged') makes the closure layer dt-independent up to 4 delta_C/v.  All results of the previous
  review used dt = 0.5 delta_C/v, where both checks agree; with the converged check the closure height moves by
  < one band thickness per factor 4 in dt.
  [A 78 s]

## B. width of the recorded eps_final(z) jump at the layer 1 -> 2 boundary (10 %-90 % of the jump)
  n_band=20              jump 0.2448 -> 0.1801, 10-90 % width  10.95 um  (delta_C,1 = 0.541 um, window L = 32.5 um)
  n_band=40 (reference)  jump 0.2447 -> 0.1798, 10-90 % width  10.84 um  (delta_C,1 = 0.541 um, window L = 32.5 um)
  n_band=80              jump 0.2447 -> 0.1797, 10-90 % width  10.84 um  (delta_C,1 = 0.541 um, window L = 32.5 um)
  n_band=160             jump 0.2447 -> 0.1796, 10-90 % width  10.84 um  (delta_C,1 = 0.541 um, window L = 32.5 um)
  stretch 1.015          jump 0.2446 -> 0.1790, 10-90 % width   8.46 um  (delta_C,1 = 0.541 um, window L = 32.5 um)
  L = 120 delta_C        jump 0.2462 -> 0.1855, 10-90 % width  16.70 um  (delta_C,1 = 0.541 um, window L = 64.9 um)
  Reading: the physical transition is ~delta_C (the band); the recorded width is set by the coarse cells near the
  window bottom (implicit upwind advection of eps over L): it shrinks with the stretching, not with n_band, and it
  does not affect per-layer medians or closure.  Figures s2_* of the previous review show this smeared width.

## C. dense birth-porosity scan for the eps_final interpretation (target 0.43), both variants
  d_c = 0.481 um hybrid                : max eps_final 0.9544 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.883, 0.999]; closed (eps -> 0) at 2 of 25 points
  d_c = 0.481 um consistent n_band 40  : max eps_final 0.9990 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.947, 0.999]; closed (eps -> 0) at 0 of 25 points
  d_c = 0.481 um consistent n_band 80  : max eps_final 0.7341 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.947, 0.999]; closed (eps -> 0) at 0 of 25 points
  d_c = 0.241 um hybrid                : max eps_final 0.9360 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.931, 0.999]; closed (eps -> 0) at 4 of 25 points
  d_c = 0.241 um consistent n_band 40  : max eps_final 0.9978 at eps_new 0.9978; eps_final >= 0.43 for eps_new in [0.998, 0.998]; closed (eps -> 0) at 5 of 25 points
  d_c = 0.241 um consistent n_band 80  : max eps_final 0.9990 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.994, 0.999]; closed (eps -> 0) at 2 of 25 points
  d_c = 0.120 um hybrid                : max eps_final 0.9097 at eps_new 0.9990; eps_final >= 0.43 for eps_new in [0.969, 0.999]; closed (eps -> 0) at 6 of 25 points
  d_c = 0.120 um consistent n_band 40  : max eps_final 0.0683 at eps_new 0.9917; NO eps_new in [0.43, 0.999] reaches 0.43; closed (eps -> 0) at 15 of 25 points
  d_c = 0.120 um consistent n_band 80  : max eps_final 0.9987 at eps_new 0.9987; eps_final >= 0.43 for eps_new in [0.999, 0.999]; closed (eps -> 0) at 10 of 25 points
  Qualification: the scan covers 25 points in 1 - eps_new from 0.57 to 1e-3 on two grids.  Up to eps_new = 0.99 the
  consistent model reaches eps_final = 0.43 on neither grid for d_c <= 0.24 um, while the linearised model does
  (eps_new ~ 0.93).  Beyond 0.99 a grid-DEPENDENT sliver appears (eps_new >= 0.994-0.998) where a_v = 5(1-eps)/d_c -> 0
  makes the band inert and eps_final ~ eps_new: Prat's packing closure has no meaning at 99 % void, and the numerics
  are not converged there.  The previous review's 'infeasible' claim is therefore restated as: no birth porosity
  <= 0.99 reaches the target; it rests on the local-eps closure and is a finite-scan statement, not a theorem.
  [C 117 s]
