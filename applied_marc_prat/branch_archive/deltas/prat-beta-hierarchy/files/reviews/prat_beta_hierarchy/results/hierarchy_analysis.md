# hierarchy_analysis.py output

## A. analytic maps
  d_c,crit (beta 2, Prat params): eps 0.43 -> 0.834 um; eps 0.20 -> 1.717 um
  Reaction-transport decides n* through d_c(n) < d_c,crit only; alpha and c_h never enter.

  closure layer n*, eps_s = 0.43 (rows r = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], cols d_sub[um] = [5.0, 10.0, 27.499999999999996, 56.0, 112.5, 250.0, 500.0, 1080.0])
  r=0.3:  1   2   3   3   4   5   5   6
  r=0.4:  1   2   3   4   5   6   6   7
  r=0.5:  1   2   4   5   6   7   8   9
  r=0.6:  1   3   4   6   7   9  10  12
  r=0.7:  1   3   6   8  10  12  14  16
  r=0.8:  1   4   8  12  15  18  21  25
  r=0.9:  1   6  16  23  29  37  44  51

  closure layer n*, eps_s = 0.2 (rows r = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], cols d_sub[um] = [5.0, 10.0, 27.499999999999996, 56.0, 112.5, 250.0, 500.0, 1080.0])
  r=0.3:  1   1   2   3   3   4   5   5
  r=0.4:  1   1   2   3   4   5   6   6
  r=0.5:  1   1   3   4   5   6   7   8
  r=0.6:  1   1   3   4   6   7   9  10
  r=0.7:  1   1   4   6   8  10  12  14
  r=0.8:  1   1   5   8  11  15  18  22
  r=0.9:  1   1   9  16  23  30  37  44

  trigger height h* = layer bottom of n* = 4 for 27.5 um, eps 0.43, r 0.5 [um] (rows alpha, cols c_h):
  alpha=1.25:    31.5     52.4    104.8    146.8    209.7    314.5
  alpha=1.50:    39.2     65.3    130.6    182.9    261.3    391.9
  alpha=2.00:    57.8     96.3    192.5    269.5    385.0    577.5
  alpha=3.00:   107.2    178.8    357.5    500.5    715.0   1072.5
  alpha=4.00:   173.2    288.8    577.5    808.5   1155.0   1732.5
  h* = c_h d_sub (alpha^(n*-1) - 1)/(alpha - 1): linear in c_h, geometric in alpha; the layer index n* is untouched.
  Exponent of h* ~ d_sub^m with n*(d_sub) from reaction-transport: m = 1 + ln(alpha)/ln(1/r) (plan/13), e.g.
   r = 0.4: alpha 1.25: m = 1.24, alpha 1.50: m = 1.44, alpha 2.00: m = 1.76, alpha 3.00: m = 2.20, alpha 4.00: m = 2.51
   r = 0.5: alpha 1.25: m = 1.32, alpha 1.50: m = 1.58, alpha 2.00: m = 2.00, alpha 3.00: m = 2.58, alpha 4.00: m = 3.00
   r = 0.6: alpha 1.25: m = 1.44, alpha 1.50: m = 1.79, alpha 2.00: m = 2.36, alpha 3.00: m = 3.15, alpha 4.00: m = 3.71
   r = 0.7: alpha 1.25: m = 1.63, alpha 1.50: m = 2.14, alpha 2.00: m = 2.94, alpha 3.00: m = 4.08, alpha 4.00: m = 4.89
  m is pure geometry (alpha, r); only the prefactor (through d_c,crit ~ D_eff(1-eps)/(beta^2 eps (1-C_eq) k_r)) is physics.

## B. numerical sweeps (27.5 um, eps_new 0.43, consistent unless stated)

  r sweep (alpha 2, c_h 1, up to 7 layers or 1.5 mm)
  r=0.4 hybrid                       closure layer 3 h    83.3 um t  0.86 h | analytic n* 3 | F plateaus 0.740 0.643 | transient [delta_C/v]: 2: 3.2
  r=0.4 consistent                   closure layer 4 h   192.7 um t  2.08 h | analytic n* 3 | F plateaus 0.755 0.674 0.600 | transient [delta_C/v]: 2: 1.1, 3: 1.6
  r=0.5 hybrid                       closure layer 4 h   193.0 um t  2.12 h | analytic n* 4 | F plateaus 0.740 0.668 0.587 | transient [delta_C/v]: 2: 3.4, 3: 1.4
  r=0.5 consistent                   closure layer 5 h   412.7 um t  4.59 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 1.2, 3: 3.0, 4: 2.0
  r=0.6 hybrid                       closure layer 4 h   196.2 um t  2.05 h | analytic n* 4 | F plateaus 0.740 0.688 0.631 0.570 | transient [delta_C/v]: 2: 3.8, 3: 1.5
  r=0.6 consistent                   closure layer 6 h   852.7 um t  9.53 h | analytic n* 4 | F plateaus 0.755 0.711 0.665 0.622 0.588 | transient [delta_C/v]: 2: 1.3, 3: 1.6, 4: 1.9, 5: 4.6
  r=0.7 hybrid                       closure layer 6 h   853.4 um t  9.58 h | analytic n* 6 | F plateaus 0.740 0.704 0.666 0.625 0.583 | transient [delta_C/v]: 2: 1.5, 3: 1.7, 4: 3.7, 5: 4.1
  r=0.7 consistent                   closure none | analytic n* 6 | F plateaus 0.755 0.724 0.693 0.661 0.631 0.604 | transient [delta_C/v]: 2: 1.5, 3: 3.4, 4: 3.9, 5: 2.2, 6: 2.6
  r=0.8 hybrid                       closure none | analytic n* 8 | F plateaus 0.740 0.718 0.695 0.671 0.646 0.620 | transient [delta_C/v]: 2: 1.7, 3: 3.7, 4: 0.0, 5: 0.0, 6: 4.7
  r=0.8 consistent                   closure none | analytic n* 8 | F plateaus 0.755 0.736 0.716 0.697 0.677 0.658 | transient [delta_C/v]: 2: 1.8, 3: 1.9, 4: 4.2, 5: 2.3, 6: 2.5

  alpha sweep (r 0.5, c_h 1, 5 layers)
  alpha=1.25 consistent              closure layer 5 h   158.8 um t  1.69 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 3.8, 4: 2.5
  alpha=1.5 consistent               closure layer 5 h   223.6 um t  2.43 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 3.8, 4: 5.0
  alpha=2.0 consistent               closure layer 5 h   412.7 um t  4.58 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 1.9, 4: 2.5
  alpha=3.0 consistent               closure layer 5 h  1100.2 um t 12.50 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 3.8, 4: 5.0
  alpha=4.0 consistent               closure layer 5 h  2337.7 um t 26.83 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 3.8, 4: 2.5

  c_h sweep (r 0.5, alpha 2, 5 layers)
  c_h=0.3 consistent                 closure layer 5 h   124.0 um t  1.38 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 1.5, 3: 3.8, 4: 2.5
  c_h=0.5 consistent                 closure layer 5 h   206.4 um t  2.29 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 1.5, 3: 1.9, 4: 5.0
  c_h=1.0 consistent                 closure layer 5 h   412.7 um t  4.58 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 1.9, 4: 2.5
  c_h=1.4 consistent                 closure layer 5 h   577.7 um t  6.42 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 3.0, 3: 1.9, 4: 5.0
  c_h=2.0 consistent                 closure layer 5 h   825.2 um t  9.17 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 1.5, 3: 1.9, 4: 2.5
  c_h=3.0 consistent                 closure layer 5 h  1237.7 um t 13.75 h | analytic n* 4 | F plateaus 0.755 0.694 0.634 0.586 | transient [delta_C/v]: 2: 1.5, 3: 3.8, 4: 5.0
  [B 238 s]
