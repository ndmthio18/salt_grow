# beta_analysis.py output

## A. beta from the growth cell
  Prat Appendix cell (A1, A2): A = 24 d^2, V = 4 d^3                     beta = 2.0000 (at d = 7.3: 2.0000)
  one free cube, six faces growing                                       beta = 2.0000 (at d = 7.3: 2.0000)
  cube attached by its base (five faces grow)                            beta = 1.6667 (at d = 7.3: 1.6667)
  cube with base + two lateral contacts (three faces grow)               beta = 1.0000 (at d = 7.3: 1.0000)
  cube with only its top face growing (layer-by-layer overgrowth)        beta = 0.3333 (at d = 7.3: 0.3333)
  hollow (hopper) cube, wall 0.2 d, inner faces wetted                   beta = 3.4694 (at d = 7.3: 3.4694)
  needle d x d x 3 d, h along the needle, only L grows                   beta = 14.0000 (= 4 L/d + 2)
  needle d x d x 10 d, h along the needle, only L grows                  beta = 42.0000 (= 4 L/d + 2)
  Dilation d -> lambda d leaves beta unchanged for every cell (A ~ d^2, V ~ d^3 => A/(dV/dd) ~ d^0).
  The gap e (porosity) does not enter A1-A4 at all: beta = 2 is independent of eps_s by construction.
  Replicating the cell laterally (N cells) multiplies A and dV/dd by N: beta unchanged.
  What changes beta: faces that cannot grow (contact, partial wetting) -> 2 - n_blocked/3; a habit whose
  growing area per volume differs from a cube (hollow/hopper cubes of Wijnhorst 2025 micro-CT, needles);
  face-dependent kinetics (effectively a different k_r per face); a layer-thickness convention other than
  h = d_t + e (A5).  None of these is a similarity transformation of the cubic cell.
  mass-balance view, eps = 0.43: a_v d_c = 2.85 growing-area per macro area in one crystal layer, i.e. beta_mb = a_v d_c/(1-eps) = 5.00 if the whole top crystal layer grew on all its wetted faces; Prat's beta = 2 instead counts cubes growing from small nuclei (left sketch of Figure 3).
  mass-balance view, eps = 0.2: a_v d_c = 4.00 growing-area per macro area in one crystal layer, i.e. beta_mb = a_v d_c/(1-eps) = 5.00 if the whole top crystal layer grew on all its wetted faces; Prat's beta = 2 instead counts cubes growing from small nuclei (left sketch of Figure 3).

## B. analytic sensitivity to beta (eqs 27, 29, 31)

  d_sub = 27.5 um, eps_s = 0.43, r = 0.5
  beta   Da*      F      v/U     C*_top    1-F    eps_final  d_c,crit[um]  n*  dlnv/dlnbeta
  0.50    0.507  0.4158  0.1447  2.634e-05  0.584  -0.371      13.351     1   0.584
  1.00    2.027  0.5874  0.2044  1.860e-05  0.413  +0.030       3.338     2   0.413
  1.41    4.054  0.6681  0.2325  1.496e-05  0.332  +0.147       1.669     3   0.332
  1.67    5.630  0.7035  0.2448  1.337e-05  0.296  +0.190       1.202     3   0.296
  2.00    8.107  0.7401  0.2575  1.172e-05  0.260  +0.230       0.834     4   0.260
  2.83   16.215  0.8011  0.2788  8.969e-06  0.199  +0.288       0.417     5   0.199
  3.00   18.242  0.8103  0.2820  8.553e-06  0.190  +0.297       0.371     5   0.190
  4.00   32.430  0.8506  0.2960  6.734e-06  0.149  +0.330       0.209     6   0.149
  5.00   50.672  0.8768  0.3051  5.553e-06  0.123  +0.350       0.134     6   0.123
  6.00   72.967  0.8952  0.3115  4.725e-06  0.105  +0.363       0.093     7   0.105
  8.00  129.719  0.9193  0.3199  3.639e-06  0.081  +0.380       0.052     8   0.081

  d_sub = 112.5 um, eps_s = 0.43, r = 0.5
  beta   Da*      F      v/U     C*_top    1-F    eps_final  d_c,crit[um]  n*  dlnv/dlnbeta
  0.50    2.073  0.5901  0.2054  3.738e-05  0.410  +0.034      13.351     2   0.410
  1.00    8.292  0.7422  0.2583  2.350e-05  0.258  +0.232       3.338     4   0.258
  1.41   16.583  0.8028  0.2794  1.798e-05  0.197  +0.290       1.669     5   0.197
  1.67   23.033  0.8276  0.2880  1.572e-05  0.172  +0.311       1.202     5   0.172
  2.00   33.167  0.8521  0.2965  1.349e-05  0.148  +0.331       0.834     6   0.148
  2.83   66.334  0.8906  0.3099  9.972e-06  0.109  +0.360       0.417     7   0.109
  3.00   74.626  0.8963  0.3119  9.461e-06  0.104  +0.364       0.371     7   0.104
  4.00  132.668  0.9201  0.3202  7.284e-06  0.080  +0.381       0.209     8   0.080
  5.00  207.293  0.9351  0.3254  5.922e-06  0.065  +0.390       0.134     8   0.065
  6.00  298.502  0.9453  0.3290  4.989e-06  0.055  +0.397       0.093     9   0.055
  8.00  530.671  0.9584  0.3335  3.794e-06  0.042  +0.405       0.052     10   0.042

  d_sub = 1080.0 um, eps_s = 0.43, r = 0.5
  beta   Da*      F      v/U     C*_top    1-F    eps_final  d_c,crit[um]  n*  dlnv/dlnbeta
  0.50   19.900  0.8169  0.2843  5.174e-05  0.183  +0.302      13.351     5   0.183
  1.00   79.601  0.8992  0.3129  2.848e-05  0.101  +0.366       3.338     7   0.101
  1.41  159.201  0.9266  0.3224  2.075e-05  0.073  +0.385       1.669     8   0.073
  1.67  221.113  0.9370  0.3261  1.780e-05  0.063  +0.392       1.202     8   0.063
  2.00  318.402  0.9469  0.3295  1.499e-05  0.053  +0.398       0.834     9   0.053
  2.83  636.805  0.9619  0.3347  1.077e-05  0.038  +0.407       0.417     10   0.038
  3.00  716.405  0.9640  0.3355  1.018e-05  0.036  +0.409       0.371     10   0.036
  4.00  1273.609  0.9727  0.3385  7.701e-06  0.027  +0.414       0.209     11   0.027
  5.00  1990.015  0.9781  0.3404  6.195e-06  0.022  +0.417       0.134     12   0.022
  6.00  2865.621  0.9817  0.3416  5.181e-06  0.018  +0.419       0.093     12   0.018
  8.00  5094.437  0.9862  0.3432  3.904e-06  0.014  +0.422       0.052     13   0.014
  Scalings: Da* ~ beta^2;  F = sqrt(Da*)/(1+sqrt(Da*));  dlnF/dlnbeta = 1 - F;  C*_top ~ 1/[beta(1-eps)(1-C_eq)k_r + const];
  d_c,crit ~ beta^-2, hence with r = 0.5 the closure layer shifts by ONE for every factor sqrt(2) in beta:
  n* = 2 + floor(ln(0.14 d_sub / d_c,crit) / ln 2) = 2 + floor(ln(0.14 d_sub beta^2 / K) / ln 2).

## C. numerical sweep (reference solver)

  homogeneous steady state, d_sub = 27.5 um, eps_new = 0.43
  beta  | eq29 v/U | hybrid v/U (rel)   C*_top   1-F_share  eps_final | consistent v/U   C*_top   1-F_share  eps_final
  1.00  | 0.2044   | 0.2044 (+1.3e-04)  1.860e-05  0.413  +0.0298   | 0.2209   2.010e-05  0.365  +0.1021
  1.41  | 0.2325   | 0.2325 (+1.1e-04)  1.496e-05  0.332  +0.1470   | 0.2418   1.556e-05  0.305  +0.1796
  2.00  | 0.2575   | 0.2576 (+8.6e-05)  1.172e-05  0.260  +0.2299   | 0.2627   1.195e-05  0.245  +0.2449
  2.83  | 0.2788   | 0.2788 (+6.6e-05)  8.969e-06  0.199  +0.2885   | 0.2816   9.060e-06  0.191  +0.2956
  4.00  | 0.2960   | 0.2960 (+5.0e-05)  6.735e-06  0.149  +0.3299   | 0.2975   6.769e-06  0.145  +0.3333
  6.00  | 0.3115   | 0.3115 (+3.5e-05)  4.725e-06  0.105  +0.3633   | 0.3122   4.736e-06  0.103  +0.3647

  homogeneous steady state, d_sub = 112.5 um, eps_new = 0.43
  beta  | eq29 v/U | hybrid v/U (rel)   C*_top   1-F_share  eps_final | consistent v/U   C*_top   1-F_share  eps_final
  1.00  | 0.2583   | 0.2583 (+1.1e-04)  2.351e-05  0.258  +0.2321   | 0.2634   2.397e-05  0.243  +0.2468
  1.41  | 0.2794   | 0.2794 (+8.7e-05)  1.798e-05  0.197  +0.2901   | 0.2822   1.816e-05  0.189  +0.2970
  2.00  | 0.2965   | 0.2965 (+6.6e-05)  1.349e-05  0.148  +0.3311   | 0.2980   1.356e-05  0.144  +0.3344
  2.83  | 0.3099   | 0.3100 (+4.9e-05)  9.972e-06  0.109  +0.3600   | 0.3107   9.997e-06  0.107  +0.3616
  4.00  | 0.3202   | 0.3202 (+3.6e-05)  7.285e-06  0.080  +0.3805   | 0.3206   7.294e-06  0.079  +0.3813
  6.00  | 0.3290   | 0.3290 (+2.5e-05)  4.989e-06  0.055  +0.3970   | 0.3291   4.992e-06  0.054  +0.3974

  layered 27.5 um (r 0.5, alpha 2, c_h 1, eps_new 0.43, 5 layers): closure layer / h / t  vs analytic n*
  beta  | analytic n* | hybrid                     | consistent
  1.00  |      2      | layer 2, h   28.6 um, t  0.34 h | layer 3, h   82.9 um, t  0.95 h
  1.41  |      3      | layer 3, h   83.3 um, t  0.94 h | layer 4, h  192.8 um, t  2.17 h
  2.00  |      4      | layer 4, h  193.1 um, t  2.12 h | layer 5, h  412.7 um, t  4.58 h
  2.83  |      5      | layer 5, h  412.9 um, t  4.47 h | none in 5 layers
  4.00  |      6      | none in 5 layers             | none in 5 layers
  6.00  |      7      | none in 5 layers             | none in 5 layers
  [C 150 s]
