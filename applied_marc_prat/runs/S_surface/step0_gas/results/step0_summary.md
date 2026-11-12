# run.py output (Step 0, gas-only single-pore evaporation)

## settings
  T          295.15
  RH         0.07
  C_brine    0.264
  pore_side  0.00054
  L_box      0.0007
  delta_eff  0.0206
  a_voxel    5e-05
  h_box      {'A': 0.002, 'B': 0.004}
  D_v        2.4486936347474e-05
  rho_g      1.1752512392969043
  W_sat      0.016621258564030007
  rho_D      2.877830228895323e-05
  a_w        0.7530474533159688
  X_wet      0.012516596432549033
  X_top      0.0011634880994821006

## the four calculations (+ supplementary shifted footprint)
  case  wet       h_box  L_ext   E [kg/s]          j_proj [kg/m2/s]  j_wet [kg/m2/s]   A_wet [m2]  area err   imbalance  residual   L_spread [um]
  A     plane      2.0    18.6   7.7715708703e-12  1.5860348715e-05  1.5860348715e-05  4.900e-07   +0.0000  +6.9e-13  5.9e-15     -0.0000
  A     pore       2.0    18.6   7.7437260730e-12  1.5803522598e-05  3.0974904292e-05  2.500e-07   -0.1427  +9.5e-13  7.9e-15     74.0732
  A     pore_shifted  2.0    18.6   7.7553998436e-12  1.5827346620e-05  2.5637685433e-05  3.025e-07   +0.0374  +8.0e-13  7.4e-15     42.9537
  B     plane      4.0    16.6   7.7715708703e-12  1.5860348715e-05  1.5860348715e-05  4.900e-07   +0.0000  +1.2e-12  8.4e-15     -0.0000
  B     pore       4.0    16.6   7.7437260730e-12  1.5803522598e-05  3.0974904292e-05  2.500e-07   -0.1427  +9.6e-13  1.1e-14     74.0732
  B     pore_shifted  4.0    16.6   7.7553998436e-12  1.5827346620e-05  2.5637685433e-05  3.025e-07   +0.0374  +1.0e-12  1.0e-14     42.9537

## V0.1 planar analytical reference
  A: j_1d 1.5860348715e-05  j_num 1.5860348715e-05  rel err +6.61e-13
  B: j_1d 1.5860348715e-05  j_num 1.5860348715e-05  rel err +1.11e-12

## V0.2 vapour mass conservation (E - E_top) / E_top and discrete residual ||Ax-b||/||b||
  plane A: imbalance +6.87e-13  residual 5.95e-15
  plane B: imbalance +1.23e-12  residual 8.36e-15
  pore  A: imbalance +9.47e-13  residual 7.93e-15
  pore  B: imbalance +9.62e-13  residual 1.13e-14

## V0.3 domain-height independence of the single-pore result (B = 4 mm as reference)
  E            A 7.7437260730e-12  B 7.7437260730e-12  (A - B)/B +7.59e-14
  j_projected  A 1.5803522598e-05  B 1.5803522598e-05  (A - B)/B +7.59e-14
  j_wet        A 3.0974904292e-05  B 3.0974904292e-05  (A - B)/B +7.59e-14
  L_spread     A 7.4073232963e-05  B 7.4073232965e-05  (A - B)/B -2.12e-11
  lateral non-uniformity of X at the box top, A: 4.84e-10  (first lateral mode exp(-2 pi h / L_box) = 1.60e-08)
  lateral non-uniformity of X at the box top, B: 2.44e-15  (first lateral mode exp(-2 pi h / L_box) = 2.55e-16)
  E_pore / E_plane = 0.996417  (delta_eff / (delta_eff + L_spread) = 0.996417); per wet area the opening evaporates 1.9530 x the plane flux

## V0.4 pore-opening rasterisation (production footprint rule of coupling.interface)
  physical side 540 um, area 2.9160e-07 m2; voxel 50 um -> 10 x 10 cells = 500 um, area 2.5000e-07 m2, relative error -0.1427
  supplementary: pore centre shifted by half a voxel -> 11 x 11 cells, area error +0.0374, E changes by +1.51e-03 (the far-field resistance delta_eff dominates the spreading resistance L_spread, so E is insensitive to the footprint area)
  edge enhancement of the base flux on the wet cells (max / mean): 1.724

[1.8 s]
