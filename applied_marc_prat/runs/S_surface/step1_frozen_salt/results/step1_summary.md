# run.py output (Step 1, prescribed wet salt pillar, gas only)

## 0. face-selective assembly against the production kernel (A2 = SALT_WET pillar; h = 0 = bottom opening)
  a = 50 um (direct): A2 E rel -2.3e-13, X max rel 5.5e-16, E_top rel -1.2e-13, kernel base flux 0.0e+00; h=0 E rel -7.8e-14, X max rel 6.9e-16  [0.7 s]
  a = 20 um (cg): A2 E rel -3.5e-13, X max rel 3.1e-14, E_top rel +5.8e-13, kernel base flux 0.0e+00; h=0 E rel +0.0e+00, X max rel 0.0e+00  [1.9 s]

## 1. primary matrix: 20 um voxels, 2 mm box (far field 20.6 mm)
  A0             A0  h 0.00 mm  a 20 um  box 2.0 mm  E_tot 7.7570212943e-12  top 7.757021e-12  side 0.000000e+00  far 7.7570212940e-12  imb +3.8e-11  res 8.8e-13  n 122500     0.5 s
  A1_h0.0        A1  h 0.00 mm  a 20 um  box 2.0 mm  E_tot 7.7570212943e-12  top 7.757021e-12  side 0.000000e+00  far 7.7570212940e-12  imb +3.8e-11  res 8.8e-13  n 122500     0.6 s
  A1_h0.1        A1  h 0.10 mm  a 20 um  box 2.0 mm  E_tot 7.7955373789e-12  top 7.795537e-12  side 0.000000e+00  far 7.7955373786e-12  imb +3.7e-11  res 9.0e-13  n 118855     0.9 s
  A1_h0.3        A1  h 0.30 mm  a 20 um  box 2.0 mm  E_tot 7.8722240774e-12  top 7.872224e-12  side 0.000000e+00  far 7.8722240772e-12  imb +3.5e-11  res 9.9e-13  n 111565     1.1 s
  A1_h0.5        A1  h 0.50 mm  a 20 um  box 2.0 mm  E_tot 7.9504126294e-12  top 7.950413e-12  side 0.000000e+00  far 7.9504126292e-12  imb +3.2e-11  res 9.2e-13  n 104275     0.8 s
  A2_h0.0        A2  h 0.00 mm  a 20 um  box 2.0 mm  E_tot 7.7570212943e-12  top 7.757021e-12  side 0.000000e+00  far 7.7570212940e-12  imb +3.8e-11  res 8.8e-13  n 122500     1.4 s
  A2_h0.1        A2  h 0.10 mm  a 20 um  box 2.0 mm  E_tot 7.8041821779e-12  top 5.729256e-12  side 2.074927e-12  far 7.8041821776e-12  imb +3.6e-11  res 8.6e-13  n 118855     0.8 s
  A2_h0.3        A2  h 0.30 mm  a 20 um  box 2.0 mm  E_tot 7.8812209805e-12  top 5.750506e-12  side 2.130715e-12  far 7.8812209803e-12  imb +3.2e-11  res 9.9e-13  n 111565     0.5 s
  A2_h0.5        A2  h 0.50 mm  a 20 um  box 2.0 mm  E_tot 7.9595892017e-12  top 5.807632e-12  side 2.151958e-12  far 7.9595892014e-12  imb +3.1e-11  res 9.0e-13  n 104275     0.6 s

  ratios (E0 = A0 bare opening)
  h [mm]   R_height = E_A1/E0   R_sides = E_A2/E_A1   R_total = E_A2/E0   1-D control d/(d-h)   A2 side share   j_top_mean A1 / A2   j_side_mean A2
  0.00     1.000000             1.000000             1.000000         1.000000            0.0000        2.6602e-05 / 2.6602e-05   0.0000e+00
  0.10     1.004965             1.001109             1.006080         1.004878            0.2659        2.6734e-05 / 1.9648e-05   9.6061e-06
  0.30     1.014851             1.001143             1.016011         1.014778            0.2704        2.6997e-05 / 1.9721e-05   3.2881e-06
  0.50     1.024931             1.001154             1.026114         1.024876            0.2704        2.7265e-05 / 1.9916e-05   1.9926e-06

## 2. grid comparison: 50 um (the Step 0 grid, 10 x 10 cells = 500 um) at A0, A1 and A2 h = 0.3 mm
  A0_50um        A0  h 0.00 mm  a 50 um  box 2.0 mm  E_tot 7.7437260730e-12  top 7.743726e-12  side 0.000000e+00  far 7.7437260730e-12  imb +8.7e-13  res 7.9e-15  n   7840     0.3 s
      vs 20 um: E_total -1.714e-03, j_top_mean +1.644e-01; A_open_raster 2.500e-07 vs 2.916e-07 (phys 2.916e-07)
  A1_h0.3_50um   A1  h 0.30 mm  a 50 um  box 2.0 mm  E_tot 7.8583944200e-12  top 7.858394e-12  side 0.000000e+00  far 7.8583944199e-12  imb +5.7e-13  res 7.7e-15  n   7240     0.2 s
      vs 20 um: E_total -1.757e-03, j_top_mean +1.644e-01; A_open_raster 2.500e-07 vs 2.916e-07 (phys 2.916e-07)
  A2_h0.3_50um   A2  h 0.30 mm  a 50 um  box 2.0 mm  E_tot 7.8770055014e-12  top 5.259655e-12  side 2.617350e-12  far 7.8770055014e-12  imb +8.5e-13  res 4.3e-15  n   7240     0.2 s
      vs 20 um: E_total -5.349e-04, j_top_mean +6.684e-02, j_side_mean +3.267e-01; A_open_raster 2.500e-07 vs 2.916e-07 (phys 2.916e-07)
  Step 0 regression: A0 50 um E = 7.7437260730e-12 vs Step 0 7.7437260730e-12 (-7.8e-14)

## 3. domain truncation: A2 h = 0.5 mm in a 4 mm box (L_ext 16.6 mm) against the 2 mm box
  A2_h0.5_box4   A2  h 0.50 mm  a 20 um  box 4.0 mm  E_tot 7.9595892019e-12  top 5.807632e-12  side 2.151958e-12  far 7.9595892014e-12  imb +6.4e-11  res 9.3e-13  n 226775     1.7 s
  (4 mm - 2 mm) / 2 mm: E_total +2.88e-11, E_top +1.98e-11, E_side +5.31e-11, j_top_mean +1.98e-11, j_side_mean +5.31e-11
  local flux, max |4 mm - 2 mm| / max: top face 4.17e-10, side faces 6.81e-10

[15 s total]
