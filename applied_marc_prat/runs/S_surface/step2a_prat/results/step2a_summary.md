# run.py output (Step 2A, Prat single-scale growth law)

## A. analytic Prat law (prat2026/core.py)
  a_v = 36 (1-eps)/d_sub (core, eq 17 rounded)
    a_v 1.824000e+05 1/m  D_eff 8.524670e-10 m2/s  delta_C 1.08958 um  Da* 32.903681  F 0.8515478
    C*_top 1.553249e-06  v 3.95112268e-09 m/s = 14.22404 um/h  v/U_evap 0.2963342  U/U_evap 1.432870
    salt delivered 5.739130e-06 kg/m2/s -> front 4.887144e-06 (F), internal 8.519868e-07 (1-F); eps_final closure 0.330631
  a_v = 5 (1-eps)/d_c, d_c = 0.14 d_sub (unrounded)
    a_v 1.809524e+05 1/m  D_eff 8.524670e-10 m2/s  delta_C 1.09393 um  Da* 33.166910  F 0.8520507
    C*_top 1.554167e-06  v 3.95345629e-09 m/s = 14.23244 um/h  v/U_evap 0.2965092  U/U_evap 1.432870
    salt delivered 5.739130e-06 kg/m2/s -> front 4.890030e-06 (F), internal 8.491004e-07 (1-F); eps_final closure 0.331026
  36 vs 5/0.14: a_v +8.000e-03, Da* -7.937e-03, F -5.903e-04, v -5.903e-04
  large-Da check at h0 = 0.3 mm: Da = 7.581e+04, Pe = 1.091e-02, tanh(sqrt Da) = 1.0000000000000000, eq 28 vs eq 29: +2.1e-16
  reference trajectory: h(t) = 0.3 mm + 14.23244 um/h x t (unrounded a_v; the reference solver uses this a_v)

## B. independent reference (V0), steady travelling wave, reaction-zone refinement
  n_band  N    v/U_evap     vs eq 29     F_prat (v/v_max)  F_num (P_i/(P_i+S))  C*_top       eps_final   closure 1-(1-eps)/F_num  res salt   res solid  [s]
    20   123  0.2965340  +8.365e-05   0.8521220         0.8521220            1.55430e-06  0.331082    0.331082               +1.3e-10  +4.7e-11  1.7
    40   146  0.2965201  +3.658e-05   0.8520819         0.8520819            1.55422e-06  0.331050    0.331050               +4.4e-10  +1.9e-10  0.4
    80   169  0.2965156  +2.146e-05   0.8520690         0.8520690            1.55420e-06  0.331040    0.331040               -5.1e-10  +3.2e-11  0.4
  profile: max |C*_num - C*_top exp(-xi/delta_C)| / C*_top = 2.16e-04 / 3.92e-05 / 1.40e-05 (n_band 20/40/80)
  porosity range in the window (n_band 40): [0.331050, 0.427557]
  Richardson (n_band 20/40/80, v): order 1.64

## B'. transient reference from h0, 24 h, against the linear trajectory h0 + v_Prat t
  dt 138.4 s, 625 steps, 0.4 s; v settles to 1e-6 of the steady value at t = 4427 s (16.0 delta_C/v)
  h(24 h) numeric 641.5748 um, linear h0 + v_Prat t 641.5786 um: growth error -1.13e-05; v(24 h) vs eq 29 +3.66e-05; offset of h from the steady line h0 + v_ss t: -16.362 nm
  cumulative 24 h: salt in 4.931499e-01 kg/m2, precipitated (front + internal) 4.931499e-01, acc 8.771e-10, res salt +3.86e-10, res solid +1.18e-07

## C. salt, solid and liquid-volume balances of the steady V0 solution (n_band 40)
  salt [kg/m2/s]: supplied rho W_L C_eq 5.739130e-06; front P_i 4.890209e-06; internal S 8.489215e-07; dissolved accumulation 0 (steady); residual +4.70e-10
      eq 29 delivery C_eq j/(1-C_eq) = 5.739130e-06 (supplied/delivered - 1 = -2.22e-10)
      front share P_i/(P_i+S) = 0.8520819 vs F = 0.8520507 (+3.66e-05); internal share 0.1479181 vs 1-F 0.1479493
  solid: front rho_cr (1-eps_new) v = 4.890209e-06 = P_i (+1.7e-16); internal volume S/rho_cr 3.912080e-10 m/s = pore volume lost (eps_new - eps_final) v 3.912080e-10 (-9.2e-10)
  liquid: W_L - W_i = 7.074346e-10 m/s = brine volume consumed S/rho 7.074346e-10 (-1.0e-09);
      W_i = U_evap + (rho_cr/rho)(1-eps_new) v = 1.3333e-08 + 4.0752e-09 = 1.740851e-08; water reaching the interface rho W_i - P_i = 1.600000e-05 = j (+0.0e+00)
      new pore space eps_new v = 1.7000e-09 m/s, filled with brine 2.0401e-06 kg/m2/s (part of the liquid advected with the solid, not of W)
  porosity: eps_final 0.331050; 1-(1-eps_new)/F_num 0.331050 (+2.8e-10); 1-(1-eps_new)/F_analytic 0.331026 (+7.4e-05)
  constant-density closure: brine volume removed per precipitated salt S/rho; the water share (1-C_eq) S = 6.248e-07 kg/m2/s = 3.91 % of j is liquid that physically stays (as diluted brine): an idealisation of eq 4, not a conservation error of the solver

[3 s]
