# kernel_crossval.py output

## column_1x1x50_uniform  (N = 50 pores, demand uniform)
  kernel             dt [s]  steps  t_sat [h]     cum salt res / gain   cum water res / evap   max |step salt| [kg]  Picard
  network_step         60.0    116  1.919560    -1.32e-13            -3.27e-10            8.8e-22              -
  conservative_step    60.0    116  1.919636    +1.03e-11            +8.99e-12            4.2e-21              5
  network_step         30.0    231  1.917636    -4.29e-15            -1.98e-09            1.0e-21              -
  conservative_step    30.0    231  1.917673    +8.78e-12            +7.66e-12            4.4e-21              5
  network_step         15.0    461  1.916674    -1.03e-13            -3.54e-09            1.0e-21              -
  conservative_step    15.0    461  1.916693    -8.07e-13            -8.86e-13            1.0e-20              5
  network_step          7.5    920  1.916193    -2.40e-13            -3.52e-09            9.9e-22              -
  conservative_step     7.5    920  1.916202    -6.27e-11            -5.49e-11            1.5e-20              5
  dt [s]   |t_sat(net) - t_sat(cons)| [s]   max |w_surf(1 h) diff|
    60.0         0.2732                    2.784e-07
    30.0         0.1366                    1.393e-07
    15.0         0.0683                    6.965e-08
     7.5         0.0341                    3.483e-08
  Richardson t_sat: network_step order 1.00, limit 1.915713 h; conservative_step order 1.00, limit 1.915713 h; limits differ by 0.0000 s
  Stage-T reference (conservative_step, dt 30 s): 1.917673 h; here 1.917673 h (diff 1.71e-03 s)
  [13 s]

## hotspot_5x5x10  (N = 250 pores, demand hotspot)
  kernel             dt [s]  steps  t_sat [h]     cum salt res / gain   cum water res / evap   max |step salt| [kg]  Picard
  network_step         60.0     90  1.484066    +1.44e-12            -4.07e-10            4.7e-20              -
  conservative_step    60.0     90  1.484155    +7.58e-12            +6.07e-12            6.2e-20              5
  network_step         30.0    178  1.481829    +1.66e-11            -2.10e-09            1.1e-19              -
  conservative_step    30.0    178  1.481873    +6.80e-12            +5.47e-12            3.0e-20              5
  network_step         15.0    356  1.480711    +1.16e-11            -2.49e-09            3.6e-20              -
  conservative_step    15.0    356  1.480733    +6.36e-12            +5.17e-12            1.6e-20              5
  network_step          7.5    711  1.480152    +6.03e-12            -2.48e-09            1.0e-20              -
  conservative_step     7.5    711  1.480163    +6.31e-12            +4.87e-12            9.4e-21              5
  dt [s]   |t_sat(net) - t_sat(cons)| [s]   max |w_surf(1 h) diff|
    60.0         0.3185                    4.049e-07
    30.0         0.1585                    2.023e-07
    15.0         0.0793                    1.011e-07
     7.5         0.0396                    5.057e-08
  Richardson t_sat: network_step order 1.00, limit 1.479594 h; conservative_step order 1.00, limit 1.479594 h; limits differ by 0.0002 s
  [56 s]

Reading.  Both paths close the salt and water books to round-off with the same yardstick (the production
salt_step does not, Stage-T README section 0).  Their t_sat differ at first order in dt because the throat
densities are taken at t^n (network_step) or t^{n+1} (conservative_step); the Richardson limits agree to
the extrapolation error.  network_step stays the production kernel; conservative_step is kept for this
cross-check only.  (network_step does not expose its Picard count; conservative_step needs 5 at tol 1e-13.)
