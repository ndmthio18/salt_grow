"""The small boxes of the test suite (they are cases like any other: every value here would
change with the experiment, so none of it belongs in a kernel or in a test helper)."""

# smallest 3-D case that still exercises every branch: 2 x 2 pores at lambda = 32, 3 voxels per
# pitch, a 0.5 mm box, the 2014 far field (18 mm, RH 7 %) and a 30 % brine
SMALL_3D = dict(lam=32, C0=0.30, Lt_res=22e-3, n_side=2, vpp=3, delta_eff_mm=18.0, rh=0.07,
                z_top_mm=0.5, tag="test", verbose=False)
# the 3-D reference case: c_k = 1e-6 seals inside 24 h, ell = d_sub, eps_min = 0.05
CASE_3D = dict(SMALL_3D, c_k=1e-6, t_end=24 * 3600, ell_over_dsub=1.0, eps_min=0.05)
# smallest 2-D case: two seeds at 34 um, 4 cells per pitch, a 0.4 mm box, 12 h
SMALL_2D = dict(d_sub=34e-6, n_seeds=2, cells_per_pitch=4, z_box_mm=0.4, t_end_h=12.0,
                snap_hours=0.0, verbose=False)
