"""The production cases of the 3-D closed loop (``saltv2.coupling.driver3d.run``).

Regime: evaporation-wicking, the conditions of Veran-Tissoires & Prat (2014) section 6.1.
"""
# ---- the production cases (2026-09-11): the 2014 section-6.1 conditions (RH 7 %, C0 = 25 %,
# L_s = 22 mm, delta = 15 mm recess -> delta_eff = 18 mm reproduces Table 1 within +-10 %),
# c_k = 1e-5 (the 2-D recalibration after front tracking, runs/E2_feedback/README.md section
# 12), 5x5 pores.  vpp 6 is the 2-D mesh at which front tracking made the result
# grid-independent (cpp 6 reproduced cpp 36 to 2-4 %); the 3-D mesh control is cases 7/8.
# Box heights sit under 0.8 * delta_eff = 14.4 mm; the coarse lam use the 2-D sealing heights
# (0.36-0.65 mm at 34-68 um, a few mm at 270 um) plus margin.
CRUST_COMMON = dict(C0=0.25, Lt_res=22e-3, n_side=5, delta_eff_mm=18.0, rh=0.07, tag="crust",
                    c_k=1e-5, eps_min=0.02, ell_over_dsub=0.25, vpp=6, t_end=150.0 * 3600)
SCAN_CRUST = [dict(lam=32, z_top_mm=3.0, **CRUST_COMMON),        # 0: 34 um  ~ 5-50 um class
              dict(lam=16, z_top_mm=4.0, **CRUST_COMMON),        # 1: 68 um  ~ 50-62 um class
              dict(lam=8, z_top_mm=6.0, **CRUST_COMMON),         # 2: 135 um ~ 100-125 um class
              dict(lam=1, z_top_mm=8.0, **CRUST_COMMON),         # 3: 1 mm
              dict(lam=16, z_top_mm=4.0, **{**CRUST_COMMON, "vpp": 12}),   # 4: mesh control for 1
              # 5, 6: seed spacing at 1 mm beads (added 2026-09-12 after case 3 merged into a slab
              # by 12 h).  pitch_factor widens the pore pitch at the same pore radius, vpp is raised
              # with it so the voxel stays 116.7 um: spacing 2.1 mm (90x90 lateral) and 3.5 mm
              # (150x150) against 0.7 mm in case 3 and the 5.3 mm Veran-Tissoires & Prat 2014
              # figure 20 gives at delta = 15 mm.  Tests whether patchy needs gap >> height.
              dict(lam=1, z_top_mm=8.0, pitch_factor=3.0, **{**CRUST_COMMON, "vpp": 18}),
              dict(lam=1, z_top_mm=8.0, pitch_factor=5.0, **{**CRUST_COMMON, "vpp": 30}),
              # 7, 8: mesh control in a box the event-driven step can afford.  Case 4 (vpp 12 on
              # the 5x5 box) advanced 0.7 h of model time in 4.8 h of wall time near the sealing
              # event -- 3600 columns sealing one voxel per step -- and was stopped.  The same
              # comparison on the 2x2 box (both meshes see the same elevated-plane artefact,
              # README section 9.10, so only their mutual agreement is read, never E itself).
              dict(lam=16, z_top_mm=1.5, **{**CRUST_COMMON, "n_side": 2, "vpp": 6, "t_end": 60.0 * 3600}),
              dict(lam=16, z_top_mm=1.5, **{**CRUST_COMMON, "n_side": 2, "vpp": 12, "t_end": 60.0 * 3600})]
TABLES = {"scancrust": (SCAN_CRUST, "vapor_crust_")}
