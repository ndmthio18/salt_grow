# DuMux–Python 5×5 PNM comparison (`x_res = 0.08`)

This run compares the project's saturated-liquid Python PNM with the fully saturated limit of
the DuMux `2pncmin_isothermal` application on the exact same heterogeneous DGF.  It is not a
validation of drainage, moving menisci, or unsaturated drying.

## Case

- Original DGF: `input/N5x5_R270um_SD10um_L3.5mmx3.5mm_n1_top_varthroat.dgf`
- DGF SHA-256: `E48FF929AEAE9149DF9B575FD2BD612C7C19E111810D023DBF7013EE5341B825`
- Bottom row (`y = 0`, file pores 0, 5, 10, 15, 20): pressure and composition reservoir.
- Top row (`y = 3.5 mm`, file pores 25–29): boundary-layer evaporation.
- NaCl reservoir/initial mole fraction: `0.08` (`w = 0.22001805096761068`).
- Temperature: `293.15 K`; boundary layer: `1.6 mm`; far-field vapour mole fraction: `0.0047`.
- Gravity off; liquid saturation one; equilibrium precipitation; pore-volume and throat-radius
  feedback on; requested end time `40000 s`.

## Layout

- `input/`: immutable DGF and DuMux inputs.
- `python_dt{10,5,2p5}/`: self-driven Python runs; `2p5` is the preregistered fallback after
  the 10 s versus 5 s concentration check narrowly failed.
- `python_replay_dt5/`: Python transport driven by DuMux's five pore-resolved evaporation
  fluxes, separating internal-network error from evaporation-closure error.
- `dumux_dt{10,5}/raw/`: synchronized raw DuMux VTP/PVD/evaporation output.
- `dumux_diagnostics_1000s_1s/`: Float64 pore inventories/states and signed phase/component
  throat fluxes used for early pressure, flow and ledger checks.
- `normalized/`: common NumPy representation and DuMux flux-replay inputs.
- `tables/`: pore- and throat-resolved values behind the pressure figures.
- `figures/`: editable PDF/SVG and 600 dpi PNG pressure figures.
- `comparison_metrics.json`: preregistered numerical checks.
- `validation_criteria.json`: machine-readable pass/fail decisions; missing evidence is not
  counted as a pass.
- `REPORT_ZH.md`: evidence-bounded interpretation and GO/NO-GO decisions.

The standard DuMux VTP writer stores absolute liquid pressure as Float32, which removes the small
early network drop from its pressure baseline.  Validation commit `ab303806…` therefore adds
Float64 pore state/inventory and signed throat/component CSV diagnostics without changing the
equations.  The 940/1000 s pressure panels use those direct values.  The formal long runs predate
the diagnostic commit, so their late gauge pressure is reconstructed from
`volumeFluxW / transmissibilityW`; the reported residual is part of the late-field uncertainty.

The Python output also contains independent dissolved-plus-solid salt and water ledgers.  The
water ledger is intentionally reported even when it fails: equilibrium precipitation shrinks
the liquid pore volume before the displaced-volume source is applied on the next transport
step, so a non-closing result is evidence about the present coupling rather than a plotting
error.

