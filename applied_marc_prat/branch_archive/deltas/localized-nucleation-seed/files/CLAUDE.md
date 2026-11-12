# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

`salt_grow_v2` is the workspace for a **new salt-growth model with the closed feedback loop**

    evaporation → capillary liquid supply → ion transport → crystallization
    → structure (geometry / porosity / permeability) change → evaporation

that is meant to be coupled to a pore-network model (PNM) later. It is the successor of
`../salt_grow` (git, has the Python venv, the 2014 lattice model, `new_model/efflorescence3d`,
and the validation datasets). This directory is a git repository (branch `master`) and has no venv
of its own. Current contents:

| Path | Role |
|---|---|
| `ideas/` | Chinese design notes that *are* the model plan. `v2_建模方案.md` is the execution plan (goal, numerical design, `saltv2/` package layout, milestones, calibration plan, risks, fixed decisions) and `plan/00_总览.md`…`08_S8_PNM接口.md` split it into eight sub-plans S1–S8, each with physical basis, implementation, acceptance tests, and what may / may not be claimed afterwards — implement in that order, one sub-plan at a time; `patchy_crusty_transition_model.md` (user's, 2026-09-02) makes the patchy–crusty transition itself the target and `plan/11_转变模型_可行性与实施阶梯.md` reviews it (five corrections: units in its growth closure, d_b→K_s scaling cannot give crust, spanning degenerates to θ=1 in 2-D, wicking transition data = 2014 bead series, fixed-dt occupancy stepping) and gives the V0–V5 ladder that now replaces S3; `plan/09_BenAmara2025_在子计划中的用法.md` says how the Ben Amara 2025 experiments define the default case (δ = 15 mm, RH 7 %, 21 °C, substrates #1/#4) and which of their data each sub-plan may use; `plan/14_预实验参数汇总.md` lists every parameter of the E1_pre pre-experiment with value, provenance tag (M/L/A/B/V/N) and the calibration freedoms (c_k, δ_eff); `v2_model_formulation.md` is the current mathematical formulation (domains, equations with source tags, gating, PNM interface, parameter list, validation sequence); `v2_盐域参数在体素模型中的实现.md` says how each salt-structure parameter enters the voxel model (birth-time properties per `structure_type`, wet/dry closures; salt–substrate detachment/gap state is deliberately out of scope for now); the others are the literature-derived reasoning behind them. Read them before proposing equations. |
| `knowledge/` | Per-paper knowledge base plus `porous_salt_parameters.md`, the cross-paper compendium of salt-structure parameters (ε_s, pore size, k_s, P_cth, a_v, crust vapour resistance) with provenance tags. **Check it before choosing any parameter or closure**; `README.md` there lists where the papers conflict. |
| `literature/` | The nine papers the model is built on (see the regime table below). |
| `saltv2/` + `tests/` | The model, in five sub-packages (modularised 2026-09-18, move only; `runs/E1_pre/README.md` §14 has the old → new map). `materials/`: the ONLY place a property value is written — `brine.py` (DuMuX port: Batzle-Wang density, D = 1.54e-9, simplified Pitzer a_w; plus MU, GAMMA, RHO_L, C_SAT), `air.py` (D_V, RHO_G, W_SAT, derived RHO_D / X_SALT), `salt.py`, `constants.py`. `pnm/`: the substrate pore network aligned with DuMuX — `network.py` (`PoreLayer`: generated lattice, or `PoreLayer.from_dgf` for a DuMuX dgf; grids in `pnm/grids/`), `flow.py` (`solve_pressure`, mass-conserving, upwind density), `transport.py` (`salt_step`), `precipitation.py` (OPT-IN in-pore precipitation, off on the production path), `boundary.py` (the DuMuX top boundary condition, reference cases only), `dgf.py`. `growth/`: the voxel domain — `grid.py`, `geometry.py`, `regularize.py`, `morphology.py` (no physics), `fv.py` + `vapor.py` (steady gas diffusion, sealed crust as a porous conductor), `crust.py` (per-voxel salt properties, reaction band, WET → SEALED), `supply.py` (hierarchical supply gate), `prat.py` (the one place that puts `prat2026/` on `sys.path`). `coupling/`: `interface.py` (where the pores open into the grid), `ledger.py` (named mass channels), `driver3d.py` (the 3-D closed loop, `run`), `driver2d.py` (the 2-D loop). `cases/`: every value that changes with the experiment — `scan_crust.py` (production table), `dumux1d.py`, `vtp2014_column.py` (a fully resolved 35 mm pore column per 2014 bead class, built with the paper's own network geometry; d_sub = d_b), `small.py` (test boxes), `pnm_column.py` (`python -m saltv2.cases.pnm_column`). `tests/conftest.py` puts the repo root on `sys.path` (nothing else is needed); `tests/test_regression.py` pins the production path bit for bit against `tests/data/ref_*.npz`; `tests/test_layering.py` enforces the dependency rules below. |
| `runs/` | Launchers, reports and execution records only — no model code. `E1_pre/`: `pre_experiment_vapor.py` (CLI launcher of the `SCAN_CRUST` cases), `run_crust.sh` (frodo), `report_3d.py` → `report_3d.html`, `pnm_dumux_compare.py` + `report_pnm_dumux.py` → `report_pnm_dumux.html` (PNM vs DuMuX), `dumux_column_growth.py` → `report_dumux_column_growth.html` (the DuMuX column coupled to the external salt structure; cases in `saltv2/cases/dumux1d.py`), `vtp2014_column.py` + `run_vtp2014_column.sh` (frodo) → `report_vtp2014_column.html` (the resolved 2014 bead column: does the structure over a pore seal, against figure 23), `dumux_reference/` (their inputs), README = execution record; `E2_feedback/`: README of the 2-D driver; `E3_grainsize/`: 2-D grain-size scans (`scan.py`, `scan_eloukabi.py`, `crust_vapour.py`, `front.py`) and `report_full.py` → `report_2d.html`; `exp_2014/`: the 2014 reference sheet; `S0_hierarchy/hier1d.py`. All compute jobs run on frodo, not locally. |
| `prat2026/` | Verbatim copy of `../salt_grow/prat2026`: validated analytic implementation of the Prat (2026) growth law. Reuse `core.py` as the growth-law closure; do not re-derive it. |

Design notes and READMEs are written in Chinese; code, docstrings, and tests are in English.

## Commands

There is no `.venv` here. Use the sibling project's interpreter by absolute path:

```bash
PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python      # numpy, scipy, matplotlib, pillow, pytest

OMP_NUM_THREADS=1 $PY -m pytest tests/ -q                     # ~9 min; 7 files run the drivers.  On frodo, clean clone of 4017046 + the re-captured ref_3d_small (2026-09-19): 154 passed, 1 skipped.  Locally: 151 passed, 1 skipped, 3 failed --
                                                              # the 3 platform-pinned references below.  test_P6::test_active_window_... is MARGINAL: window vs full box differ by ~1e-6 in E_total against rtol 1e-6, so it failed (1.5e-6) before the 2026-09-19 sink fix and passes after it without its cause having been looked at
                                                              # test_regression's references are pinned to frodo's numpy/platform; locally the
                                                              # three small cases differ at the ulp level (run the suite in a clean clone on frodo)
SALTV2_MEDIUM=1 $PY -m pytest tests/test_regression.py -q      # + the medium frodo case (crust_out/ref_scancrust_7.npz vs vapor_crust_7.npz)
$PY runs/E1_pre/pre_experiment_vapor.py scancrust 0 snap      # one production case (hours; run it on frodo: bash runs/E1_pre/run_crust.sh, K="0" selects)
$PY runs/E1_pre/report_3d.py                                  # rebuild runs/E1_pre/report_3d.html from crust_out/*.npz
$PY -m tests.npz_compare ref.npz new.npz                      # exact comparison of two run records
$PY -m saltv2.cases.pnm_column                                # stand-alone pore-network column check (no vapour solve)
$PY runs/E1_pre/pnm_dumux_compare.py --run-dir <dumux run> --top-bc boundary_layer   # PNM vs DuMuX, closed loop
$PY runs/E1_pre/dumux_column_growth.py growth                 # DuMuX column + external salt (check | farfield | growth | growth_shell | report)
$PY runs/E1_pre/vtp2014_column.py run 1mm                     # resolved 2014 bead column (5-50um | 50-62um | 100-125um | 1mm; --c-k, --vpp; report).  Fine classes: bash runs/E1_pre/run_vtp2014_column.sh (frodo, hours)
cd prat2026
$PY -m pytest tests/ -q                                       # 37 tests, all pass (verified 2026-09-01)
$PY -m pytest tests/test_anchors.py::test_growth_rate_anchor_25um -q   # single test
$PY -m pytest tests/ -q -k "Da_star"                          # by keyword
for d in fig0*/; do $PY $d/run.py; done                       # regenerate every figure PNG
$PY fig04_DaStar_vs_dsub/investigate.py                       # Da* / beta discrepancy analysis
```

`prat2026/README.md` writes `../.venv/bin/python`; that path resolves only inside `../salt_grow`.

### Scripts in `prat2026/` whose inputs live only in `../salt_grow`

The copy was made without its neighbours, so these fail here unless run from
`../salt_grow/prat2026` or given corrected paths. Do not silently add fallbacks; fix the path
explicitly or run them in the original tree.

| Script | Expects |
|---|---|
| `digitize.py` | `../reference/Prat - 2026 - ….pdf` (here the PDF is in `literature/`). Its CSV outputs in `data/` are already present, so figure scripts run without it. |
| `multi_c/analyze.py` | `../../new_model/validation/data/nr2013_*.csv` |
| `table1_2014_check.py` | `paper.txt` (plain-text Veran-Tissoires & Prat 2014) and `new_model/` |
| `fig09_diouf_heave/heave_per_water.py` | `new_model/validation` anchors |

## Where code goes (rules of the 2026-09-18 modularisation)

Dependency direction, enforced by `tests/test_layering.py`:
`materials ← pnm`, `materials ← growth`, `{materials, pnm, growth} ← coupling ← cases ← runs/, tests/`.
`pnm` and `growth` do not import each other.

- **R1** `saltv2/` never imports `runs/`, and a kernel's docstring does not cite a file under `runs/` as its evidence (`cases/` may: a case value's provenance is the execution record).
- **R2** Case values do not live in kernels. Test: "would this value change if the experiment changed?" Yes → `cases/`; no → `materials/` or the kernel.
- **R3** Anything imported by two case directories or by `tests/` belongs in `saltv2/`.
- **R4** A property value is written once, in `materials/`; a default elsewhere reads `= materials.<X>`; derived quantities (`RHO_D`, `X_SALT`) are computed once, never stored twice.
- **R5** One case = one module in `cases/`, consistent within itself. No config-file format and no layered defaults (a layered default is what produced the DuMuX 10 °C trap).

Known residue, deliberately NOT fixed by the move because each changes run records (batch them into one production rerun): the production table still takes `d_sub = 4 r` (the driver's default; a case can state it since 2026-09-20, `run(d_sub=)`, and 4 r puts a pore-BODY radius into a relation that holds for throats) and `Lt_res = 22 mm`, which is the 2014 L_s (height of the bead surface over the reservoir level) where the salt transport length is the 35 mm bead column; `materials/air.py` constants are not at one temperature (W_SAT ≈ 20.2 °C, RHO_G ≈ 21.6 °C vs `brine.T_REF` = 21 °C, Pitzer parameters 25 °C) and should become f(T) with T declared by the case; `growth/crust.py` `M_W = 0.018` vs `materials.constants.M_H2O`; `driver2d.run` still carries case defaults (`rh`, `delta_eff_mm`); the λ=32 lateral throat length is clamped at 10 µm (geometric 5.5 µm); the lumped 22 mm reservoir throat is 2.19× more resistive than a resolved column; the pressure solve has no pore storage term and `salt_step` is in non-conservative form, so V·w·dρ of pore salt has no flux behind it (bounded: −0.21 % of the reservoir salt on the DuMuX column, `runs/E1_pre/mass_audit.py`, README E1 §16); a wet cluster no pore feeds keeps β(C_sat) with `salt_unsourced` and its water has no source either (it would really dry out); no global water inventory exists. OPEN modelling decision, not residue: the supply gate keeps Prat eq (33) as written, so the brine of non-growing faces counts for its water only; counting its salt as well costs `ref_3d_small` one voxel layer (−13 %). Salt closure since 2026-09-19 (README E1 §15.5–15.8): a pore under salt is held at C_sat (`pnm.salt_step(crystal=, w_eq=)`, implicit, released when salt diffuses away faster than it arrives); what crystallises there is the structure's ONLY salt, handed to the faces ∝ evaporated water through a per-pore interface buffer, at most β(C_sat) per kg — Prat eq (29) is the upper bound, reached by every reservoir at or above saturation; `salt_crystallised = salt_supplied + salt_nucleated + salt_interface_buffer`. No nucleation supersaturation since 2026-09-19 (`S_NUC = 1.0`, user's decision, §15.9): the first crystal appears at C_sat; the 1.2 of Rodríguez-Rodríguez 2026 / VT&P 2014 needed a one-step reset 1.2 C_sat → C_sat that no source describes. Not yet assessed: what this does to nucleation-site selectivity in the production table (C0 = 0.25).

## The physics you must keep straight: drying vs evaporation-wicking

Every experiment, equation, and validation in this project belongs to one of two regimes.
State the regime explicitly whenever you write a model, compare to data, or cite a paper.

- **Evaporation-wicking**: porous medium connected to a brine reservoir, stays liquid-saturated,
  continuous salt + water supply. No drying front. This is the **planned starting point**
  (simpler; isolates the external salt morphology problem).
- **Drying**: finite initial brine, saturation and concentration fields evolve inside the
  substrate, hydraulic disconnection ends growth. Planned as stage 2.

| Paper (`literature/`) | Regime | What it supplies |
|---|---|---|
| Eloukabi et al. 2013, IJHMT 56 (`1-s2.0-S0017931012007478-main.pdf`) | drying | crusty–patchy transition; capillary-feeding vs viscous-loss argument |
| Veran-Tissoires & Prat 2014, JFM 749 | **wicking** | nucleation sites from local Pe, evaporation screening, patch competition; continuum + PNM + lattice growth |
| Prat 2026, Langmuir (`prat2026/` implements it) | local closure, regime-agnostic | growth law `V_n = α(Da*) · j` (eqs 29/32/36) |
| Ben Amara, Sghaier & Kharaghani 2025, Drying Technology 43(13) | **drying** | quantitative data: mass-loss curves (Figs 4, 9), efflorescence height vs time from structured-light scans (Figs 7, 12), sample tube open at the top with the bead surface recessed 15 mm below the rim (no lid; δ = recess depth), height–layer-thickness power law (Fig 13), Table 1 (bead 1–50 µm / 400–600 µm, pore 6.25 / 125 µm, k = 5.4e-14 / 2.2e-11 m², a_w ≈ 0.75), initial w_s = 25 %, T = 21 °C |
| Shokri-Kuehni, Sahimi & Shokri 2022, Drying Technology 40(4) | perspective | context for prediction approaches |
| Piotrowski et al. 2020, Geosciences 10:423 (`geosciences-10-00423-v2.pdf`) | drying | only direct measurement of salt-layer permeability: NaCl efflorescent layer on 240 µm sand, gas-permeameter k = 4.77e-12 m² (0.23 × sand). This is the *patchy-type* value; the nanoporous crust on 1–50 µm beads is ~10 orders lower |
| Wijnhorst, Van der Sloot, Pel & Shahidzadeh 2024, Phys. Rev. Applied 21:064055 | drying (sandstone cones) | same material, Pe 1.5 → non-blocking thin film, Pe 2.7 → cauliflowers coalesce into a blocking crust; after the crust top loses hydraulic contact, crystallization moves to the crust base and salt "legs" lift the crust (mechanism behind the deferred detachment state); wet/dry cycles of surface crystals |
| Rodríguez-Rodríguez et al. 2026, arXiv 2607.11587 | capillary, no substrate | salt creeping as a relaxation oscillator: heterogeneous nucleation threshold φ_eff/φ_sat ≈ 1.20 vs deliquescence at saturation; episodic imbibition into the salt deposit (rewetting with hysteresis); creeping deposits are dense fine crystals, free radial growth is loose coarse crystals |
| Wijnhorst, Prat & Shahidzadeh 2025, arXiv 2508.18779 (`2508.18779v1 (1).pdf`) | drying (sandstone) / beaker | only quantitative salt-structure microstructure: crystal and pore size halve per layer, pore ≈ crystal/4, ε_s ≈ 0.20, fractal D 1.71 (2D) / 2.56 (3D); hierarchical model showing capillary–viscous competition does **not** limit cauliflower height |

Parameter rule that follows from these papers (verified against Licsandru 2022 in `../salt_grow`): the salt structure's
permeability and capillary threshold must be **independent inputs**. k_s is set by the coarsest internal pores, P_cth by
the finest outer pores; deriving both from one pore size via Carman–Kozeny overestimates k_s by ~3 orders of magnitude
for crusts.

The Ben Amara 2025 data are valuable for validating morphology (crust vs patch), the
evaporation-rate drop caused by a crust (their Fig 4, 1–50 µm beads: a kink at ≈ 20–30 h (its width is not resolved by the sampling, ≲ 8 h),
≈ 1.3 g evaporated, rate down by ≈ 15–20, read off the figure; the "factor 14 for 5–50 µm beads"
in their §3.4 is a sentence of Veran-Tissoires & Prat 2014 §6.3 about THAT paper's figure 23, not
a measurement of this one), and height kinetics (Fig 7: ≈ 0.26 mm by 22 h, then ≈ 2.5 µm/h),
but they are **drying** experiments. A wicking-regime model may be compared to them only for
the early, still-saturated phase, and that limitation must be stated in any comparison.
`ideas/external_salt_deposition_model_notes.md` §17 lists what the paper cannot validate.

## Model architecture (from `ideas/`)

The design decomposes the problem into four modules; each has a designated source paper:

1. **Substrate salt transport** (2014 JFM): advection–diffusion in the porous medium; local
   Pe sets where surface concentration peaks. Nucleation threshold `C > C_nuc > C_sat`.
2. **Nucleation → seed** (2014 JFM): seeds emerge from the concentration field, not by hand.
3. **Gas-side evaporation** (2014 JFM + Laplacian growth): solve vapour transport in the air
   domain with `c_v = a_w(C_eq) c_v,sat` on wet salt; `j = -D_v ∇c_v · n`. Must be re-solved
   every geometry update (geometry ↔ evaporation is two-way).
4. **External porous salt growth** (Prat 2026 via `prat2026/core.growth_rate`): after
   nucleation the structure interior is at `C ≈ C_eq`; growth is `V_n ∝ j` with prefactor
   `F(Da*) = √Da*/(1+√Da*)`. Add a **hydraulic feasibility gate** from the 2013 picture:
   `Π = P_c / (ρ g H + μ U H / k_s)`; growth arrests when the required capillary pressure
   exceeds what the salt pores can supply.

Key modelling commitments already decided in the notes:

- External salt is a **porous domain** (`Ω_PM + Ω_salt + Ω_gas`) with its own `ε_s, k_s,
  p_c,s, D_v,s^eff`, never an impermeable solid layer.
- The active evaporation surface moves to the outer salt surface as the structure grows.
- Numerical representation target: voxel/cellular states `AIR / POROUS_SALT / SUBSTRATE` (or
  level-set) with a continuum vapour solve, sized for millimetre aggregates, not single crystals.
  This representation is chosen for PNM coupling.
- Staged delivery: **V1** conservative wicking model (mass conservation of water and salt,
  growing porous salt layer, moving evaporation surface) before **V2** morphology selection
  (screening, ion redirection, crusty/patchy criterion) and the drying extension.

## `prat2026/` conventions and known traps

- `core.py` is the single source of truth: pure, stateless, SI, numpy-broadcasting functions,
  each docstring names its paper equation. `params.py` carries provenance on every value.
- `tests/test_anchors.py` pins numbers quoted in the paper's text; `tests/test_consistency.py`
  checks limits and numerical stability with no paper numbers. Keep that split when adding tests.
- `RHO = 1200 kg/m³` is **back-calibrated** (paper never states it); `test_growth_rate_anchor_25um`
  is its sole support. `BETA = 2` is the paper's stated value, while its plotted curves imply
  ≈ 1.8; effect on `dh/dt` < 1 %. Any result depending on these must say so.
- Eq (19) in `sinh` form overflows at the paper's own `Da`; use `core.conc_profile` (exponential
  form). Fig 5b needs the unstated `h = 8.4 mm`.
- Compare in the space of physical quantities (`F`, `dh/dt`), never through the ill-conditioned
  inversion `Da* = (F/(1-F))²`.
- Model lines in Fig 8 start at the first experimental point (nucleation delay); Fig 9b passes
  through the origin. Treat them differently.

## Related directories outside this repo

- `../salt_grow/`: git repo with `.venv`, `reference/` (PDFs + `prat2026.txt` + the 2014 paper
  derivation notes), `paper.txt`, `new_model/efflorescence3d` (3-D evaporation + Prat-law
  V0–V3), `validation/data` (Norouzi Rad 2013 digitisations), the 2014 lattice model.
- `../validation/`: QSPNM × Wijnhorst 2-D pure-water validation project.

## Housekeeping

- A Codex CLI config exists at `~/.codex/config.toml`. If you want its MCP servers / commands /
  instructions available here, reply `/import` to scan it, then `/import --yes=<digest>` to
  apply (or run `claude import` from a terminal if `/import` is unavailable).
