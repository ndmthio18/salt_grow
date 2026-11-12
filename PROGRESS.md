# Project progress

## 2026-11-07 — Phase 2A.3b-4R front-error attribution

- Status: complete; stopped at Phase 2A.3b-4R.
- Branch: `codex/phase2a3b4r-front-error-audit`.
- Base: `2e6eb17800311c5a966c5ee5127f2152bbf6fe56`.
- Previous locked decision preserved: `NO_GO_DYNAMIC_PLANAR_FRONT`.
- New decision: conservative mass accounting is retained, but binary-voxel geometry is not sufficient for quantitative interface advance. The minimum next architecture is a conservative subvoxel volume-fraction front with PLIC or verified height-function reconstruction; no replacement was implemented in this phase.
- Principal locked 5 µm attribution: binary quantization contributes −1.173 to −3.169 µm, measurement definition −0.470 to −1.361 µm, and evolving-area allocation −1.447 to +0.975 µm across `<100>/<110>/<111>`.
- Finite-window result: enlarging the growth radius from 40 to 70 µm changes the fixed-interior deposited-field advance by only 0.036–0.057 µm at 5 µm resolution and does not materially reduce orientation spread.
- Numerical integrity: maximum new mass-ledger residual `5.048709793414476e-29 m3`; deferred unallocated volume `0 m3`; 84 tests passed; 20 output hashes verified.
- Deliverable: `docs/phase2a_3b4r_front_error_attribution.md`.
- Outputs: `phase2_ca/outputs/phase2a3b4r_front_error_audit/`.
- Scope held: no multigrain growth, film hydraulics, permeability, vapour transport, PNM coupling, or level-set implementation.

## 2026-11-12 — liqin Marc–Prat implementation migration

- Status: complete; curated historical snapshot copied to `applied_marc_prat/`.
- Primary source: liqin `/home/zhixin/Documents/Humboldt/salt/salt_grow_v2-step2a`, branch `verify/step2a-prat-single-scale`, commit `f7dbf3876df9a563c3572844d48afd886395614f`.
- Copied scope: `prat2026/`, `saltv2/`, tests, the Step 0/1/2A surface-first verification chain, knowledge/design records, and the later Prat applicability review from commit `b7419a24e2151467a3f4757e06cceb22df939237`.
- Integrity: 264 archived files matched the remote Git tree exactly; the separate applicability-review SHA-256 also matched liqin.
- Local core verification: `48 passed` with UTF-8 mode. Full unmodified suite: `269 passed, 1 skipped, 6 failed`; three failures are Windows source-decoding issues and three are bit-for-bit regression mismatches at floating-point roundoff scale. Details are in `applied_marc_prat/PROVENANCE.md`.
- Scientific decision: retain the analytical Prat reproduction, independent 1-D Step 2A reference, conservation-ledger architecture, and gas/PNM/3-D coupling code as reusable evidence and implementation material. Do not treat the imported 3-D morphology, sealing, hydraulic gate, or patchy–crusty behaviour as physically validated; their conclusions depend on unverified porosity, hierarchy, reaction-band, and supply closures.
- Relationship to current work: the migrated snapshot complements rather than replaces `prat_dynamic_growth/`; no historical closure was merged into the current ALE model.

## 2026-11-12 — Complete liqin branch archive

- Status: complete; all 14 local `salt_grow_v2` branch refs and complete reachable history preserved in `applied_marc_prat/branch_archive/salt_grow_v2_all_branches.bundle`.
- Bundle integrity: `git bundle verify` passed; size 48,921,084 bytes; SHA-256 `0F8F44BC80C921269112920793B515D0011145333B8947640B42C87E98C3CA88`.
- Browsable organisation: exported exact path inventories, commits, patches, diffstats, hashes and changed files for the seven tips outside the Step 2A ancestry: hierarchical, single-scale, localized-nucleation records, Prat–PNM audit, layered validation, beta/hierarchy audit and Stage-T substrate wicking.
- Dirty-worktree preservation: `structure-hierarchical` had four modified tracked files plus untracked `tests/test_layer_av.py`; these were captured separately under `branch_archive/uncommitted/structure-hierarchical/` and were not represented as a committed or registered result.
- Restore verification: the bundle restored representative worktrees successfully. Dedicated tests passed for the hierarchical overlay (`2 passed`), layered reference (`12 passed`), beta/hierarchy review (`11 passed`) and Stage-T (`11 passed`).
- Scientific classification: the branches represent competing assumptions, reviews and staged verification—not independent replications. Stage-T remains conditional; Prat–PNM and hierarchy branches preserve identified defects and non-identifiability rather than validating 3-D morphology or patchy–crusty transition physics.

## 2026-11-12 — 5×5 DGF DuMux–Python PNM validation

- Status: implemented and computed; overall `NO-GO`. Deliverables are in `applied_marc_prat/runs/dumux_python_5x5_x008/`; the primary report is `REPORT_ZH.md`.
- Input: 30-pore/45-throat DGF archived byte-for-byte, SHA-256 `E48FF929AEAE9149DF9B575FD2BD612C7C19E111810D023DBF7013EE5341B825`; bottom five pores are the `x=0.08` reservoir and top five pores use the 1.6 mm activity-coupled evaporation boundary.
- Implementation: added a multi-pore Dirichlet DGF reader with 25 dynamic pores, 36 internal throats and five reservoir throats; retained original pore/throat indices; added 5×5 Python driver, per-pore/per-throat output, DuMux-flux replay, water/salt ledgers, normalized CSV/NPZ output and pore-plus-throat pressure plotting. Final combined new/affected plus legacy one-dimensional/conservation/column verification: `40 passed`.
- DuMux provenance: isolated liqin worktree branch `codex/dumux-5x5-x008-validation`; baseline `6dd6d62ae80a0d78d8bae0f652b4d06438427fea`; validation-only Float64 state/inventory and signed phase/component-flux diagnostics commit `ab303806040028cdf9d137dd5b67238b436f2a6f`. The diagnostic patch changes only CSV output, not governing equations. Formal-run binary SHA-256 `A064BC5549C6781B246FF5C6191CB2FEF49AE69323511384D287F5EB64C5C5F0`; diagnostic binary `94865A37D73E04D61D1CB845EB9DBE752F79D84E1E21E3C38F5AD729C34CAF3B`.
- Runtime outcome: Python 10/5/2.5 s runs reached 40000 s. DuMux 5 s aborted near 27340 s after a singular linear system; the 10 s run stalled near 27374 s with `dt=3.08e-173 s` and no forward progress. Raw failures were retained and archived.
- Primary scientific failure: DuMux did not remain in the all-liquid limit—minimum `Sliq=0.970683` at 100 s and about `0.74395` at 900 s—whereas Python is liquid-filled. Its nominally fixed bottom `x=0.08` Dirichlet values also drifted by as much as `3.069e-4` by 1000 s. At 900 s the self-driven and DuMux-flux-replay maximum mole-fraction differences were about `3.56e-3` and `3.72e-3`; refined DuMux solubility crossing was about 934.515 s versus Python 1279.317 s. These results cannot validate a single-phase PNM.
- Conservation: Python salt residual was `7.36e-15` relative at 5 s, but the post-precipitation water residual was 6.150% and did not vanish at 2.5 s; the no-solid 900 s replay water residual was `1.15e-12`. This isolates a non-closing precipitation/pore-volume feedback. The new 1 s DuMux diagnostic reconstructed 0–1000 s maximum relative residuals of 0.196% for water and 3.367% for NaCl, also failing the `1e-6` criterion; the failed formal runs predate these signed diagnostics.
- Figures: full 30-pore/45-throat pressure maps and difference panels are supplied at 940, 1000 and the last common saved time 27300 s, plus a Python-only 40000 s map and pressure-range history, in SVG/PDF/600 dpi PNG with pore/throat CSV. No DuMux 40000 s difference map is claimed because no such state exists.
- Scope: this is a failed all-liquid-limit validation attempt and a useful mismatch diagnosis, not a validation of drainage, moving menisci, unsaturated drying or precipitation feedback.
