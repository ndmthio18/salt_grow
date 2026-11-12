# Branch index and scientific classification

Reference base: `verify/step2a-prat-single-scale` at
`f7dbf3876df9a563c3572844d48afd886395614f`.

## Branches already contained in the Step 2A ancestry

The main migrated snapshot already contains the committed history of:

- `master` (`cf81ea1`)
- `baseline-selfconsistent` (`72ee682`)
- `surface-first` (`158353d`)
- `verify/step0-single-pore-evaporation` (`a23cb78`)
- `verify/step1-frozen-salt-pillar` (`57903db`)
- `verify/step1-grid-refinement-10um` (`ea48aac`)
- `verify/step2a-prat-single-scale` (`f7dbf38`)

## Independent branch tips

| Route | Tip | Unique commits | Changed files | Classification | Evidence boundary |
|---|---|---:|---:|---|---|
| [`structure-hierarchical`](deltas/structure-hierarchical/BRANCH.md) | `6faf96f` | 2 | 2 | configuration branch | Changes the launcher default and documentation; it is not a separate validated solver. The substantial `layer_av` experiment exists only in the uncommitted overlay. |
| [`structure-single-scale`](deltas/structure-single-scale/BRANCH.md) | `33379c4` | 2 | 2 | configuration/control branch | Changes the launcher default to the single-scale structure; no new physical validation. |
| [`codex/localized-nucleation-seed`](deltas/localized-nucleation-seed/BRANCH.md) | `7179312` | 2 | 5 | records and mechanism audit | Adds the Prat applicability review, framework-defect plan, resolved 2014-column report, and kinetic crust-layer check. It does not repair production code. |
| [`review/prat-pnm-audit`](deltas/prat-pnm-audit/BRANCH.md) | `facdfff` | 3 | 11 | critical implementation/physics audit | Demonstrates bookkeeping consistency but identifies unsourced salt, dead hydraulic fields, hierarchy inconsistency, and a supply-failure/continued-evaporation contradiction. No production fix was made. |
| [`review/prat-layered-validation`](deltas/prat-layered-validation/BRANCH.md) | `d80e791` | 1 | 26 | independent 1-D numerical review | Tests homogeneous and layered variants and exposes strong dependence on the `epsilon_new` versus `epsilon_final` interpretation. This is equation-level verification, not experimental validation. |
| [`review/prat-beta-hierarchy`](deltas/prat-beta-hierarchy/BRANCH.md) | `a43a0a8` | 2 | 49 | beta/hierarchy sensitivity and solver audit | Shows that beta depends on assumed crystal growth geometry and that predicted closure layer is highly closure-dependent. It also identifies a Picard-iterate closure defect for large time steps. |
| [`codex/substrate-wicking-stage-t`](deltas/substrate-wicking-stage-t/BRANCH.md) | `f7b0355` | 5 | 80 | pre-crystallisation transport study | Stage-T is only conditionally passed. The opt-in conservative path closes, but the production path has large water/salt budget defects; lateral-domain and evaporation-demand closure remain open. No crystallisation or morphology Stage-G run was launched. |

## Cross-branch scientific conclusion

The branch history does not provide seven independent confirmations of one
model. It records competing structural assumptions, increasingly critical
audits, and staged verification problems. The strongest reusable evidence is:

- analytical Prat reproduction and one-dimensional reference verification;
- explicit identification of mass-bookkeeping versus constitutive assumptions;
- conservation-aware PNM and coupling test patterns; and
- documented failure modes and non-identifiability.

The branch history does **not** validate the three-dimensional morphology,
hierarchical sealing mechanism, hydraulic growth stop, or patchy–crusty
transition against independent experiments.

