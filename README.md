# Salt growth: dynamic internal precipitation

This repository contains a reproducible one-dimensional comparison of the
Prat (2026) analytical porous-salt growth law with a conservative moving-domain
model that permits internal precipitation and porosity evolution.

- Research instructions: [`instruction/1.md`](instruction/1.md)
- Reproduction guide: [`prat_dynamic_growth/README.md`](prat_dynamic_growth/README.md)
- Results: [`prat_dynamic_growth/RESULTS.md`](prat_dynamic_growth/RESULTS.md)
- Assumptions and limitations:
  [`prat_dynamic_growth/MODEL_ASSUMPTIONS.md`](prat_dynamic_growth/MODEL_ASSUMPTIONS.md)
- GO/NO-GO assessment:
  [`prat_dynamic_growth/GO_NO_GO.md`](prat_dynamic_growth/GO_NO_GO.md)
- M2 critical review and follow-up design:
  [`prat_dynamic_growth/M2_CRITICAL_REVIEW.md`](prat_dynamic_growth/M2_CRITICAL_REVIEW.md)
- M3 controlled attribution and corrected ALE equations:
  [`prat_dynamic_growth/M3_METHODS.md`](prat_dynamic_growth/M3_METHODS.md)
- Current numerical evidence:
  [`prat_dynamic_growth/outputs/m3/`](prat_dynamic_growth/outputs/m3/)
- M4 fixed-material and travelling-profile analysis:
  [`prat_dynamic_growth/M4_MATERIAL_HISTORY.md`](prat_dynamic_growth/M4_MATERIAL_HISTORY.md)
- Independent Phase 1 evaporation and salt mass-budget validation:
  [`docs/phase1_feasibility_and_validation.md`](docs/phase1_feasibility_and_validation.md)
- Phase 1 source, tests, and reproducible reference case:
  [`phase1_mass_budget/`](phase1_mass_budget/)
- Phase 2A minimal CPU 3D CA and deterministic demonstration:
  [`phase2_ca/`](phase2_ca/)
- Phase 2A design, numerical evidence, and limitations:
  [`docs/phase2a_minimal_ca.md`](docs/phase2a_minimal_ca.md)
- Phase 2A.1 stochastic growth-rule and morphology validation:
  [`docs/phase2a_1_growth_rule_validation.md`](docs/phase2a_1_growth_rule_validation.md)
- Phase 2A.3a surface-normal and cubic-orientation verification:
  [`docs/phase2a_3a_normals_orientation_verification.md`](docs/phase2a_3a_normals_orientation_verification.md)

The repository records numerical verification separately from physical
validation. The reported reaction layer is smaller than the estimated crystal
scale, so no claim of physical validation is made.
