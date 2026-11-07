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
