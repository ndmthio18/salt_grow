# Run from the repository root with the sibling venv:  PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/<this file>   (check4 needs PYTHONPATH=prat2026)
"""Counterexample 3: a seed placed over inert substrate (no pore under it) evaporates and
precipitates for the whole run with neither water nor salt source; booked as salt_unsourced."""
import sys; sys.path.insert(0, ".")
import numpy as np
from saltv2.coupling import driver3d as v
from saltv2.coupling import interface
from saltv2.cases.small import SMALL_3D

orig = interface.PoreOpenings.nucleation_seed
def off_pore_seed(self, i):
    """the lateral cell farthest from every pore footprint (pore_id < 0)"""
    free = self.pore_id < 0
    if not free.any():
        return orig(self, i)
    out = np.zeros_like(free)
    idx = np.argwhere(free)
    out[tuple(idx[i % len(idx)])] = True
    return out
interface.PoreOpenings.nucleation_seed = off_pore_seed
v.MassLedger.check = lambda self, name=None: None   # the off-pore seed is not debited from any pore (per_pore sums footprints), which the interface balance catches; silenced ONLY for this counterexample
r = v.run(c_k=1.0, ell_over_dsub=0.25, eps_min=0.05, t_end=4 * 3600.0, pitch_factor=10.0, **SMALL_3D)
h = r["hist"]; c = lambda n: v.col(h, n); led = r["ledger"]
print(f"steps={len(h)}  t_end={c('t_h')[-1]:.2f} h  n_clusters end {int(c('n_clusters')[-1])}  n_frozen end {int(c('n_frozen')[-1])}")
print(f"h_max {c('h_max')[-1]*1e6:.1f} um  V_sealed {c('V_sealed')[-1]:.3e}  E_total end {c('E_total')[-1]:.3f}")
for k in ("water_evaporated", "water_evaporated_on_salt", "salt_released_by_evaporation",
          "salt_grown", "salt_front_external", "salt_crust_front", "salt_base",
          "salt_supplied", "salt_unsourced", "salt_crystallised", "salt_nucleated", "salt_interface_buffer"):
    print(f"  {k:30s} {led[k]:.4e}")
print(f"residuals: network {led['relative_residual_salt_network']:.2e} interface {led['relative_residual_salt_interface']:.2e} precipitation {led['relative_residual_salt_precipitation']:.2e}")
print("M_liquid_in (reservoir supply, kg/s) max:", float(np.max(r['M_liquid_in'])), " vs water_evaporated rate ~", led['water_evaporated']/(c('t_h')[-1]*3600))
