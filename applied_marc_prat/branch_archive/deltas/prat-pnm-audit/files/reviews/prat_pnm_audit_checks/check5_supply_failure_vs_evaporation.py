# Run from the repository root with the sibling venv:  PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/<this file>   (check4 needs PYTHONPATH=prat2026)
"""Counterexample 5: a supply-limited cluster keeps evaporating at the full vapour-field flux."""
import sys; sys.path.insert(0, ".")
import numpy as np
from saltv2.coupling import driver3d as v
from saltv2.cases.small import SMALL_3D
kw = dict(SMALL_3D, ell_over_dsub=0.25, eps_min=0.05, t_end=12 * 3600.0)
for ck in (1e-4, 1e-7, 1e-9):
    r = v.run(c_k=ck, **kw)
    h = r["hist"]; c = lambda n: v.col(h, n)
    led = r["ledger"]
    i_fr = int(np.argmax(c("n_frozen") > 0)) if (c("n_frozen") > 0).any() else -1
    print(f"\nc_k={ck:.0e}: steps={len(h)} t_end={c('t_h')[-1]:.2f} h  first frozen row {i_fr} "
          f"(t={c('t_h')[i_fr]:.2f} h)" if i_fr >= 0 else f"\nc_k={ck:.0e}: never frozen")
    print(f"  h_freeze (reach, m) at end: {r['h_freeze']}   a={r['a']*1e6:.2f} um  h_max={c('h_max')[-1]*1e6:.1f} um")
    print(f"  E_total: at freeze {c('E_total')[max(i_fr,0)]:.3f}  end {c('E_total')[-1]:.3f}  max {c('E_total').max():.3f}")
    print(f"  V_sealed end {c('V_sealed')[-1]:.3e}  f_blocked {c('f_blocked')[-1]:.2f}")
    for k in ("water_evaporated", "water_evaporated_on_salt", "salt_released_by_evaporation",
              "salt_grown", "salt_front_external", "salt_crust_front", "salt_supplied", "salt_unsourced"):
        print(f"  {k:30s} {led[k]:.4e}")
    print(f"  residual salt_network {led['relative_residual_salt_network']:.2e}")
