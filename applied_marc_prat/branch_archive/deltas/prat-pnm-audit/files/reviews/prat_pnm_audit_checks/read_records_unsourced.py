# Run from the repository root with the sibling venv:  PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/<this file>   (check4 needs PYTHONPATH=prat2026)
import numpy as np, glob, json, ast
for f in sorted(glob.glob("runs/E1_pre/crust_out/vapor_crust_[0-9].npz")):
    d = np.load(f, allow_pickle=True)
    led = d["ledger"]; led = led.item() if led.dtype == object else ast.literal_eval(str(led))
    p = d["params"]; p = p.item() if p.dtype == object else ast.literal_eval(str(p))
    keys = ["salt_grown","salt_front_external","salt_crust_front","salt_base","salt_unsourced","salt_supplied","water_evaporated_on_salt"]
    print(f.split('/')[-1], {k: f"{led.get(k, float('nan')):.3e}" for k in keys})
    print("   ", {k: p.get(k) for k in ("lam","c_k","eps_min","d_sub","t_end","tag")}, "unsourced/grown =", f"{led.get('salt_unsourced',0)/max(led.get('salt_grown',1e-300),1e-300):.2e}")
