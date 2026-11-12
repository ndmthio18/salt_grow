# Run from the repository root with the sibling venv:  PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/<this file>   (check4 needs PYTHONPATH=prat2026)
"""Counterexample 1: eps_birth changes no transport property and no supply decision."""
import sys; sys.path.insert(0, ".")
import numpy as np
from saltv2.growth.crust import BirthClosure, CrustField
from saltv2.growth.geometry import GeometryEngine
from saltv2.growth.grid import Grid, SALT_WET
from saltv2.growth.supply import first_failure_height
from saltv2.growth.prat import prat

d_sub = 34e-6
F = float(prat.F(prat.Da_star(eps_s=0.43, d_sub=d_sub)))
dC = float(prat.boundary_layer_thickness(eps_s=0.43, d_sub=d_sub))
print(f"d_sub={d_sub*1e6:.0f} um  F={F:.4f}  eps_birth={1-(1-0.43)/F:.4f}  delta_C={dC*1e6:.2f} um")
rows = {}
for label, Fk in (("self-replication (F=None)", None), ("task B (F)", F)):
    c = BirthClosure(d_sub=d_sub, eps_s0=0.43, kinetic_efficiency=Fk, reaction_band=dC)
    g = Grid((2, 2, 4), 2.5e-6); eng = GeometryEngine(g); eng.state[..., :2] = SALT_WET
    cf = CrustField(g, eng, c)
    i = (0, 0, 0)
    rows[label] = dict(eps=cf.eps_s[i], k_s=cf.k_s[i], P_cth=cf.P_cth[i], r_p=cf.r_p[i],
                       D_eff=cf.D_eff[i], cap=cf.skin_capacity(0.02)[i])
for k in rows[list(rows)[0]]:
    a, b = rows["self-replication (F=None)"][k], rows["task B (F)"][k]
    print(f"  {k:6s}  F=None {a:.4e}   F={F:.3f} {b:.4e}   ratio {b/a:.4f}")
# supply gate: has no F argument at all
U = 1e-6
for ck in (3e-3, 1e-7):
    print(f"  first_failure_height(U={U}, c_k={ck}) = {first_failure_height(U, d_sub, 0.0, ck, 1e-3)*1e6:.2f} um  (signature has no F, no eps)")
# permeability of the closure does not contain eps at all
c = BirthClosure(d_sub=d_sub, eps_s0=0.43)
print("  k_s(n=1..4) [m^2]:", c.permeability(np.arange(1, 5)))
c2 = BirthClosure(d_sub=d_sub, eps_s0=0.20)
print("  same with eps_s0=0.20:", c2.permeability(np.arange(1, 5)))
