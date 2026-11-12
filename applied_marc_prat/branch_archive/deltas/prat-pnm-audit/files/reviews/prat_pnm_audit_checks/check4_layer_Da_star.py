# Run from the repository root with the sibling venv:  PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python; OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/<this file>   (check4 needs PYTHONPATH=prat2026)
"""Counterexample 4: Prat's one-scale Da* applied to the Wijnhorst hierarchy."""
import sys; sys.path.insert(0, ".")
import numpy as np
from saltv2.growth.prat import prat
from saltv2.growth.crust import BirthClosure, C_H, ALPHA, R_HALVING
import params as P
eps = 0.43
print(f"prat2026 params: Ds={P.DS_T21}, kr={P.KR_T21}, beta={P.BETA}, C_eq={P.C_EQ}; crust limit F<=1-eps={1-eps:.2f}")
# crust-limit grain size
Da_lim = ((1 - eps) / eps) ** 2
d_lim = Da_lim / float(prat.Da_star(eps_s=eps, d_sub=1.0))
print(f"Da*_crit = {Da_lim:.3f}  ->  d_sub,crit = {d_lim*1e6:.2f} um (one-scale, layer 1)")
for d in (27.5e-6, 34e-6, 56e-6, 112.5e-6, 1080e-6):
    c = BirthClosure(d_sub=d, eps_s0=eps)
    print(f"\nd_sub = {d*1e6:.1f} um")
    print("  n   d_c[um]  a_v[1/m]   Da*_n     F_n    eps_birth_n  delta_C_n[um]  a_1(n)[um]  h_n[um]  Da(h_n)")
    for n in range(1, 8):
        dc = 0.14 * d * R_HALVING ** (n - 1)
        av = 5 * (1 - eps) / dc
        Da = float(prat.Da_star(eps_s=eps, av=av)); F = float(prat.F(Da))
        dC = float(prat.boundary_layer_thickness(eps_s=eps, av=av))
        hn = C_H * d * ALPHA ** (n - 1)
        eb = 1 - (1 - eps) / F
        a1 = float(c.pore_radius(n))
        print(f"  {n}  {dc*1e6:7.3f}  {av:9.3e}  {Da:8.3f}  {F:6.3f}  {eb:+8.3f}      {dC*1e6:7.3f}      {a1*1e6:7.3f}   {hn*1e6:8.1f}  {hn**2/dC**2:9.2e}")
