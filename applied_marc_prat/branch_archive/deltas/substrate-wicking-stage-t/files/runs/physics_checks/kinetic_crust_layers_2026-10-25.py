"""Kinetic crust limit of the Prat 2026 law under hierarchical refinement.

F = sqrt(Da*)/(1+sqrt(Da*)) is the fraction of the supplied salt that advances the front;
(1-F) precipitates on the internal surface a_v of the reaction band.  Birth porosity
eps_birth = eps_s - (1-F)(1-eps_s)/F  (project closure, driver3d) reaches eps_min when
F = F_min = (1-eps_s)/(1-eps_min).  With a_v of layer n = a_v1 / r^(n-1) (crystals r x per
layer), Da*_n = Da*_1 r^(n-1).  Layer n spans height h_1 alpha^(n-1).
"""
import sys, numpy as np
sys.path.insert(0, 'prat2026')
import core, params

EPS_MIN = 0.02
D_S = 1.54e-9                      # project value (network D_NACL)
classes = {'5-50um': 27.5e-6, '50-62um': 56e-6, '100-125um': 112.5e-6, '1mm': 1e-3}
growth_um_h = {'5-50um': 13.3, '50-62um': 15.1, '100-125um': 15.9, '1mm': 19.2}   # baseline hier. h_max/t
t_nuc_h = {'5-50um': 0.65, '50-62um': 0.66, '100-125um': 0.68, '1mm': 1.0}
exp_drop_h = {'5-50um': 12.5, '50-62um': 67.0, '100-125um': None, '1mm': None}

def dc_crit(eps_s):
    """Crystal size at which F = 1 - eps_s (kinetic crust limit), from Da* with a_v = 5(1-eps)/d_c."""
    Da_c = ((1 - eps_s) / eps_s) ** 2
    Ds_eff = core.Ds_eff(eps_s, D_S)
    kr = params.KR_T21 if hasattr(params, 'KR_T21') else core.KR_T21
    # Da* = beta^2 (1-eps)^2 (1-C) kr / (eps Ds_eff a_v);  a_v = 5(1-eps)/d_c
    a_v = params.BETA ** 2 * (1 - eps_s) ** 2 * (1 - core.C_EQ) * kr / (eps_s * Ds_eff * Da_c)
    return 5 * (1 - eps_s) / a_v

for eps_s in (0.20, 0.43):
    print(f"eps_s = {eps_s}: F_crust = {1-eps_s:.2f}, Da*_c = {((1-eps_s)/eps_s)**2:.2f}, "
          f"critical crystal size d_c = {dc_crit(eps_s)*1e6:.2f} um  (d_sub = d_c/0.14 = {dc_crit(eps_s)/0.14*1e6:.0f} um)")

print()
for eps_s, label in ((0.20, 'HIERARCHICAL'), ):
    F_min = (1 - eps_s) / (1 - EPS_MIN)
    Da_min = (F_min / (1 - F_min)) ** 2
    print(f"{label}: F_min = {F_min:.3f}, Da*_min = {Da_min:.1f}")
    for alpha, r, c_h, where in ((2.0, 0.5, 1.0, 'bottom'), (2.0, 0.5, 1.0, 'top'), (2.0, 0.5, 5.0, 'bottom'), (2.0, 0.6, 1.0, 'bottom'), (3.0, 0.5, 1.0, 'bottom'), (1.5, 0.5, 1.0, 'bottom')):
        expo = 1 + np.log(alpha) / np.log(1 / r)
        print(f"\n  alpha = {alpha}, r = {r}, c_h = {c_h}, trigger at layer {where}: h_trig ~ d^{expo:.2f}")
        print(f"  {'class':>10} {'Da*_1':>7} {'F_1':>6} {'n_trig':>6} {'h_trig[mm]':>10} {'t_trig[h]':>9} {'exp drop[h]':>11}")
        for name, d in classes.items():
            Da1 = float(core.Da_star(eps_s=eps_s, d_sub=d, Ds=D_S))
            n = 1
            while Da1 * r ** (n - 1) > Da_min:
                n += 1
            # bottom of layer n
            m = n - 1 if where == 'bottom' else n
            h_trig = c_h * d * (alpha ** m - 1) / (alpha - 1)
            t = t_nuc_h[name] + h_trig * 1e6 / growth_um_h[name]
            e = exp_drop_h[name]
            print(f"  {name:>10} {Da1:7.1f} {core.F(Da1):6.3f} {n:6d} {h_trig*1e3:10.3f} {t:9.1f} {str(e) if e else '>142':>11}")
