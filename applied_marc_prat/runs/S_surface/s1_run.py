"""S1 launcher (plan/19 section 4.1): the validation matrix of saltv2.cases.surface_s1.

    OMP_NUM_THREADS=1 $PY runs/S_surface/s1_run.py        # ~40 min (the four 24 h growth cases dominate)

Writes results/s1_<case>.json (summary + thinned series) and results/s1_summary.md.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from saltv2.cases import surface_s1 as s1                                      # noqa: E402
from saltv2.materials.brine import C_SAT                                       # noqa: E402

RES = Path(__file__).resolve().parent / "results"
LINES = ["# s1_run.py output", ""]


def say(s=""):
    print(s, flush=True); LINES.append(s)


def dump(name, summary, series):
    out = dict(summary=summary, series={k: np.asarray(v).tolist() for k, v in series.items()})
    (RES / f"s1_{name}.json").write_text(json.dumps(out, indent=1, default=float))


def main():
    RES.mkdir(exist_ok=True)
    say("## 1. vapour benchmark: full-plane wet bottom at a_w(C_sat) against the 1-D flux")
    pc = s1.plane_check()
    say(f"  j_1d {pc['j_1d']:.6e}  j_num {pc['j_num']:.6e}  rel err {pc['rel_err']:.1e}  imbalance {pc['imbalance']:.1e}  "
        f"(2014 mean flux J_MEAN {pc['j_mean_2014']:.4e}: the legacy air constants; here moist_air at 21 C)")
    (RES / "s1_plane_check.json").write_text(json.dumps({k: float(v) for k, v in pc.items()}, indent=1))

    say("\n## 2. first saturation without growth: prototype depth vs the 50-layer reference (G2 head 22 mm)")
    say("  n_layers  demand      t_sat [h]   vs D-50     E0 [kg/s]    E(t_sat) [kg/s]  p_req G2 / G1 [Pa]  P_open  res salt / mass / water [kg]")
    ref = {}
    for fb in (False, True):
        for nl in (50, 10, 15, 25):                      # the reference first, so 'vs D-50' is defined
            summ, ser = s1.run(n_layers=nl, feedback=fb, growth=False, t_end_h=4.0, record_s=30.0)
            ts = s1.t_sat_of(ser); ref[(fb, nl)] = ts
            dump(f"nogrowth_{'feedback' if fb else 'prescribed'}_D{nl}", dict(summ, t_sat=ts), ser)
            rel = (ts - ref[(fb, 50)]) / ref[(fb, 50)] if (fb, 50) in ref else np.nan
            say(f"  {nl:4d}      {'feedback ' if fb else 'prescribed'}  {ts/3600:.6f}   {rel:+.4%}   {ser['E'][0]:.4e}   {ser['E'][-1]:.4e}      "
                f"{ser['p_req_G2'][-1]:6.1f} / {ser['p_req_G1'][-1]:6.1f}     {summ['P_open']:.0f}   {summ['res_salt']:+.1e} / {summ['res_mass']:+.1e} / {summ['res_water']:+.1e}")
        for nl in (10, 15, 25):
            ts = ref[(fb, nl)]; say(f"    D-{nl} vs D-50 ({'feedback' if fb else 'prescribed'}): {(ts - ref[(fb, 50)])/ref[(fb, 50)]:+.4%}")
    say("  Reading: the prescribed D-50 value is the kernel cross-validation anchor (1.917636 h with network_step).  The")
    say("  lumped prototype keeps the series resistance of the unresolved throats but drops their pore storage, so it")
    say("  enriches EARLIER (less downward capacity); the Stage-T shallow 10-layer column (reservoir at 7 mm, a perfect")
    say("  sink) was later by +0.25 %.  Convergence with the resolved depth is tabulated above.")

    say("\n## 3. growth: single-scale Prat structure after nucleation at C_sat")
    say("  case                  t_nuc [h]  m_seed [kg]  h(end) [um]   v_end [nm/s]  v/v_prat  X/E (end)  conv C_eq/(1-C_eq)  E(end)/E0  res salt / mass / water  wall [s]")
    conv = C_SAT / (1.0 - C_SAT)
    for name, kw in (("growth_D50_C025", dict(n_layers=50, c0=0.25, t_end_h=24.0)),
                     ("growth_D25_C025", dict(n_layers=25, c0=0.25, t_end_h=24.0)),
                     ("growth_D15_C025", dict(n_layers=15, c0=0.25, t_end_h=24.0)),
                     ("growth_D10_C025", dict(n_layers=10, c0=0.25, t_end_h=24.0)),
                     ("growth_D50_Csat", dict(n_layers=50, c0=C_SAT, t_end_h=6.0)),
                     ("growth_D10_Csat", dict(n_layers=10, c0=C_SAT, t_end_h=6.0))):
        t0 = time.time()
        summ, ser = s1.run(feedback=True, growth=True, record_s=600.0, **kw)
        dump(name, summ, ser)
        ratio = ser["v"][-1] / ser["v_prat"][-1] if ser["v_prat"][-1] > 0 else np.nan
        say(f"  {name:20s}  {summ['t_nuc']/3600 if summ['t_nuc'] else np.nan:8.4f}   {summ['m_seed'] if summ['m_seed'] else np.nan:.3e}   "
            f"{summ['h_end']*1e6:9.2f}    {ser['v'][-1]*1e9:8.4f}    {ratio:.5f}   {ser['X_rate'][-1]/ser['E'][-1]:.5f}   {conv:.5f}            "
            f"{ser['E'][-1]/ser['E'][0]:.4f}   {summ['res_salt']:+.1e} / {summ['res_mass']:+.1e} / {summ['res_water']:+.1e}  {summ['wall_s']:.0f}")
    say("  Reading: with the reservoir at C_sat the column delivers exactly C_eq/(1-C_eq) per unit of evaporated water and")
    say("  dh/dt equals prat2026.core.growth_rate (eq 29) to round-off; with C0 = 0.25 the column delivers less (X/E below")
    say("  the conversion) and the structure grows slower than eq 29 by that ratio -- a supply statement, not a kinetic one.")
    say("  E(end)/E0 > 1 is the far-field approach of a growing wet structure (plan/19 section 3.4), not screening.")
    (RES / "s1_summary.md").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
