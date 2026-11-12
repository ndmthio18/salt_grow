"""The self-consistent baseline (``saltv2/cases/baseline.py``): launcher, the six consistency
criteria per run, figures and report.  No model code.

Regime: evaporation-wicking (Dirichlet reservoir under a saturated, fully resolved pore column;
one surface pore in a periodic box).  The question is NOT whether the runs match figure 23 but
whether the coupled model is consistent with itself:

  1  the global water and salt books close                       (ledger balances of the run)
  2  k, P_c, D_eff, eps and a_v of the salt come from ONE microstructure   (closure record)
  3  every state transition is triggered by model variables      (freeze / seal events, gate margin)
  4  closures are used inside their validity, or the violation is declared
  5  limits return known solutions                                (quasi-steady column profile)
  6  what is reported converges with the mesh                     (vpp 7 against 14)

    python runs/E1_pre/baseline.py run 1mm single_scale
    python runs/E1_pre/baseline.py run 50-62um hierarchical --vpp 14
    python runs/E1_pre/baseline.py run dumux1d single_scale
    python runs/E1_pre/baseline.py report

Records and figures go to the ignored ``simulation_outputs/baseline_<case>_<structure>_vpp<vpp>/``.
The fine classes take hours: ``bash runs/E1_pre/run_baseline.sh`` (frodo).
"""
import argparse
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (str(ROOT), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

from saltv2.cases import baseline as case, dumux1d, vtp2014_column as vtp  # noqa: E402
from saltv2.coupling.driver3d import col, run  # noqa: E402
from saltv2.growth.crust import HIERARCHICAL, SINGLE_SCALE, BirthClosure  # noqa: E402
from saltv2.growth.prat import prat  # noqa: E402
from saltv2.growth.supply import supply_need  # noqa: E402
from saltv2.materials.brine import C_SAT, D_NACL, RHO_L  # noqa: E402
from saltv2.materials.constants import G  # noqa: E402
from dumux_column_growth import _img, frames_by_salt_volume  # noqa: E402
from pre_experiment_vapor import _code_version as code_version  # noqa: E402
from report_3d import experiment  # noqa: E402
from single_pore_growth_figures import plot_fluxes, plot_growth  # noqa: E402
from vtp2014_column import LABEL as BEAD_LABEL, plot_column_concentration  # noqa: E402

OUT = HERE / "simulation_outputs"
RECORD = "baseline.npz"
STRUCTURES = {s.name: s for s in (SINGLE_SCALE, HIERARCHICAL)}
CASES = list(case.BEADS) + ["dumux1d"]
LABEL = dict(BEAD_LABEL, dumux1d="DuMuX 柱 (600 µm)")
STRUCT_LABEL = {"single_scale": "单尺度 r=1, ε=0.43", "hierarchical": "分层 r=0.5, ε=0.20"}
N_TOP = 8
PROFILE_EVERY_H = 1.0


def outdir(name, structure, vpp, layer_av=False):
    return OUT / f"baseline_{name}_{structure}{'_layerav' if layer_av else ''}_vpp{vpp}"


def quasi_steady_profile(layer, Q, C0):
    """Closed form of the settled column under a surface pore held at C_sat: with identical
    throats of diffusive conductance g carrying the volume flow Q, upwind advection-diffusion
    gives a geometric profile, C_k = C0 + (C_sat - C0) (r^k - 1) / (r^n - 1), r = 1 + Q/g, pore
    k counted from the reservoir.  (Density is taken constant: Q varies ~1 % along the column.)"""
    d_lat, d_res = layer.diffusive()
    g = float(d_lat[0])
    if not (np.allclose(d_lat, g) and np.isclose(d_res[layer.bottom][0], g)):
        return None
    k = np.arange(1, layer.N + 1)
    log_r = np.log1p(Q / g)
    return C0 + (C_SAT - C0) * np.expm1(k * log_r) / np.expm1(layer.N * log_r)


def run_case(name, structure, vpp, snap_hours, t_end_h, layer_av=False):
    """``layer_av``: B' (2026-09-22) -- per-layer F and dense birth at the clamp; the record goes
    to its own directory (suffix ``_layerav``) so that the baseline records stay what they are."""
    st = STRUCTURES[structure]
    out = outdir(name, structure, vpp, layer_av); out.mkdir(parents=True, exist_ok=True)
    if name == "dumux1d":
        kw = case.dumux_column(st); kw["vpp"] = vpp
        d_sub = dumux1d.D_SUB
    else:
        kw = case.case(name, vtp.write_dgf(name, out), st, vpp=vpp)
        d_sub = case.BEADS[name]
    if t_end_h is not None:
        kw["t_end"] = t_end_h * 3600.0
    tic = time.time()
    if layer_av:
        kw["layer_av"] = True; kw["tag"] += "_layerav"
    res = run(**kw, snap_hours=snap_hours, verbose=True)
    wall = time.time() - tic
    h = res["hist"]; t_h = col(h, "t_h"); led = res["ledger"]; layer = res["layer"]
    order = np.argsort(-res["pore_z"])                       # surface pore first
    rows = np.unique(np.r_[np.searchsorted(t_h, np.arange(0.0, t_h[-1], PROFILE_EVERY_H)), len(t_h) - 1])

    # criterion 3: how far the supply gate is from failing, at the end of the run.  The Darcy
    # demand is the liquid the reservoir delivers over the opening; the footprint suction is the
    # viscous drop of that flow over the column (single column: throats in series).
    g_lat, g_res = layer.conductances()
    Q_end = float(res["M_liquid_in"][-1]) / RHO_L            # [m^3/s]
    s_foot = Q_end * (float(np.sum(1.0 / g_lat)) + float(1.0 / g_res[layer.bottom][0]))
    A_open = float(layer.A_open[layer.surface].sum())
    h_end = float(col(h, "h_max")[-1])
    need, Pc = supply_need(max(h_end, 1e-9), Q_end / A_open, d_sub, s_foot, kw["c_k"], st)
    # criterion 4: the reaction band against the crystal it is supposed to sit in, and F against
    # the crust limit of the closure
    kr = res["params"]["kr"]
    F = float(prat.F(prat.Da_star(eps_s=st.eps_s, d_sub=d_sub, Ds=D_NACL, kr=kr)))
    band = float(prat.boundary_layer_thickness(eps_s=st.eps_s, d_sub=d_sub, Ds=D_NACL, kr=kr))
    # criterion 5: the settled column against its closed form
    closed = quasi_steady_profile(layer, Q_end, kw["C0"])
    C_end = res["pore_C"][-1][np.argsort(res["pore_z"])]     # reservoir end first
    profile_dev = (float(np.max(np.abs(C_end - closed))) / (C_SAT - kw["C0"])
                   if closed is not None and np.isfinite(res["t_nuc"]).all() else np.nan)
    closure = BirthClosure(d_sub=d_sub, structure=st, c_k=kw["c_k"])
    n_end = int(closure.generation(np.array(h_end)))
    rec = dict(hist=h, hist_columns=np.asarray(res["hist_columns"]),
               pore_C_top=res["pore_C"][:, order[:N_TOP]], pore_z_top=res["pore_z"][order[:N_TOP]],
               profile_t_h=t_h[rows], profile_C=res["pore_C"][rows][:, order], profile_z=res["pore_z"][order],
               M_liquid_in=res["M_liquid_in"], t_nuc=res["t_nuc"], a=res["a"], M_ref=res["M_ref"],
               snapshots_t_h=np.asarray([s[0] for s in res["snaps"]]),
               snapshots_state=np.asarray([s[1] for s in res["snaps"]]),
               C0=kw["C0"], d_sub=d_sub, n_pores=layer.N, n_steps=res["n_steps"], wall_s=wall,
               stop_reason=np.asarray(res["stop_reason"]),
               gate_need=need, gate_Pc=Pc, gate_s_foot=s_foot, gate_gravity=RHO_L * G * h_end,
               gate_structure=need - s_foot - RHO_L * G * max(h_end, 1e-9), generation_end=n_end,
               k_first=float(closure.permeability(np.array(1.0))), k_end=float(closure.permeability(np.array(float(n_end)))),
               Pc_first=float(closure.capillary_threshold(np.array(1.0))),
               F=F, F_limit=1.0 - st.eps_s, band=band, crystal=0.14 * d_sub,
               eps_birth_mean=float(np.mean(res["crust"].eps_birth[np.isfinite(res["crust"].t_birth)]))
               if np.isfinite(res["crust"].t_birth).any() else np.nan,
               profile_dev=profile_dev,
               ledger=np.asarray(repr(led)), closure=np.asarray(repr(res["closure"])),
               resolution=np.asarray(repr(res["resolution"])), params=np.asarray(repr(res["params"])),
               command=np.asarray(" ".join(sys.argv)), code_version=np.asarray(code_version()))
    np.savez_compressed(out / RECORD, **rec)
    title = f"{LABEL[name]}, {st.name}"
    plot_column_concentration(rec, out / "column_concentration.png", title)
    plot_fluxes(res, out / "evaporation_paths.png", title=f"{title}: evaporation by source path")
    frames = frames_by_salt_volume(res["snaps"])
    if frames:
        k_top = max(int(np.argwhere(s.any(axis=(0, 1))).max()) for _, s in frames) + 3
        plot_growth(dict(res, snaps=[(t, s[:, :, :k_top]) for t, s in frames]), out / "growth_3d.png",
                    title=f"{title}: salt above the surface pore")
    print(f"{name} {structure} vpp={vpp}: {res['n_steps']} steps, {wall:.0f} s, stop = {res['stop_reason']}, "
          f"t_nuc = {res['t_nuc'][0]:.0f} s, h_max = {h_end * 1e3:.3f} mm (layer {n_end}), "
          f"E/E0 = {col(h, 'E_total')[-1]:.4f}, frozen max {col(h, 'n_frozen').max():.0f}, "
          f"sealed cells {col(h, 'n_sealed_cells').max():.0f}, gate need/Pc = {need / Pc:.2e}")
    print("  residuals:", {k[18:]: f"{v:+.1e}" for k, v in led.items() if k.startswith("relative_residual")})


# ---- report -------------------------------------------------------------------------------------

def _load(name, structure, vpp):
    f = outdir(name, structure, vpp) / RECORD
    return dict(np.load(f, allow_pickle=False)) if f.exists() else None


def _hist(rec):
    return {str(c): rec["hist"][:, i] for i, c in enumerate(rec["hist_columns"])}


def _led(rec):
    return eval(str(rec["ledger"]), {"nan": np.nan, "inf": np.inf})


def plot_E(recs, output):
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4))
    style = {("single_scale", case.VPP): ("#1f6f8b", "-"), ("hierarchical", case.VPP): ("#c0392b", "-"),
             ("single_scale", case.VPP_FINE): ("#1f6f8b", ":"), ("hierarchical", case.VPP_FINE): ("#c0392b", ":")}
    for ax, name in zip(axes.ravel(), case.BEADS):
        ex = experiment(name)
        ax.plot(ex["t"], ex["E"], "o", ms=3, mfc="none", color="0.25", label="experiment (fig. 23, local rate)")
        for (structure, vpp), (color, ls) in style.items():
            rec = recs.get((name, structure, vpp))
            if rec is not None:
                hh = _hist(rec)
                ax.plot(hh["t_h"], hh["E_total"], color=color, ls=ls, lw=1.6, label=f"{structure}, vpp {vpp}")
        ax.axhline(1.0, color="0.7", lw=0.8)
        ax.set_title(f"{BEAD_LABEL[name]} (paper: {vtp.OBSERVED[name]})"); ax.set_xlabel("time [h]")
        ax.set_ylabel("evaporation rate / bare-surface rate [-]"); ax.grid(alpha=0.25); ax.legend(fontsize=7)
        ax.set_xlim(0, 150); ax.set_ylim(0, None)
    fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)


def plot_height(recs, output):
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    colors = dict(zip(CASES, plt.cm.viridis(np.linspace(0.0, 0.9, len(CASES)))))
    for (name, structure, vpp), rec in recs.items():
        if vpp != case.VPP and name != "dumux1d":
            continue
        hh = _hist(rec)
        ax.plot(hh["t_h"], hh["h_max"] / float(rec["d_sub"]), color=colors[name],
                ls="-" if structure == "single_scale" else "--", lw=1.5, label=f"{LABEL[name]}, {structure}")
    ax.set_xlabel("time [h]"); ax.set_ylabel("h_max / d_sub [-]"); ax.set_yscale("log"); ax.grid(alpha=0.25)
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(output, dpi=180); plt.close(fig)


def delta_mm(name):
    """Far-field distance of the case [mm]."""
    return dumux1d.DELTA_BL * 1e3 if name == "dumux1d" else case.DELTA_EFF_MM


def build_report(output):
    recs = {}
    for name in CASES:
        for structure in STRUCTURES:
            for vpp in sorted({case.VPP, case.VPP_FINE, dumux1d.RUN_3D_GROWTH["vpp"]}):
                rec = _load(name, structure, vpp)
                if rec is not None:
                    recs[(name, structure, vpp)] = rec
    if not recs:
        sys.exit("no record found: run the cases first")
    fig_dir = OUT / "baseline_report"; fig_dir.mkdir(parents=True, exist_ok=True)
    plot_E(recs, fig_dir / "E_vs_experiment.png"); plot_height(recs, fig_dir / "height.png")
    key = lambda kv: (CASES.index(kv[0][0]), kv[0][1], kv[0][2])

    books = micro = trans = valid = ""
    for (name, structure, vpp), rec in sorted(recs.items(), key=key):
        led = _led(rec); hh = _hist(rec); head = f"<td>{LABEL[name]}</td><td>{structure}</td><td>{vpp}</td>"
        orphan = led["salt_unsourced"] * (1.0 + (1.0 - C_SAT) / C_SAT)
        used = (led["water_evaporated"] + led["salt_released_by_evaporation"] + led["salt_released_by_base"]
                + led["salt_brine_in_pores"] / C_SAT)
        books += (f"<tr>{head}<td>{led['relative_residual_salt_pores']:+.1e}</td><td>{led['relative_residual_salt_interface']:+.1e}</td>"
                  f"<td>{led['relative_residual_salt_network']:+.1e}</td><td>{led['relative_residual_salt_precipitation']:+.1e}</td>"
                  f"<td>{np.abs(hh['res_salt']).max():.1e}</td><td>{led['relative_residual_liquid_network']:+.1e}</td>"
                  f"<td>{(led['liquid_out'] - used + orphan) / led['liquid_out']:+.1e}</td>"
                  f"<td>{led['salt_unsourced'] / max(led['salt_reservoir_in'], 1e-300):.2%}</td>"
                  f"<td>{led['salt_interface_buffer'] / max(led['salt_crystallised'], 1e-300):+.2%}</td>"
                  f"<td>{led['liquid_nucleation_jump'] / led['liquid_in']:+.1e}</td></tr>")
        micro += (f"<tr>{head}<td>{float(rec['crystal']) * 1e6:.2f}</td><td>{float(rec['crystal']) / 4 * 1e6:.3f}</td>"
                  f"<td>{float(rec['k_first']):.2e}</td><td>{float(rec['k_end']):.2e}</td><td>{float(rec['Pc_first']) / 1e3:.1f}</td>"
                  f"<td>{float(rec['eps_birth_mean']):.3f}</td><td>{int(rec['generation_end'])}</td></tr>")
        trans += (f"<tr>{head}<td>{float(np.atleast_1d(rec['t_nuc'])[0]):.0f}</td><td>{hh['n_frozen'].max():.0f}</td>"
                  f"<td>{hh['n_sealed_cells'].max():.0f}</td><td>{float(rec['gate_s_foot']):.2e}</td>"
                  f"<td>{float(rec['gate_gravity']):.2e}</td><td>{float(rec['gate_structure']):.2e}</td><td>{float(rec['gate_Pc']):.2e}</td>"
                  f"<td>{float(rec['gate_need']) / float(rec['gate_Pc']):.1e}</td><td>{hh['h_max'][-1] * 1e3:.3f}</td>"
                  f"<td>{hh['E_total'][-1]:.3f}</td><td>{hh['t_h'][-1]:.0f} h · {str(rec['stop_reason'])}</td></tr>")
        valid += (f"<tr>{head}<td>{float(rec['F']):.4f}</td><td>{float(rec['F_limit']):.2f}</td>"
                  f"<td>{(float(rec['F']) / float(rec['F_limit']) - 1):+.1%}</td><td>{float(rec['band']) * 1e6:.3f}</td>"
                  f"<td>{float(rec['band']) / float(rec['crystal']):.3f}</td><td>{float(rec['band']) / float(rec['a']):.3f}</td>"
                  f"<td>{float(rec['profile_dev']):.2%}</td>"
                  f"<td>{hh['E_total'][-1]:.3f} / {delta_mm(name) / (delta_mm(name) - hh['h_max'][-1] * 1e3):.3f}</td></tr>")

    conv = ""
    for name in CASES:
        for structure in STRUCTURES:
            a, b = recs.get((name, structure, case.VPP)), recs.get((name, structure, case.VPP_FINE))
            if a is None or b is None:
                continue
            ha, hb = _hist(a), _hist(b); t_c = min(ha["t_h"][-1], hb["t_h"][-1])
            at = lambda hh, c: float(np.interp(t_c, hh["t_h"], hh[c]))
            la, lb = _led(a), _led(b)
            same_end = abs(ha["t_h"][-1] - hb["t_h"][-1]) < 1e-6        # whole-run integrals compare
            integral = (lambda k: f"{lb[k] / la[k] - 1:+.2%}") if same_end else (lambda k: "—")
            conv += (f"<tr><td>{LABEL[name]}</td><td>{structure}</td><td>{t_c:.0f}</td>"
                     f"<td>{float(np.atleast_1d(a['t_nuc'])[0]):.0f} / {float(np.atleast_1d(b['t_nuc'])[0]):.0f}</td>"
                     f"<td>{at(ha, 'E_total'):.4f} / {at(hb, 'E_total'):.4f}</td>"
                     f"<td>{integral('water_evaporated')}</td><td>{integral('salt_supplied')}</td>"
                     f"<td>{at(ha, 'h_max') * 1e3:.3f} / {at(hb, 'h_max') * 1e3:.3f} ({at(hb, 'h_max') / at(ha, 'h_max') - 1:+.1%})</td>"
                     f"<td>{at(hb, 'V_wet') / at(ha, 'V_wet') - 1:+.1%}</td></tr>")

    figs = ""
    for (name, structure, vpp), rec in sorted(recs.items(), key=key):
        if vpp == case.VPP_FINE:
            continue
        d = outdir(name, structure, vpp)
        figs += f"<h3>{LABEL[name]}，{STRUCT_LABEL[structure]}</h3>"
        for f in ("growth_3d.png", "evaporation_paths.png", "column_concentration.png"):
            if (d / f).exists():
                figs += _img(d / f)
    any_rec = next(iter(recs.values()))
    default_vpp = lambda name: dumux1d.RUN_3D_GROWTH["vpp"] if name == "dumux1d" else case.VPP
    missing = [f"{LABEL[n]} · {st}" for n in CASES for st in STRUCTURES if (n, st, default_vpp(n)) not in recs]
    pending = (f'<p class="note"><b>中间版本</b>：以下算例还没有记录，表和图里没有它们：{"；".join(missing)}。</p>'
               if missing else "")
    html = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>自洽基线</title>
<style>
:root {{ --bg:#fbfbf9; --fg:#1d2327; --mut:#5d666d; --line:#d9dcd6; --acc:#1f6f8b; --warn:#9a5b00; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#15191c; --fg:#e6e8e6; --mut:#9aa4aa; --line:#2e363b; --acc:#6fb6cf; --warn:#e0a458; }} }}
body {{ background:var(--bg); color:var(--fg); font:16px/1.65 "Noto Sans CJK SC","Source Han Sans SC",system-ui,sans-serif; margin:0; }}
main {{ max-width:66rem; margin:0 auto; padding:2rem 1rem 4rem; }}
h1 {{ font-size:1.6rem; margin:0 0 .3rem; }} h2 {{ font-size:1.2rem; margin:2.4rem 0 .6rem; border-bottom:1px solid var(--line); padding-bottom:.3rem; }}
h3 {{ font-size:1.02rem; margin:1.6rem 0 .3rem; }}
.sub {{ color:var(--mut); font-size:.9rem; }} code {{ font-size:.88em; }}
img {{ max-width:100%; background:#fff; border:1px solid var(--line); margin:.6rem 0; }}
.tbl {{ overflow-x:auto; }} table {{ border-collapse:collapse; font-size:.84rem; font-variant-numeric:tabular-nums; }}
th,td {{ border-bottom:1px solid var(--line); padding:.3rem .55rem; text-align:right; white-space:nowrap; }}
th:first-child,td:first-child,th:nth-child(2),td:nth-child(2) {{ text-align:left; }}
.note {{ border-left:3px solid var(--warn); padding:.2rem 0 .2rem .9rem; }} .key {{ border-left:3px solid var(--acc); padding:.2rem 0 .2rem .9rem; }}
</style></head><body><main>
<h1>自洽基线：蒸发–芯吸–结晶–生长耦合，不放任何为成壳而设的东西</h1>
<p class="sub">体系：蒸发-芯吸（柱底 Dirichlet 储液池 C<sub>0</sub> = {case.C0:.2f}，全程液相饱和，全解析孔柱，盒内一个表面孔）·
T = {case.T - 273.15:.0f} °C · RH {case.RH:.0%} · δ<sub>eff</sub> = {case.DELTA_EFF_MM:g} mm [B] · c<sub>k</sub> = {case.C_K:g} ·
<code>{str(any_rec['code_version'])}</code> · 生成于 {time.strftime('%Y-%m-%d %H:%M')}</p>
<p class="key">目的不是与 Fig 23 吻合，而是检查耦合模型<b>自身</b>是否自洽。盐结构是算例输入，两种取法并排：
<b>单尺度</b>（孔径不随高度变，ε<sub>s</sub> = 0.43 [M, Prat 2026 基线]）与<b>分层</b>（每层晶粒与孔径减半，ε<sub>s</sub> = 0.20 [L, Wijnhorst 2025]）。
旧记录用的"分层孔径 + 0.43 孔隙率"混搭不在这里出现。</p>
{pending}

<h2>判据 ①　全局水、盐账闭合</h2>
<div class="tbl"><table><tr><th>算例</th><th>盐结构</th><th>vpp</th><th>孔内盐（全局）</th><th>界面</th><th>网络→结构</th><th>蒸发交盐→落盐</th><th>固相存量 max</th><th>液体 进=出+存</th><th>出孔液体 − 结构所用（扣除无源团簇后；诊断量）</th><th>无源团簇的盐 / 进盐</th><th>界面缓冲 / 结晶</th><th>成核跳变 / 进液</th></tr>{books}</table></div>
<p class="sub">前七列是相对残差，应在舍入量级。后三列是<b>具名的缺口</b>，不是残差：没有孔供液的湿团簇（它的水和盐都没有来源，真实情况会干掉——模型没有干化态）；
结晶了但尚未落到面上的盐（孔–盐体接触处的晶体，无体积 [A]）；成核瞬间把孔重置到 C<sub>sat</sub> 带来的 V·Δρ（无通量）。</p>

<h2>判据 ②　盐的 k、P<sub>c</sub>、ε 出自同一种微结构</h2>
<div class="tbl"><table><tr><th>算例</th><th>盐结构</th><th>vpp</th><th>首层晶粒 d<sub>c</sub> [µm]</th><th>首层孔 a₁ [µm]</th><th>k<sub>s</sub>(第 1 层) [m²]</th><th>k<sub>s</sub>(末层) [m²]</th><th>P<sub>c</sub>(第 1 层) [kPa]</th><th>出生孔隙率（均值）</th><th>末态层号</th></tr>{micro}</table></div>
<p class="sub">k<sub>s</sub> = c<sub>k</sub>·2a³/(12 d<sub>c</sub>)、P<sub>c</sub> = 2γ/a、D<sub>eff</sub> = ε<sup>1.5</sup>(1/D<sub>v</sub>+1/D<sub>K</sub>)<sup>−1</sup> 都读同一个 a(n) 与同一个 ε；c<sub>k</sub> = 1 即不再把渗透率单独压低。
<b>仍未统一的一项（声明）</b>：生长律里的比表面 a<sub>v</sub> = 5(1−ε)/d<sub>c</sub> 按第 1 层晶粒取，分层取法下外层更细的晶粒没有进 Da*——逐层 a<sub>v</sub> 另行立项。</p>

<h2>判据 ③　状态转换由模型变量触发</h2>
<div class="tbl"><table><tr><th>算例</th><th>盐结构</th><th>vpp</th><th>首晶 [s]</th><th>冻结团簇数 max</th><th>封闭体素 max</th><th>柱内黏性 [Pa]</th><th>重力 ρgh [Pa]</th><th>盐体内黏性 [Pa]</th><th>末层 P<sub>c</sub> [Pa]</th><th>需求合计 / P<sub>c</sub></th><th>h<sub>max</sub> [mm]</th><th>末态 E/E₀</th><th>结束</th></tr>{trans}</table></div>
<p class="sub">首晶 = 表面孔到 C<sub>sat</sub>；冻结 = 供液判据失败；封闭 = 反应薄层填满。基线里后两者都应为 0：没有任何转换是"规定"出来的。
需求 = 柱内黏性压降 + ρgh + 盐体内黏性压降，U 取末态储液池供液量 / 开口面积。三项分列：柱内黏性压降在各档都可忽略（≤ 0.5 Pa）。<b>单尺度</b>取法下盐体内黏性项也可忽略，需求就是重力 ρgh，判据失败的高度是毛细上升高度 P<sub>c</sub>/ρg，远在盒高之外。
<b>分层</b>取法下每多一层黏性压降 ×8 而 P<sub>c</sub> 只 ×2：珠径越细、层数越多，盐体内黏性项越大（50–62 µm 档第 6 层已是重力的 10 倍），但到盒顶为止需求 / P<sub>c</sub> 仍在 1e−4——
这就是 Wijnhorst 2025 的结论（毛细–黏性竞争不限制高度，n* ≈ 12）在本模型里的复现，不是新结果。</p>
{_img(fig_dir / 'height.png')}

<h2>判据 ④　闭合是否在适用范围内</h2>
<div class="tbl"><table><tr><th>算例</th><th>盐结构</th><th>vpp</th><th>F(Da*)</th><th>结壳极限 1−ε</th><th>F 的余量</th><th>反应薄层 δ<sub>C</sub> [µm]</th><th>δ<sub>C</sub> / 晶粒</th><th>δ<sub>C</sub> / 体素</th><th>判据 ⑤：末态柱剖面对闭式解的最大偏差</th><th>末态 E/E₀：模型 / 平面极限 δ/(δ−h<sub>max</sub>)</th></tr>{valid}</table></div>
<p class="note"><b>E(t) 的上升不是预测</b>：末列里模型的 E/E₀ 始终略低于把远场面当成固定平面时的几何极限 δ/(δ−h)，即 E(t) 基本只是 h(t)/δ<sub>eff</sub> 的函数——
盐体把扩散路程缩短了。而 δ<sub>eff</sub> 是标定量 [B]，吸收了箱内自由对流；真实的对流边界层不是钉在 δ<sub>eff</sub> 处的一块平面。所以"蒸发缓慢上升"是这个闭合的几何后果，
不能当作与实验（1 mm 档升到 2.2）相比较的预测；DuMuX 柱（δ = 1.6 mm，盒顶 1.5 mm）的 E ≈ 5 是同一效应的极端情形，不应引用。</p>
<p class="note"><b>越界声明</b>：δ<sub>C</sub> / 晶粒 ≪ 1 意味着 Prat 的连续介质反应薄层比一个晶粒还薄——这来自理论本身，不是网格问题；δ<sub>C</sub> / 体素 ≪ 1 说明薄层未被网格分辨，
所以模型用的是它的积分结果 F(Da*)，不解薄层内部。F 离结壳极限的余量在细珠 + 分层取法下只有几个百分点：越界时代码<b>报错而不裁剪</b>。</p>

<h2>判据 ⑥　网格收敛（vpp 7 / 14，同一时刻）</h2>
<div class="tbl"><table><tr><th>算例</th><th>盐结构</th><th>比较时刻 [h]</th><th>首晶 [s] 7 / 14</th><th>E/E₀ 7 / 14</th><th>全程蒸发水量 14 对 7</th><th>全程落盐量 14 对 7</th><th>h<sub>max</sub> [mm] 7 / 14</th><th>已翻转体素体积 14 对 7</th></tr>{conv or '<tr><td colspan="9">尚无 vpp 14 的记录</td></tr>'}</table></div>
<p class="note"><b>只部分满足</b>：积分量（蒸发水量、落盐量、E）在 vpp 7 → 14 之间变化 ≲ 1 %，是收敛的；<b>几何量不是</b>——最高点高 6 % 左右，已翻转体素的体积多 15 % 左右。
同样的盐量在细网格上占了更多的"已翻转"体积，是因为占据度方法里表面那一层未填满的体素不计入：盐体是细柱，表面体素与体内体素同量级，粗网格上"在途"的盐占比大。
所以 vpp 7 报告的高度和体积是<b>下界</b>，对 a 一阶收敛；要报几何量需要外推或更细的网格。</p>

<h2>对照实验（不是目标，只是位置）</h2>
{_img(fig_dir / 'E_vs_experiment.png')}
<p>基线<b>不成壳</b>：任何珠径、两种盐结构下，供液需求都低于毛细阈值（上表末列），盐体持续生长。这与 1 mm 档的实验一致，与细珠档（5–50、50–62 µm 成壳、蒸发塌到 0.04–0.10）<b>不一致</b>。
这是"Prat 生长律 + 单一微结构毛细供液"框架内的<b>无参数陈述</b>：细珠成壳需要这个框架之外的机制（旧链条用 c<sub>k</sub> ~ 1e−5 人为压低渗透率来得到它，见 README §17）。</p>

<h2>各算例的图（vpp = {case.VPP}）</h2>
{figs}
</main></body></html>"""
    Path(output).write_text(html, encoding="utf-8")
    print(f"{output}  ({Path(output).stat().st_size / 1e3:.0f} kB)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="what", required=True)
    r = sub.add_parser("run"); r.add_argument("case", choices=CASES); r.add_argument("structure", choices=list(STRUCTURES))
    r.add_argument("--vpp", type=int, default=None, help=f"default: {case.VPP} (the DuMuX column: its own)")
    r.add_argument("--snap-hours", type=float, default=1.0)
    r.add_argument("--t-end-h", type=float, default=None, help="shorten the run (pilot); default: the case's")
    r.add_argument("--layer-av", action="store_true", help="B': per-layer F, dense birth at the clamp (own record dir)")
    p = sub.add_parser("report"); p.add_argument("--out", type=Path, default=HERE / "report_baseline.html")
    args = ap.parse_args()
    if args.what == "report":
        build_report(args.out)
    else:
        vpp = args.vpp if args.vpp is not None else (dumux1d.RUN_3D_GROWTH["vpp"] if args.case == "dumux1d" else case.VPP)
        run_case(args.case, args.structure, vpp, args.snap_hours, args.t_end_h, layer_av=args.layer_av)


if __name__ == "__main__":
    main()
