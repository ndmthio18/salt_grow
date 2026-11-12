"""Step 2A figures from results/step2a_results.json (run.py first).

    $PY runs/S_surface/step2a_prat/plot.py

fig1_growth_analytic_vs_numerical.{png,svg}   h(t) and v(t) of the reference against h0 + v_Prat t
fig2_reaction_zone_partition.{png,svg}        C*(xi) against eq 19; precipitation partition and porosity
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; FIG = HERE / "figures"


def save(fig, name):
    FIG.mkdir(exist_ok=True); fig.savefig(FIG / f"{name}.png", dpi=160); fig.savefig(FIG / f"{name}.svg"); plt.close(fig)


def fig1(r):
    T = r["B_transient"]; A = r["A"]["a_v_exact"]; c = r["case"]
    t = np.array(T["t"]); h = np.array(T["h"]); v = np.array(T["v"]); h_lin = np.array(T["h_linear"])
    tau = A["delta_C"] / A["v"]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(t / 3600, h * 1e6, "-", color="tab:red", lw=1.8, label="reference solver (eqs 1-3, 6, 7, 20; V0)")
    ax[0].plot(t / 3600, h_lin * 1e6, "--", color="k", lw=1.0, label=r"Prat eq 29: $h_0 + v_{Prat}\,t$")
    ax[0].set_xlabel("t [h]"); ax[0].set_ylabel("structure height h [um]")
    ax[0].set_title(f"(a) h(t) from h0 = {c['h0']*1e3:.1f} mm; 24 h growth error {T['rel_err_h_growth_end']:+.1e}", fontsize=9); ax[0].legend(fontsize=8)
    ax[1].plot(t / tau, v / A["v"], "-", color="tab:red", lw=1.5)
    ax[1].axhline(1.0, color="k", lw=0.8, ls="--"); ax[1].axvline(T["t_transient_1e6"] / tau, color="0.5", lw=0.8, ls=":")
    ax[1].set_xscale("log"); ax[1].set_xlabel(r"t / ($\delta_C$ / v)  [-]"); ax[1].set_ylabel(r"v / v$_{Prat}$  [-]")
    ax[1].set_title(f"(b) start-up transient: v within 1e-6 of steady after {T['t_transient_1e6']/tau:.0f} delta_C/v = {T['t_transient_1e6']/60:.0f} min", fontsize=9)
    ax[1].set_ylim(0.9, 1.01)
    err = (h - h_lin) * 1e9
    ax[2].plot(t / 3600, err, "-", color="tab:red", lw=1.5, label=r"$h - (h_0 + v_{Prat} t)$")
    ax[2].axhline(T["h_offset_vs_steady_line"] * 1e9, color="0.5", lw=0.8, ls=":", label=f"start-up offset {T['h_offset_vs_steady_line']*1e9:.1f} nm")
    ax[2].set_xlabel("t [h]"); ax[2].set_ylabel("h - linear [nm]")
    ax[2].set_title(f"(c) deviation; steady v vs eq 29: {T['rel_err_v_end']:+.2e} (n_band 40)", fontsize=9); ax[2].legend(fontsize=8)
    fig.suptitle("Fig. 1  Prat (2026) analytic growth law against the independent 1-D reaction-transport reference (V0 linearisation), "
                 f"d_sub = {c['d_sub']*1e6:.1f} um, eps_new = {c['eps_new']}, j = {c['j_w']:.1e} kg/m2/s, 21 C parameters", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, "fig1_growth_analytic_vs_numerical")


def fig2(r):
    B = r["B"]["40"]; A = B["analytic"]; C = r["C"]; c = r["case"]; closure = r["A"]["a_v_exact"]["eps_final_closure"]
    xi = np.array(B["profile"]["xi"]); Cs = np.array(B["profile"]["C_star"]); eps = np.array(B["profile"]["eps"])
    dC = A["delta_C"]; sel = xi < 8 * dC
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw=dict(width_ratios=[1.3, 0.9, 1.1]))
    ax[0].plot(xi[sel] / dC, Cs[sel] / A["C_top_star"], "o", color="tab:red", ms=3.5, mfc="none", label="reference solver, cell centres")
    xx = np.linspace(0, 8, 300)
    ax[0].plot(xx, np.exp(-xx), "k-", lw=1, label=r"eq 19 (Da >> 1): $C^*_{top}\,e^{-\xi/\delta_C}$")
    ax[0].set_xlabel(r"depth below the interface $\xi / \delta_C$  ($\delta_C$ = %.3f um)" % (dC * 1e6)); ax[0].set_ylabel(r"$C^* / C^*_{top}$  [-]")
    ax[0].set_title(f"(a) supersaturation profile\nC*_top num {B['C_top_star']:.4e} vs eq 31 {A['C_top_star']:.4e} ({B['rel_C_top']:+.1e}); max |num - eq 19| / C*_top = {B['profile_max_dev_over_Ctop']:.1e}", fontsize=8.5)
    ax[0].legend(fontsize=8)
    num = [C["salt_front_over_total"], C["salt_internal_over_total"]]; ana = [A["F"], 1 - A["F"]]
    x = np.arange(2)
    ax[1].bar(x - 0.18, ana, 0.34, color="0.6", label="Prat: F, 1 - F (eq 29)")
    ax[1].bar(x + 0.18, num, 0.34, color="tab:red", label="reference solver")
    for xi_, a_, n_ in zip(x, ana, num):
        ax[1].text(xi_ - 0.18, a_, f"{a_:.5f}", ha="center", va="bottom", fontsize=7.5); ax[1].text(xi_ + 0.18, n_, f"{n_:.5f}", ha="center", va="bottom", fontsize=7.5)
    ax[1].set_xticks(x); ax[1].set_xticklabels(["front P_i", "internal S"]); ax[1].set_ylabel("share of precipitated salt [-]"); ax[1].set_ylim(0, 1.1)
    ax[1].legend(fontsize=8, loc="upper right"); ax[1].set_title(f"(b) precipitation partition\nfront share vs F: {C['partition_vs_F']:+.1e}", fontsize=8.5)
    ax[2].plot(xi[sel] / dC, eps[sel], "-", color="tab:red", lw=1.5, label="eps(xi), reference solver")
    ax[2].axhline(c["eps_new"], color="0.4", lw=0.8, ls=":", label=f"birth porosity eps_new = {c['eps_new']}")
    ax[2].axhline(closure, color="k", lw=0.8, ls="--", label=f"1 - (1 - eps_new) / F = {closure:.4f}")
    ax[2].set_xlabel(r"$\xi / \delta_C$"); ax[2].set_ylabel("porosity [-]"); ax[2].legend(fontsize=8)
    ax[2].set_title(f"(c) porosity through the reaction zone\neps_final num {C['eps_final']:.4f} vs closure {closure:.4f} ({C['eps_final_vs_closure_analytic']:+.1e})", fontsize=8.5)
    fig.suptitle("Fig. 2  Reaction zone of the V0 reference (n_band 40): the exponential supersaturation layer, the front / internal precipitation partition "
                 "and the final porosity it sets\n(prescribed flux, saturated supply, single scale; the material enters at eps_new and leaves the zone at eps_final)", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90)); save(fig, "fig2_reaction_zone_partition")


def main():
    r = json.loads((RES / "step2a_results.json").read_text()); fig1(r); fig2(r); print("figures written to", FIG)


if __name__ == "__main__":
    main()
