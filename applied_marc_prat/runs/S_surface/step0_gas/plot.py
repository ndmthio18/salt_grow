"""Step 0 figures from results/step0_results.json and results/step0_fields.npz (run.py first).

    $PY runs/S_surface/step0_gas/plot.py

figures/fig1_planar_vs_analytical.{png,svg}   V0.1: fully wet plane against the 1-D solution
figures/fig2_pore_domain_height.{png,svg}     V0.3: single pore, 2 mm vs 4 mm numerical domain
figures/fig3_vapour_field.{png,svg}           cross-section of X through the pore, both domains
Every number plotted is read from the result files; nothing is recomputed here.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
FIG = HERE / "figures"
CASES = {"A": "2 mm domain", "B": "4 mm domain"}


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=160)
    fig.savefig(FIG / f"{name}.svg")
    plt.close(fig)


def fig1(r):
    s = r["settings"]; pl = r["plane"]
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))
    # (a) flux: analytical and the two numerical domains
    labels = ["analytical 1-D", "numerical, " + CASES["A"], "numerical, " + CASES["B"]]
    vals = [pl["A"]["j_1d"], pl["A"]["j_projected"], pl["B"]["j_projected"]]
    bars = ax[0].bar(range(3), np.array(vals) * 1e5, color=["0.3", "tab:blue", "tab:orange"], width=0.6)
    for b, v in zip(bars, vals):
        ax[0].text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v:.6e}", ha="center", va="bottom", fontsize=7)
    ax[0].set_xticks(range(3)); ax[0].set_xticklabels(labels, fontsize=8, rotation=12)
    ax[0].set_ylabel(r"evaporation flux $j$ [$10^{-5}$ kg m$^{-2}$ s$^{-1}$]")
    ax[0].set_ylim(0, max(vals) * 1e5 * 1.15)
    ax[0].set_title("(a) fully wet plane: flux per projected area", fontsize=9)
    # (b) relative error of the numerical flux and the conservation residuals
    names = ["rel. error\nvs 1-D", "imbalance\n(E-E_top)/E_top", "discrete\nresidual"]
    x = np.arange(3); w = 0.36
    for k, (case, c) in enumerate((("A", "tab:blue"), ("B", "tab:orange"))):
        e = [abs(pl[case]["rel_err_1d"]), abs(pl[case]["imbalance"]), pl[case]["discrete_residual"]]
        ax[1].bar(x + (k - 0.5) * w, e, w, color=c, label=CASES[case])
        for xi, v in zip(x + (k - 0.5) * w, e):
            ax[1].text(xi, v, f"{v:.1e}", ha="center", va="bottom", fontsize=6.5, rotation=90)
    ax[1].set_yscale("log"); ax[1].set_ylim(1e-16, 1e-8)
    ax[1].axhline(1e-10, color="0.6", ls=":", lw=0.8); ax[1].text(2.45, 1.3e-10, "1e-10", fontsize=7, color="0.4")
    ax[1].set_xticks(x); ax[1].set_xticklabels(names, fontsize=8); ax[1].legend(fontsize=8)
    ax[1].set_title("(b) |relative error| and residuals (log)", fontsize=9)
    # (c) vertical profile against the linear 1-D solution
    for case, c, m in (("A", "tab:blue", "o"), ("B", "tab:orange", "s")):
        X = np.array(pl[case]["X_profile"]); a = pl[case]["a_voxel"]
        z = (np.arange(len(X)) + 0.5) * a
        ax[2].plot(z * 1e3, X, m, color=c, ms=3, mfc="none", label=f"numerical, {CASES[case]}")
    z1 = np.linspace(0, s["delta_eff"], 200)
    ax[2].plot(z1 * 1e3, s["X_wet"] - (s["X_wet"] - s["X_top"]) * z1 / s["delta_eff"], "k-", lw=0.9,
               label=r"1-D: $X_{wet} - (X_{wet}-X_\infty)\,z/\delta_{eff}$")
    ax[2].axvline(s["delta_eff"] * 1e3, color="0.6", ls="--", lw=0.8)
    ax[2].text(s["delta_eff"] * 1e3, s["X_top"], r"  $\delta_{eff}$", fontsize=7, va="bottom")
    ax[2].set_xlabel("height above the plane z [mm]"); ax[2].set_ylabel("vapour mass fraction X [-]")
    ax[2].legend(fontsize=7); ax[2].set_title("(c) profile; far field at z = 20.6 mm", fontsize=9)
    fig.suptitle("Fig. 1  V0.1: the vapour solver on a fully wet plane reproduces the 1-D diffusion solution "
                 "(T = 22 C, RH = 0.07, a_w = 0.753, delta_eff = 20.6 mm; box 2 mm or 4 mm + 1-D segment)",
                 fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); save(fig, "fig1_planar_vs_analytical")


def fig2(r):
    po = r["pore"]; pl = r["plane"]; cmp_ = r["comparison"]
    panels = (("E", "total evaporation rate E [kg/s]", 1.0, "E_rel_diff_AB"),
              ("j_projected", r"E / A$_{box}$ [kg m$^{-2}$ s$^{-1}$]", 1.0, "j_projected_rel_diff_AB"),
              ("j_wet", r"E / A$_{wet}$ [kg m$^{-2}$ s$^{-1}$]", 1.0, "j_wet_rel_diff_AB"))
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))
    for k, (key, lab, scale, dkey) in enumerate(panels):
        vals = [po["A"][key] * scale, po["B"][key] * scale]
        bars = ax[k].bar([0, 1], vals, color=["tab:blue", "tab:orange"], width=0.55)
        for b, v in zip(bars, vals):
            ax[k].text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v:.8e}", ha="center", va="bottom", fontsize=7)
        if key == "j_projected":
            ax[k].axhline(pl["A"]["j_projected"], color="0.3", ls="--", lw=0.9)
            ax[k].text(1.3, pl["A"]["j_projected"], "fully wet\nplane", fontsize=7, va="center")
        ax[k].set_xticks([0, 1]); ax[k].set_xticklabels([CASES["A"], CASES["B"]], fontsize=9)
        ax[k].set_xlim(-0.6, 1.9); ax[k].set_ylim(0, max(vals) * 1.18); ax[k].set_ylabel(lab, fontsize=8.5)
        ax[k].set_title(f"({'abc'[k]})  (A - B)/B = {cmp_[dkey]:+.2e}", fontsize=9)
    q = po["A"]
    fig.suptitle("Fig. 2  V0.3: single wet 540 um square opening (rasterised 500 um) in a 700 um periodic box, "
                 f"physical transfer distance 20.6 mm in both cases; L_spread = {q['L_spread']*1e6:.1f} um, "
                 f"E_pore / E_plane = {cmp_['E_pore_over_E_plane']:.4f}", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, "fig2_pore_domain_height")


def fig3(r, f):
    s = r["settings"]; L = s["L_box"]; a = s["a_voxel"]
    n_lat = f["X_pore_A"].shape[0]; jrow = n_lat // 2        # row whose centre is half a voxel from the pore centre
    y_sec = (jrow + 0.5) * a
    foot_row = {c: f[f"foot_{c}"][:, jrow] for c in CASES}
    sec = {c: f[f"X_pore_{c}"][:, jrow, :] for c in CASES}            # (n_lat, n_z)
    z_zoom = L                                                             # one lateral period
    kz = int(round(z_zoom / a))
    # colour scales: full height shared between A and B; zoom shared; flux magnitude shared
    vmin_full = min(np.nanmin(v) for v in sec.values()); vmax_full = max(np.nanmax(v) for v in sec.values())
    vmin_z = min(np.nanmin(v[:, :kz]) for v in sec.values()); vmax_z = max(np.nanmax(v[:, :kz]) for v in sec.values())
    grads = {}
    for c in CASES:
        X = sec[c][:, :kz]
        gx = -(np.roll(X, -1, 0) - np.roll(X, 1, 0)) / (2 * a)           # periodic in x
        gz = -np.gradient(X, a, axis=1)
        grads[c] = (gx, gz, np.hypot(gx, gz) * s["rho_D"])
    vmax_q = max(gq[2].max() for gq in grads.values())
    fig = plt.figure(figsize=(14.5, 8.0))
    gs = fig.add_gridspec(2, 4, width_ratios=[0.75, 0.75, 1.6, 1.6], wspace=0.55, hspace=0.42,
                          left=0.05, right=0.93, top=0.86, bottom=0.07)
    axf = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[:, 1])]
    axz = [fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 3])]
    axq = [fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[1, 3])]
    xe = np.arange(n_lat + 1) * a * 1e3
    xc = (np.arange(n_lat) + 0.5) * a * 1e3; zc = (np.arange(kz) + 0.5) * a * 1e3
    for k, c in enumerate(CASES):
        n_z = sec[c].shape[1]; ze = np.arange(n_z + 1) * a * 1e3
        im = axf[k].pcolormesh(xe, ze, sec[c].T, vmin=vmin_full, vmax=vmax_full, cmap="viridis", shading="flat")
        axf[k].set_aspect("equal"); axf[k].set_xlabel("x [mm]"); axf[k].set_title(f"{CASES[c]}\nfull height", fontsize=9)
        axf[k].set_ylim(0, 4.0 + 0.0); _mark_opening(axf[k], foot_row[c], a, -0.08)
        axf[k].set_ylabel("z [mm]")
        imz = axz[k].pcolormesh(xe, ze[:kz + 1], sec[c][:, :kz].T, vmin=vmin_z, vmax=vmax_z, cmap="viridis", shading="flat")
        cs = axz[k].contour(xc, zc, sec[c][:, :kz].T, levels=8, colors="w", linewidths=0.5)
        axz[k].clabel(cs, fontsize=5.5, fmt="%.5f")
        axz[k].set_aspect("equal"); axz[k].set_title(f"{CASES[c]}: X, zoom z < {z_zoom*1e3:.1f} mm", fontsize=9)
        axz[k].set_xlabel("x [mm]"); axz[k].set_ylabel("z [mm]")
        _mark_opening(axz[k], foot_row[c], a, -0.03)
        gx, gz, mag = grads[c]
        imq = axq[k].pcolormesh(xe, ze[:kz + 1], mag.T * 1e5, vmin=0.0, vmax=vmax_q * 1e5, cmap="magma", shading="flat")
        axq[k].quiver(xc, zc, gx.T, gz.T, color="w", width=0.004, pivot="mid")
        axq[k].set_aspect("equal"); axq[k].set_xlabel("x [mm]"); axq[k].set_ylabel("z [mm]")
        axq[k].set_title(f"{CASES[c]}: flux magnitude and direction, z < {z_zoom*1e3:.1f} mm", fontsize=9)
        _mark_opening(axq[k], foot_row[c], a, -0.03)
    cb = fig.colorbar(im, ax=axf, fraction=0.08, pad=0.08, shrink=0.6); cb.set_label("X [-] (full-height panels, shared)", fontsize=8)
    cbz = fig.colorbar(imz, ax=axz, fraction=0.025, pad=0.02); cbz.set_label("X [-] (zoom panels, shared)", fontsize=8)
    cbq = fig.colorbar(imq, ax=axq, fraction=0.025, pad=0.02); cbq.set_label(r"$|\rho_g D_v \nabla X|$ [kg m$^{-2}$ s$^{-1}$] (shared)", fontsize=8)
    fig.suptitle(f"Fig. 3  vapour mass fraction X in the section y = {y_sec*1e6:.0f} um through the opening (pore centre y = 350 um, "
                 "half a voxel away); blue bar = wet opening (10 cells = 500 um), hatched = impermeable solid.\n"
                 "Equal aspect in every panel; the far field X_inf = 1.16e-3 lies at z = 20.6 mm, beyond both boxes "
                 "(1-D segment in series).  Left: full numerical heights 2 mm and 4 mm on one colour scale.", fontsize=8.5)
    save(fig, "fig3_vapour_field")


def _mark_opening(ax, foot_row, a, y_off):
    n = len(foot_row); xe = np.arange(n + 1) * a * 1e3
    for i in range(n):
        if foot_row[i]:
            ax.add_patch(plt.Rectangle((xe[i], y_off), xe[i + 1] - xe[i], -y_off, color="tab:blue", lw=0, clip_on=False))
        else:
            ax.add_patch(plt.Rectangle((xe[i], y_off), xe[i + 1] - xe[i], -y_off, facecolor="0.7", hatch="////",
                                       edgecolor="k", lw=0, clip_on=False))
    ax.set_ylim(y_off, ax.get_ylim()[1])


def main():
    r = json.loads((RES / "step0_results.json").read_text())
    f = np.load(RES / "step0_fields.npz")
    fig1(r); fig2(r); fig3(r, f)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
