"""Step 1 figures from results/step1_results.json and results/step1_fields.npz (run.py first).

    $PY runs/S_surface/step1_frozen_salt/plot.py

fig1_geometry.{png,svg}        the prescribed pillar and its wet / inert faces (A1, A2)
fig2_E_vs_height.{png,svg}     E / E0 against h; E_top, E_side; mean flux per wet area
fig3_fields_h0.3.{png,svg}     X and local flux around the pillar at h = 0.3 mm, A1 vs A2
Flux arrows: -grad X by central differences of the saved cell-centre section (periodic in x,
one-sided at the pillar and the plane); magnitudes from the solver's face fluxes where shown.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

HERE = Path(__file__).resolve().parent
RES = HERE / "results"; FIG = HERE / "figures"


def save(fig, name):
    FIG.mkdir(exist_ok=True); fig.savefig(FIG / f"{name}.png", dpi=160); fig.savefig(FIG / f"{name}.svg"); plt.close(fig)


def fig1(r):
    s = r["settings"]; L = s["L_box"] * 1e3; w = s["pore_side"] * 1e3; h = 0.3; H = s["h_box"] * 1e3
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, (cond, title) in enumerate((("A1", "A1: wet top, inert sides"), ("A2", "A2: wet top and wet sides"))):
        a = ax[k]; x0 = (L - w) / 2
        a.add_patch(Rectangle((0, -0.12), L, 0.12, facecolor="0.75", hatch="////", edgecolor="k", lw=0.8))   # substrate
        a.add_patch(Rectangle((x0, 0), w, h, facecolor="0.9", edgecolor="k", lw=0.8))                          # pillar
        a.plot([x0, x0 + w], [h, h], color="tab:blue", lw=4, solid_capstyle="butt")                            # wet top
        side_c = "tab:blue" if cond == "A2" else "k"; side_lw = 4 if cond == "A2" else 1.5; ls = "-" if cond == "A2" else "--"
        for xs in (x0, x0 + w):
            a.plot([xs, xs], [0, h], color=side_c, lw=side_lw, ls=ls)
        a.plot([x0, x0 + w], [0.0, 0.0], color="k", lw=1.5, ls=":")                                            # covered opening
        a.plot([0, L], [H, H], color="0.4", lw=0.8, ls="-.")
        a.annotate("", xy=(x0 + w + 0.03, 0), xytext=(x0 + w + 0.03, h), arrowprops=dict(arrowstyle="<->", lw=0.8))
        a.text(x0 + w + 0.05, h / 2, f"h = {h:.1f} mm\n(0.1 / 0.3 / 0.5)", fontsize=8, va="center")
        a.annotate("", xy=(x0, -0.2), xytext=(x0 + w, -0.2), arrowprops=dict(arrowstyle="<->", lw=0.8))
        a.text(L / 2, -0.26, f"{w*1e3:.0f} um opening = pillar footprint\nA_top = {w*w:.4f} mm2", ha="center", va="top", fontsize=8)
        a.annotate("", xy=(0, 0.62), xytext=(L, 0.62), arrowprops=dict(arrowstyle="<->", lw=0.8))
        a.text(L / 2, 0.64, f"periodic cell {L*1e3:.0f} um (images on both sides)", ha="center", fontsize=8)
        a.text(0.02, H - 0.05, f"box top z = {H:.0f} mm -> 1-D segment to the far field at z = 20.6 mm (X = RH W_sat)", fontsize=7.5, va="top")
        a.text(L / 2, 0.42, "gas: div(rho_g D_v grad X) = 0", ha="center", fontsize=8, color="0.3")
        if cond == "A2":
            a.text(x0 - 0.02, h / 2, "wet\nside", ha="right", va="center", fontsize=7.5, color="tab:blue")
            a.text(x0 + w / 2, h / 2, f"A_side = 4 x {w:.2f} x h\n= {4*w*h:.3f} mm2 at h = 0.3", ha="center", va="center", fontsize=7.5)
        else:
            a.text(x0 - 0.02, h / 2, "inert\nside", ha="right", va="center", fontsize=7.5)
            a.text(x0 + w / 2, h / 2, "prescribed solid\n(no transport inside)", ha="center", va="center", fontsize=7.5)
        a.text(x0 + w / 2, h + 0.03, "wet top: X = a_w(C_sat) W_sat", ha="center", fontsize=7.5, color="tab:blue")
        a.text(x0 + w / 2, -0.06, "covered opening: no gas face", ha="center", fontsize=7, color="k")
        a.set_xlim(-0.05, L + 0.05); a.set_ylim(-0.4, 0.75); a.set_aspect("equal")
        a.set_xlabel("x [mm]"); a.set_ylabel("z [mm]"); a.set_title(title, fontsize=10)
    fig.suptitle("Fig. 1  Prescribed wet external geometry over the single surface pore (vertical section through the pore centre;\n"
                 "A0 = the bare opening, h = 0).  Not a model of the salt: faces are prescribed wet or inert.", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); save(fig, "fig1_geometry")


def fig2(r):
    p = r["primary"]; rat = r["ratios"]; hs = np.array(r["settings"]["heights"]) * 1e3
    E0 = p["A0"]["E_total"]
    keys = [f"h{h:.1f}" for h in hs]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    ax[0].plot(hs, [rat[k]["R_height"] for k in keys], "o-", color="tab:blue", label="A1: wet top, inert sides  (R_height)")
    ax[0].plot(hs, [rat[k]["R_total"] for k in keys], "s-", color="tab:red", label="A2: wet top + sides  (R_total)")
    ax[0].plot(hs, [rat[k]["control_1d"] for k in keys], ":", color="0.4", label=r"1-D planar control $\delta_{eff}/(\delta_{eff}-h)$")
    ax[0].axhline(1.0, color="0.7", lw=0.7)
    for k, h in zip(keys, hs):
        ax[0].text(h, rat[k]["R_total"], f" {rat[k]['R_total']:.4f}", fontsize=7, va="bottom", color="tab:red")
        ax[0].text(h, rat[k]["R_height"], f" {rat[k]['R_height']:.4f}", fontsize=7, va="top", color="tab:blue")
    ax[0].set_xlabel("prescribed pillar height h [mm]"); ax[0].set_ylabel(r"$E_{total} / E_0$  [-]")
    ax[0].set_title(f"(a) total evaporation relative to the bare opening\nE0 = {E0:.4e} kg/s", fontsize=9); ax[0].legend(fontsize=7)
    ax[1].plot(hs, [p[f"A1_{k}"]["E_top"] for k in keys], "o-", color="tab:blue", label="A1 E_top (= E_total)")
    ax[1].plot(hs, [p[f"A2_{k}"]["E_top"] for k in keys], "s-", color="tab:red", label="A2 E_top")
    ax[1].plot(hs, [p[f"A2_{k}"]["E_side"] for k in keys], "s--", color="tab:orange", label="A2 E_side")
    ax[1].plot(hs, [p[f"A2_{k}"]["E_total"] for k in keys], "s:", color="tab:red", mfc="none", label="A2 E_total")
    for k, h in zip(keys, hs):
        ax[1].text(h, p[f"A2_{k}"]["E_side"], f" {rat[k]['side_share']*100:.1f} %", fontsize=7, va="bottom", color="tab:orange")
    ax[1].set_xlabel("h [mm]"); ax[1].set_ylabel("E [kg/s]"); ax[1].set_title("(b) top and side contributions (A2 side share annotated)", fontsize=9)
    ax[1].legend(fontsize=7)
    ax[2].plot(hs, [p[f"A1_{k}"]["j_top_mean"] for k in keys], "o-", color="tab:blue", label="A1 j_top_mean = E/A_top")
    ax[2].plot(hs, [p[f"A2_{k}"]["j_top_mean"] for k in keys], "s-", color="tab:red", label="A2 j_top_mean")
    ax[2].plot(hs[1:], [p[f"A2_{k}"]["j_side_mean"] for k in keys[1:]], "s--", color="tab:orange", label="A2 j_side_mean = E_side/A_side")
    ax[2].plot(hs, [p[f"A2_{k}"]["j_wet_mean"] for k in keys], "s:", color="tab:red", mfc="none", label="A2 j_wet_mean = E/A_total_wet")
    ax[2].plot(hs, [p[f"A2_{k}"]["j_projected"] for k in keys], "-", color="0.4", lw=0.8, label="A2 j_projected = E/A_box")
    ax[2].set_xlabel("h [mm]"); ax[2].set_ylabel(r"mean flux [kg m$^{-2}$ s$^{-1}$]"); ax[2].set_title("(c) mean flux per wet area", fontsize=9)
    ax[2].legend(fontsize=7)
    fig.suptitle("Fig. 2  Prescribed pillar height and wet side faces against the bare opening (20 um voxels, 2 mm box, far field 20.6 mm).\n"
                 "E_total is set by the far-field resistance; the wet area enters only through the near-field spreading resistance.", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92)); save(fig, "fig2_E_vs_height")


def fig3(r, f):
    s = r["settings"]; a = s["a"]; p = r["primary"]
    names = ("A1_h0.3", "A2_h0.3"); titles = ("A1: wet top, inert sides", "A2: wet top + wet sides")
    n_lat = f["X_A1_h0.3"].shape[0]; kz = int(round(0.9e-3 / a))
    xe = np.arange(n_lat + 1) * a * 1e3; ze = np.arange(kz + 1) * a * 1e3
    xc = (np.arange(n_lat) + 0.5) * a * 1e3; zc = (np.arange(kz) + 0.5) * a * 1e3
    secs = {n: f[f"X_{n}"][:, :kz] for n in names}
    vmin = min(np.nanmin(v) for v in secs.values()); vmax = max(np.nanmax(v) for v in secs.values())
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.4), gridspec_kw=dict(width_ratios=[1, 1, 1.15]), constrained_layout=True)
    for k, n in enumerate(names):
        X = secs[n]; pil = f[f"pillar_{n}"][:, :kz]
        im = ax[0, k].pcolormesh(xe, ze, X.T, vmin=vmin, vmax=vmax, cmap="viridis", shading="flat")
        cs = ax[0, k].contour(xc, zc, np.where(pil, np.nan, X).T, levels=10, colors="w", linewidths=0.5)
        ax[0, k].clabel(cs, fontsize=5.5, fmt="%.5f")
        _draw_pillar(ax[0, k], pil, a, wet_sides=(n.startswith("A2")))
        ax[0, k].set_aspect("equal"); ax[0, k].set_title(f"{titles[k]}: X (section y = 350 um)", fontsize=9)
        ax[0, k].set_xlabel("x [mm]"); ax[0, k].set_ylabel("z [mm]")
        # flux magnitude and direction from the section
        Xn = np.where(pil, np.nan, X)
        gx = -(np.roll(Xn, -1, 0) - np.roll(Xn, 1, 0)) / (2 * a)
        gz = -np.gradient(Xn, a, axis=1)
        mag = np.hypot(gx, gz) * s["delta_eff"] * 0 + np.hypot(gx, gz) * r["primary"]["A0"]["rho_D"]
        imq = ax[1, k].pcolormesh(xe, ze, mag.T * 1e5, vmin=0, vmax=np.nanmax([np.hypot(*_grad(secs[m], f[f"pillar_{m}"][:, :kz], a)) for m in names]) * r["primary"]["A0"]["rho_D"] * 1e5,
                                  cmap="magma", shading="flat")
        st = 2
        ax[1, k].quiver(xc[::st], zc[::st], gx.T[::st, ::st], gz.T[::st, ::st], color="w", width=0.004, pivot="mid")
        _draw_pillar(ax[1, k], pil, a, wet_sides=(n.startswith("A2")))
        ax[1, k].set_aspect("equal"); ax[1, k].set_title(f"{titles[k]}: |rho_g D_v grad X| and direction", fontsize=9)
        ax[1, k].set_xlabel("x [mm]"); ax[1, k].set_ylabel("z [mm]")
    fig.colorbar(im, ax=ax[0, :2].tolist(), fraction=0.03, pad=0.02, shrink=0.9).set_label("X [-] (shared)", fontsize=8)
    fig.colorbar(imq, ax=ax[1, :2].tolist(), fraction=0.03, pad=0.02, shrink=0.9).set_label(r"$|\rho_g D_v \nabla X|$ [$10^{-5}$ kg m$^{-2}$ s$^{-1}$] (shared)", fontsize=8)
    # local flux on the faces (solver face fluxes)
    jc = n_lat // 2
    for n, c, lab in ((names[0], "tab:blue", "A1"), (names[1], "tab:red", "A2")):
        jt = f[f"jtop_{n}"][:, jc]
        ax[0, 2].plot(xc, jt * 1e5, "-", color=c, label=f"{lab} top face, row y = 350 um (mean {p[n]['j_top_mean']*1e5:.3f})")
    ax[0, 2].axhline(p["A0"]["j_top_mean"] * 1e5, color="0.4", ls=":", label=f"A0 bare opening mean {p['A0']['j_top_mean']*1e5:.3f}")
    ax[0, 2].set_xlabel("x [mm]"); ax[0, 2].set_ylabel(r"top-face flux [$10^{-5}$ kg m$^{-2}$ s$^{-1}$]"); ax[0, 2].legend(fontsize=7)
    ax[0, 2].set_title("local flux on the wet top face (solver face flux)", fontsize=9)
    js = f["jside_A2_h0.3_0_+1"]           # gas cells at x < pillar looking +x into the pillar: the pillar's -x side
    col = np.where(np.isfinite(js).any(axis=1))[0]
    for ci in col:
        zz = np.where(np.isfinite(js[ci]))[0]
        ax[1, 2].plot(js[ci, zz] * 1e5, (zz + 0.5) * a * 1e3, "s-", color="tab:orange", label=f"A2 side face x = {(ci+1)*a*1e3:.2f} mm (mean {p['A2_h0.3']['j_side_mean']*1e5:.3f})")
    ax[1, 2].axvline(p["A2_h0.3"]["j_top_mean"] * 1e5, color="tab:red", ls=":", label="A2 top mean")
    ax[1, 2].set_xlabel(r"side-face flux [$10^{-5}$ kg m$^{-2}$ s$^{-1}$]"); ax[1, 2].set_ylabel("z on the side face [mm]"); ax[1, 2].legend(fontsize=7)
    ax[1, 2].set_title("local flux on a wet side face, A2 (section row)", fontsize=9)
    fig.suptitle("Fig. 3  Vapour field and local flux around the prescribed pillar, h = 0.3 mm, 20 um voxels, section through the pore centre "
                 "(y = 350 um);\ngrey = prescribed solid, blue = wet face, dashed black = inert face, hatched = substrate.", fontsize=8.5)
    save(fig, "fig3_fields_h0.3")


def _grad(X, pil, a):
    Xn = np.where(pil, np.nan, X)
    return -(np.roll(Xn, -1, 0) - np.roll(Xn, 1, 0)) / (2 * a), -np.gradient(Xn, a, axis=1)


def _draw_pillar(ax, pil, a, wet_sides):
    n_lat, kz = pil.shape
    cols = np.where(pil.any(axis=1))[0]; n_h = int(pil.sum(axis=1).max())
    if n_h == 0:
        return
    x0, x1 = cols[0] * a * 1e3, (cols[-1] + 1) * a * 1e3; h = n_h * a * 1e3
    ax.add_patch(Rectangle((x0, 0), x1 - x0, h, facecolor="0.85", edgecolor="none", zorder=3))
    ax.plot([x0, x1], [h, h], color="tab:blue", lw=3, zorder=4)
    for xs in (x0, x1):
        ax.plot([xs, xs], [0, h], color="tab:blue" if wet_sides else "k", lw=3 if wet_sides else 1.2, ls="-" if wet_sides else "--", zorder=4)
    ax.add_patch(Rectangle((0, -0.04), n_lat * a * 1e3, 0.04, facecolor="0.75", hatch="////", edgecolor="k", lw=0.5, clip_on=False, zorder=3))
    ax.set_ylim(-0.04, kz * a * 1e3)


def main():
    r = json.loads((RES / "step1_results.json").read_text()); f = np.load(RES / "step1_fields.npz")
    fig1(r); fig2(r); fig3(r, f); print("figures written to", FIG)


if __name__ == "__main__":
    main()
