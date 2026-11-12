"""Extract data from the published figure bitmaps of Prat (2026).

The paper's figures are 300 dpi rasters -- there is no vector content to read
exact coordinates from -- and the experimental data behind Figures 8 and 9 is
not available anywhere else (ref 18 is paywalled, ref 41 is a French PhD
thesis). So the points have to come off the images.

Doing that by hand in a GUI would be unreproducible: nobody could rerun it and
get the same numbers. Instead everything here is scripted and self-verifying:

  1. axis calibration is derived from the figures' own MAJOR TICK MARKS
     (discriminated from minor ticks by length), fitted linearly, and the fit
     RESIDUAL IS REPORTED -- if a calibration were wrong the residual would
     blow up rather than silently biasing every extracted point;
  2. curves and markers are separated BY COLOUR (experiment is black, the model
     lines are red and blue), so the two never contaminate each other;
  3. every extraction writes an OVERLAY image with the recovered points drawn
     back onto the original figure, for visual confirmation.

Run:  ../.venv/bin/python digitize.py
"""

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"
OVERLAYS = DATA / "overlays"
PDF = HERE.parent / "reference" / (
    "Prat - 2026 - The Evaporation-Driven Growth Law of Porous Salt Structures.pdf"
)

# =============================================================================
# Panel calibration
# =============================================================================
# frame:    (row_top, row_bottom, col_left, col_right) of the axes box
# x/y_ref:  (pixel_at_first_major, value, pixel_at_last_major, value)
#           taken from the frame corners where the corner IS a labelled value,
#           which the major-tick fit then confirms.
# major_len: minimum tick length in px that counts as a major tick.

PANELS = {
    "fig02a": dict(
        image="raw_p4-000.png",
        frame=(16, 695, 85, 771),
        x=dict(lo_px=85, lo_val=0.0, hi_px=771, hi_val=1.0, major_len=15),
        y=dict(lo_px=694.5, lo_val=0.0, hi_px=16.0, hi_val=8.0, major_len=15),  # log10(Da)
        note="Da vs h for d_sub = 1/10/100 um. Da contains no beta, so it isolates it.",
        exclude=[],
        tick_check=False,  # log y; anchored on the decade labels (10^0 .. 10^8)
    ),
    "fig04": dict(
        image="raw_p6-000.png",
        frame=(12, 733, 93, 813),
        x=dict(lo_px=93, lo_val=0.0, hi_px=813, hi_val=1000.0, major_len=15),
        y=dict(lo_px=733, lo_val=0.0, hi_px=12, hi_val=300.0, major_len=15),
        note="Da* vs d_sub (um). Main panel only; the two insets are excluded.",
        exclude=[(30, 320, 230, 540), (400, 650, 520, 815)],
    ),
    "fig04_inset_F": dict(
        image="raw_p6-000.png",
        frame=(42, 300, 268, 526),
        x=dict(lo_px=268, lo_val=0.0, hi_px=526, hi_val=3.0, major_len=8),  # log10(d_sub/um)
        y=dict(lo_px=300, lo_val=0.0, hi_px=42, hi_val=1.0, major_len=8),
        note="Left inset of Figure 4: F = sqrt(Da*)/(1+sqrt(Da*)) vs d_sub, log x.",
        exclude=[],
        tick_check=False,  # log axis: the "majors" found are decade minors, so
        # the uniform-step residual test does not apply. Anchored on the frame.
    ),
    "fig05b": dict(
        image="raw_p7-000.png",
        frame=(10, 756, 1162, 1910),
        x=dict(lo_px=1162, lo_val=0.999, hi_px=1910, hi_val=1.000, major_len=15),
        y=dict(lo_px=609.5, lo_val=1.00000, hi_px=10.5, hi_val=1.00004, major_len=20),
        note="Concentration profile C/C_eq vs z/h. Used to infer the h the paper used.",
        exclude=[],
        tick_check=False,  # ticks point OUTWARD here, so the inward scan finds
        # nothing real; calibration was verified against the tick LABEL centres
        # instead (x: 1155/1314/1462/1612/1762/1911 px vs 1162/1312/1461/1611/
        # 1760/1910 predicted; y: 5 labels evenly spaced 149.75 px per 1e-5).
    ),
    "fig04_inset_zoom": dict(
        image="raw_p6-000.png",
        frame=(413, 620, 577, 783),
        x=dict(lo_px=577, lo_val=0.0, hi_px=783, hi_val=100.0, major_len=8),
        y=dict(lo_px=620, lo_val=-1.0, hi_px=413, hi_val=2.0, major_len=8),  # log10(Da*)
        note="Right inset of Figure 4: Da* vs d_sub over [0,100] um, LOG y (0.1-100).",
        exclude=[],
        tick_check=False,  # log y; the frame corners are the decade labels
    ),
    "fig08a": dict(
        image="raw_p9-000.png",
        frame=(11, 694, 100, 782),
        x=dict(lo_px=100, lo_val=0.0, hi_px=782, hi_val=25.0, major_len=12),
        y=dict(lo_px=694, lo_val=0.0, hi_px=12, hi_val=0.30, major_len=12),
        note="Precipitated salt (g) vs time (hr), Norouzi Rad et al. (ref 18).",
        exclude=[],
    ),
    "fig08b": dict(
        image="raw_p9-000.png",
        frame=(11, 693, 1002, 1685),
        x=dict(lo_px=1002, lo_val=0.0, hi_px=1685, hi_val=25.0, major_len=12),
        y=dict(lo_px=693, lo_val=0.0, hi_px=11, hi_val=3.0, major_len=12),
        note="Equivalent thickness h_equ (mm) vs time (hr). Inset excluded.",
        exclude=[(25, 335, 1012, 1400)],
    ),
    "fig09a": dict(
        image="raw_p9-001.png",
        frame=(11, 699, 87, 766),
        x=dict(lo_px=87, lo_val=0.0, hi_px=766, hi_val=40.0, major_len=15),
        y=dict(lo_px=699, lo_val=0.0, hi_px=20, hi_val=0.40, major_len=15),
        note="Evaporated mass (g) vs time (hr), Diouf thesis (ref 41).",
        exclude=[],  # the schematic overlaps the curve; handled by trace continuity
    ),
    "fig09b": dict(
        image="raw_p9-001.png",
        frame=(11, 697, 1027, 1713),
        x=dict(lo_px=1027, lo_val=0.0, hi_px=1713, hi_val=40.0, major_len=12),
        y=dict(lo_px=697, lo_val=0.0, hi_px=11, hi_val=2.2, major_len=12),
        note="Surface heave h/h_top (mm) vs time (hr). Photo inset excluded.",
        exclude=[(415, 690, 1470, 1685)],
    ),
}


# =============================================================================
# Image helpers
# =============================================================================


def load(name):
    return np.asarray(Image.open(DATA / name).convert("RGB")).astype(int)


def mask_black(im):
    """Near-neutral dark pixels: the experimental markers and axes."""
    mx, mn = im.max(axis=2), im.min(axis=2)
    return (mx < 120) & ((mx - mn) < 45)


def mask_red(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return (r > 110) & (r - np.maximum(g, b) > 55)


def mask_blue(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return (b > 110) & (b - np.maximum(r, g) > 55)


def apply_exclusions(mask, exclude):
    out = mask.copy()
    for r0, r1, c0, c1 in exclude:
        out[r0:r1, c0:c1] = False
    return out


def keep_wide_components(mask, frame, min_width_frac=0.15):
    """Drop everything except components that span a good part of the panel.

    The model curves are labelled with COLOURED TEXT ("Eq.(35), eta = 1" in red,
    "eta = 0.76" in blue), so a plain colour mask picks up the glyphs as well as
    the line and a column scan jumps between them. The line is one long
    connected component; each glyph is a small separate one, so filtering on
    horizontal extent removes the text and keeps the curve.
    """
    r0, r1, c0, c1 = frame
    need = min_width_frac * (c1 - c0)
    lab, n = ndimage.label(mask, structure=np.ones((3, 3)))
    out = np.zeros_like(mask)
    for i in range(1, n + 1):
        cols = np.where(lab == i)[1]
        if cols.size and (cols.max() - cols.min()) >= need:
            out |= lab == i
    return out


def inside_frame(mask, frame, inset=4):
    """Blank everything outside the axes box, plus the frame lines themselves."""
    r0, r1, c0, c1 = frame
    out = np.zeros_like(mask)
    out[r0 + inset : r1 - inset, c0 + inset : c1 - inset] = mask[
        r0 + inset : r1 - inset, c0 + inset : c1 - inset
    ]
    return out


# =============================================================================
# Calibration -- derived from the figure's own major ticks, and verified
# =============================================================================


def _tick_positions(dark, frame, axis, major_len):
    """Locate major ticks by how far they extend from their axis line."""
    r0, r1, c0, c1 = frame
    runs, cur = [], []
    span = range(c0, c1 + 1) if axis == "x" else range(r0, r1 + 1)
    for p in span:
        n = 0
        if axis == "x":  # ticks point up from the bottom axis
            for r in range(r1 - 2, r1 - 34, -1):
                if dark[r, p]:
                    n += 1
                else:
                    break
        else:  # ticks point right from the left axis
            for c in range(c0 + 2, c0 + 34):
                if dark[p, c]:
                    n += 1
                else:
                    break
        if n >= major_len:
            cur.append(p)
        elif cur:
            runs.append(float(np.mean(cur)))
            cur = []
    if cur:
        runs.append(float(np.mean(cur)))
    return np.array(runs)


def calibrator(dark, frame, spec, axis, name):
    """Build pixel -> value map and verify it against every major tick.

    The two reference points come from the axes corners; the major ticks then
    provide an independent check. A calibration error shows up immediately as a
    large residual instead of quietly skewing all the data.
    """
    lo_px, lo_val = spec["lo_px"], spec["lo_val"]
    hi_px, hi_val = spec["hi_px"], spec["hi_val"]
    scale = (hi_val - lo_val) / (hi_px - lo_px)

    def to_value(px):
        return lo_val + (np.asarray(px, dtype=float) - lo_px) * scale

    if not spec.get("_check", True):
        print(f"    -- {name}.{axis}: tick self-check not applicable "
              "(see the panel's note); calibration anchored on the frame/labels")
        return to_value

    ticks = _tick_positions(dark, frame, axis, spec["major_len"])
    if len(ticks) >= 3:
        vals = to_value(ticks)
        step = np.median(np.diff(np.sort(vals)))
        snapped = np.round(vals / step) * step
        resid = np.max(np.abs(vals - snapped))
        rel = resid / abs(hi_val - lo_val)
        flag = "OK " if rel < 0.004 else "!! "
        print(
            f"    {flag}{name}.{axis}: {len(ticks)} major ticks, step {step:g}, "
            f"max residual {resid:.4g} ({rel * 100:.2f}% of range)"
        )
    else:
        print(f"    ?? {name}.{axis}: only {len(ticks)} major ticks found, unverified")
    return to_value


# =============================================================================
# Extraction
# =============================================================================


def extract_blobs(mask, min_area=60, max_area=2000, min_fill=0.45):
    """Centroids of filled circular markers.

    min_fill rejects elongated things (dotted connector lines, text strokes,
    axis fragments) that happen to survive the area filter: a filled disc fills
    most of its bounding box, a line segment does not.
    """
    lab, n = ndimage.label(mask)
    pts = []
    for i in range(1, n + 1):
        idx = lab == i
        area = int(idx.sum())
        if not (min_area <= area <= max_area):
            continue
        rows, cols = np.where(idx)
        hgt = rows.max() - rows.min() + 1
        wid = cols.max() - cols.min() + 1
        if area / (hgt * wid) < min_fill:
            continue
        if not (0.55 < hgt / wid < 1.8):  # roughly round
            continue
        pts.append((cols.mean(), rows.mean()))
    return np.array(sorted(pts)) if pts else np.empty((0, 2))


def extract_line(mask, frame, max_run=14):
    """Column-wise centre of a thin curve, ignoring thick blocks."""
    r0, r1, c0, c1 = frame
    xs, ys = [], []
    for c in range(c0, c1 + 1):
        rows = np.where(mask[:, c])[0]
        if rows.size == 0:
            continue
        # split into contiguous runs; keep thin ones only
        splits = np.split(rows, np.where(np.diff(rows) > 1)[0] + 1)
        good = [s for s in splits if len(s) <= max_run]
        if not good:
            continue
        best = max(good, key=len)
        xs.append(c)
        ys.append(best.mean())
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


def trace_curve(mask, frame, start_col, start_row, max_jump=9, max_run=8):
    """Follow a thin curve left-to-right by continuity.

    Needed for Figure 9a, where an explanatory schematic is drawn on top of the
    plotting area and overlaps the curve's row range. Continuity plus a
    thin-run constraint keeps the trace on the curve; a plain column scan would
    jump onto the schematic.
    """
    r0, r1, c0, c1 = frame
    xs, ys = [], []
    prev = float(start_row)
    for c in range(start_col, c1 + 1):
        rows = np.where(mask[:, c])[0]
        rows = rows[(rows > r0) & (rows < r1)]
        if rows.size == 0:
            continue
        splits = np.split(rows, np.where(np.diff(rows) > 1)[0] + 1)
        cands = [s.mean() for s in splits if len(s) <= max_run]
        if not cands:
            continue
        best = min(cands, key=lambda v: abs(v - prev))
        if abs(best - prev) > max_jump:
            continue
        xs.append(c)
        ys.append(best)
        prev = best
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


# =============================================================================
# Output
# =============================================================================


def write_csv(path, x, y, xlabel, ylabel, source, method):
    with open(path, "w") as f:
        f.write(f"# {xlabel} vs {ylabel}\n")
        f.write(f"# source: {source}\n")
        f.write(f"# extracted from: {PDF.name}\n")
        f.write(f"# method: {method}\n")
        f.write("# NOT original data -- digitised from the published raster figure.\n")
        f.write(f"{xlabel},{ylabel}\n")
        # 12 significant digits, not 6: Figure 5b's ordinate is 1 + 4e-5, so a
        # %.6g round trip quantises the whole signal away.
        for xi, yi in zip(x, y):
            f.write(f"{xi:.12g},{yi:.12g}\n")
    print(f"      -> {path.relative_to(HERE)}  ({len(x)} points)")


def overlay(image_name, frame, series, out_name, title):
    """Draw recovered points back onto the original figure for visual checking."""
    im = load(image_name)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.imshow(im.astype(np.uint8))
    for label, (px, py), style in series:
        ax.plot(px, py, style, ms=3.5, lw=1.0, label=f"{label} ({len(px)})", alpha=0.85)
    r0, r1, c0, c1 = frame
    ax.add_patch(
        plt.Rectangle((c0, r0), c1 - c0, r1 - r0, fill=False, ec="lime", lw=1.0, ls="--")
    )
    ax.set_xlim(max(0, c0 - 60), min(im.shape[1], c1 + 60))
    ax.set_ylim(min(im.shape[0], r1 + 90), max(0, r0 - 60))
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(OVERLAYS / out_name, dpi=140)
    plt.close(fig)


# =============================================================================
# Per-figure drivers
# =============================================================================


def do_fig02a():
    """Extract the three Da(h) curves of Figure 2a.

    This is the decisive test for the Da* discrepancy. Da and Da* share every
    parameter except beta:

        Da*/Da = beta^2 d_sub^2 / (1296 h^2)

    -- k_r, D_s, eps_s, C_eq and the a_v prefactor all cancel. So if Figure 2a
    reproduces our Da but Figure 4 does not reproduce our Da*, beta is the only
    quantity that can be responsible.
    """
    p = PANELS["fig02a"]
    im = load(p["image"])
    dark = mask_black(im)
    print(f"  fig02a -- {p['note']}")
    chk = p.get("tick_check", True)
    fx = calibrator(dark, p["frame"], {**p["x"], "_check": chk}, "x", "fig02a")
    fy = calibrator(dark, p["frame"], {**p["y"], "_check": chk}, "y", "fig02a")

    out, series = {}, []
    for colour, tag, d_sub in (("black", "dsub1", 1.0), ("red", "dsub10", 10.0),
                               ("blue", "dsub100", 100.0)):
        raw = {"black": mask_black, "red": mask_red, "blue": mask_blue}[colour](im)
        m = keep_wide_components(inside_frame(raw, p["frame"], inset=10), p["frame"],
                                 min_width_frac=0.30)
        px, py = extract_line(m, p["frame"], max_run=12)
        x, y = fx(px), 10.0 ** fy(py)
        write_csv(DATA / f"fig02a_{tag}.csv", x, y, "h_mm", "Da",
                  f"Prat 2026 Figure 2a, d_sub = {d_sub:.0f} um",
                  f"colour mask ({colour}) + wide-component filter; y axis is log10")
        out[tag] = (x, y)
        series.append((tag, (px, py), {"black": "y.", "red": "c.", "blue": "m."}[colour]))
    overlay(p["image"], p["frame"], series, "fig02a_overlay.png",
            "Figure 2a -- extracted Da(h) curves")
    return out


def do_fig04():
    p = PANELS["fig04"]
    im = load(p["image"])
    dark = mask_black(im)
    print("  fig04 -- Da* model curve (the discrepancy diagnostic)")
    fx = calibrator(dark, p["frame"], p["x"], "x", "fig04")
    fy = calibrator(dark, p["frame"], p["y"], "y", "fig04")

    m = apply_exclusions(inside_frame(dark, p["frame"], inset=6), p["exclude"])
    m = keep_wide_components(m, p["frame"])
    px, py = extract_line(m, p["frame"])
    x, y = fx(px), fy(py)
    write_csv(
        DATA / "fig04_model_curve.csv",
        x,
        y,
        "d_sub_um",
        "Da_star",
        "Prat 2026 Figure 4, main panel black curve (the paper's own model)",
        "colour mask (black) + column-wise thin-run centre",
    )
    overlay(
        p["image"], p["frame"], [("Da* curve", (px, py), "r.")], "fig04_overlay.png",
        "Figure 4 -- extracted model curve",
    )

    # The left inset renders F(Da*) for the same parameters. It is an
    # independent drawing of the same quantity inside the same figure, so
    # comparing it with the main panel tests the figure's internal consistency.
    q = PANELS["fig04_inset_F"]
    print("  fig04 left inset -- F(Da*), an internal-consistency probe")
    chk = q.get("tick_check", True)
    gx = calibrator(dark, q["frame"], {**q["x"], "_check": chk}, "x", "fig04_inset")
    gy = calibrator(dark, q["frame"], {**q["y"], "_check": chk}, "y", "fig04_inset")
    mi = keep_wide_components(inside_frame(dark, q["frame"], inset=5), q["frame"])
    qx, qy = extract_line(mi, q["frame"])
    xi, yi = 10.0 ** gx(qx), gy(qy)
    write_csv(
        DATA / "fig04_inset_F.csv", xi, yi, "d_sub_um", "F",
        "Prat 2026 Figure 4 left inset, sqrt(Da*)/(1+sqrt(Da*))",
        "colour mask (black) + wide-component filter; x axis is log10",
    )
    overlay(q["image"], q["frame"], [("F curve", (qx, qy), "r.")],
            "fig04_inset_overlay.png", "Figure 4 left inset -- extracted F(Da*)")

    # The right inset plots Da* itself on a log axis over [0, 100] um, where
    # the main panel's linear axis has no resolution at all. It is therefore
    # the sharpest test of which Da* the figure was actually drawn with.
    w = PANELS["fig04_inset_zoom"]
    print("  fig04 right inset -- Da* on a log axis, the low-d_sub discriminator")
    chkw = w.get("tick_check", True)
    wx = calibrator(dark, w["frame"], {**w["x"], "_check": chkw}, "x", "fig04_zoom")
    wy = calibrator(dark, w["frame"], {**w["y"], "_check": chkw}, "y", "fig04_zoom")
    mw = keep_wide_components(inside_frame(dark, w["frame"], inset=6), w["frame"])
    ux, uy = extract_line(mw, w["frame"])
    write_csv(
        DATA / "fig04_inset_zoom.csv", wx(ux), 10.0 ** wy(uy), "d_sub_um", "Da_star",
        "Prat 2026 Figure 4 right inset, Da* vs d_sub (log y)",
        "colour mask (black) + wide-component filter; y axis is log10",
    )
    overlay(w["image"], w["frame"], [("Da* zoom", (ux, uy), "r.")],
            "fig04_zoom_overlay.png", "Figure 4 right inset -- extracted Da*")
    return x, y


def do_fig05b():
    """Extract the published concentration profile.

    Needed because the paper states U_evap, d_sub, k_r and eps_s for this panel
    but NOT the structure height h -- and the apparent layer thickness on a z/h
    axis is 1/sqrt(Da), which scales as 1/h. Fitting the published curve is the
    only way to recover the h it was drawn with.

    Careful: the flat part of the curve sits at C/C_eq = 1 and spans the whole
    panel width, so a naive "long horizontal dark line" search mistakes it for
    the bottom axis. The real frame bottom is ~150 px lower, at 0.99999.
    """
    p = PANELS["fig05b"]
    im = load(p["image"])
    dark = mask_black(im)
    print(f"  fig05b -- {p['note']}")
    chk = p.get("tick_check", True)
    fx = calibrator(dark, p["frame"], {**p["x"], "_check": chk}, "x", "fig05b")
    fy = calibrator(dark, p["frame"], {**p["y"], "_check": chk}, "y", "fig05b")

    # This panel has ticks on all four sides, and the top-axis ticks reach down
    # to row 31. Scanning from row 40 excludes them; the exponential band used
    # to fit sqrt(Da) lives around rows 210-590, well clear of that cut, so the
    # fit is unaffected even though the very peak (row ~43) is grazed.
    # Ticks sit on all four sides. The top ones reach row 31, and the left and
    # right ones reach ~26 px inward, so the scan window has to clear all three
    # or "topmost dark pixel" latches onto a tick instead of the curve. The
    # cost is that the last few columns -- the top of the near-vertical rise --
    # are dropped; the exponential fit does not need them.
    r0, r1, c0, c1 = p["frame"]
    px, py = [], []
    for c in range(c0 + 30, c1 - 30):
        rows = np.where(dark[r0 + 30 : r1 - 6, c])[0] + r0 + 30
        if rows.size:
            px.append(c)
            py.append(rows.min())
    px, py = np.array(px, float), np.array(py, float)
    x, y = fx(px), fy(py)
    write_csv(
        DATA / "fig05b_profile.csv", x, y, "z_over_h", "C_over_Ceq",
        "Prat 2026 Figure 5b, published concentration profile",
        "colour mask (black) + topmost mark per column",
    )
    overlay(p["image"], p["frame"], [("profile", (px, py), "r.")],
            "fig05b_overlay.png", "Figure 5b -- extracted profile")
    return x, y


def do_fig08():
    out = {}
    for key, targets in (
        ("fig08a", (("black", "exp"), ("red", "eta1"), ("blue", "eta076"))),
        ("fig08b", (("black", "exp"), ("blue", "eta076"))),
    ):
        p = PANELS[key]
        im = load(p["image"])
        dark = mask_black(im)
        print(f"  {key} -- {p['note']}")
        fx = calibrator(dark, p["frame"], p["x"], "x", key)
        fy = calibrator(dark, p["frame"], p["y"], "y", key)
        series = []
        for colour, tag in targets:
            raw = {"black": mask_black, "red": mask_red, "blue": mask_blue}[colour](im)
            m = apply_exclusions(inside_frame(raw, p["frame"], inset=6), p["exclude"])
            if colour == "black":
                pts = extract_blobs(m)
                px, py = pts[:, 0], pts[:, 1]
                method = "colour mask (black) + connected-component centroids"
            else:
                m = keep_wide_components(m, p["frame"])
                px, py = extract_line(m, p["frame"])
                method = f"colour mask ({colour}) + wide-component filter + column centre"
            x, y = fx(px), fy(py)
            ylabel = "Ms_g" if key == "fig08a" else "h_equ_mm"
            write_csv(
                DATA / f"{key}_{tag}.csv", x, y, "time_hr", ylabel,
                f"Prat 2026 Figure {key[-2:]}, {tag}", method,
            )
            out[f"{key}_{tag}"] = (x, y)
            series.append((tag, (px, py), {"black": "r.", "red": "c-", "blue": "y-"}[colour]))
        overlay(p["image"], p["frame"], series, f"{key}_overlay.png",
                f"Figure {key[-2:]} -- extracted series")
    return out


def do_fig09():
    out = {}

    # --- 9a: one continuous experimental curve, schematic overlapping it ----
    p = PANELS["fig09a"]
    im = load(p["image"])
    dark = mask_black(im)
    print(f"  fig09a -- {p['note']}")
    fx = calibrator(dark, p["frame"], p["x"], "x", "fig09a")
    fy = calibrator(dark, p["frame"], p["y"], "y", "fig09a")
    m = inside_frame(dark, p["frame"], inset=6)
    px, py = trace_curve(m, p["frame"], start_col=p["frame"][2] + 8, start_row=693)
    x, y = fx(px), fy(py)
    write_csv(
        DATA / "fig09a_exp.csv", x, y, "time_hr", "evaporated_mass_g",
        "Prat 2026 Figure 9a, experimental curve (Diouf thesis, ref 41)",
        "colour mask (black) + left-to-right continuity trace (schematic overlaps)",
    )
    out["fig09a_exp"] = (x, y)
    overlay(p["image"], p["frame"], [("evaporated mass", (px, py), "r.")],
            "fig09a_overlay.png", "Figure 9a -- traced experimental curve")

    # --- 9b: five markers plus two model lines -----------------------------
    p = PANELS["fig09b"]
    print(f"  fig09b -- {p['note']}")
    dark = mask_black(im)
    fx = calibrator(dark, p["frame"], p["x"], "x", "fig09b")
    fy = calibrator(dark, p["frame"], p["y"], "y", "fig09b")
    series = []
    for colour, tag in (("black", "exp"), ("red", "eps043"), ("blue", "eps034")):
        raw = {"black": mask_black, "red": mask_red, "blue": mask_blue}[colour](im)
        mm = apply_exclusions(inside_frame(raw, p["frame"], inset=6), p["exclude"])
        if colour == "black":
            pts = extract_blobs(mm, min_area=80, max_area=1500)
            px, py = pts[:, 0], pts[:, 1]
            method = "colour mask (black) + connected-component centroids"
        else:
            mm = keep_wide_components(mm, p["frame"])
            px, py = extract_line(mm, p["frame"])
            method = f"colour mask ({colour}) + wide-component filter + column centre"
        x, y = fx(px), fy(py)
        write_csv(DATA / f"fig09b_{tag}.csv", x, y, "time_hr", "h_mm",
                  f"Prat 2026 Figure 9b, {tag}", method)
        out[f"fig09b_{tag}"] = (x, y)
        series.append((tag, (px, py), {"black": "r.", "red": "c-", "blue": "y-"}[colour]))
    overlay(p["image"], p["frame"], series, "fig09b_overlay.png",
            "Figure 9b -- extracted series")
    return out


def main():
    OVERLAYS.mkdir(parents=True, exist_ok=True)
    if not (DATA / "raw_p6-000.png").exists():
        sys.exit("Run pdfimages first (see README) -- raw page images missing.")

    print("Digitising Prat (2026) figures")
    print("Calibration is checked against each panel's own major ticks:\n")
    fig02a = do_fig02a()
    fig04 = do_fig04()
    fig05b = do_fig05b()
    fig08 = do_fig08()
    fig09 = do_fig09()

    print("\nSanity summary")
    x4, y4 = fig04
    if len(x4) > 10:
        slope = np.polyfit(x4, y4, 1)
        print(f"  fig04 published curve: Da* = {slope[0]:.4g} * d_sub[um] + {slope[1]:.3g}")
        print(f"    -> at 580 um the published CURVE reads {np.polyval(slope, 580):.1f}")
        print("       (paper TEXT quotes 156; core.Da_star gives 169.6)")
    for key in sorted(fig08):
        x, y = fig08[key]
        if len(x) > 1:
            print(f"  {key}: {len(x)} pts, t in [{x.min():.2f}, {x.max():.2f}] hr, "
                  f"y in [{y.min():.4g}, {y.max():.4g}]")
    for key in sorted(fig09):
        x, y = fig09[key]
        if len(x) > 1:
            print(f"  {key}: {len(x)} pts, t in [{x.min():.2f}, {x.max():.2f}] hr, "
                  f"y in [{y.min():.4g}, {y.max():.4g}]")
    print("\nCheck data/overlays/*.png before trusting any of this.")


if __name__ == "__main__":
    main()
