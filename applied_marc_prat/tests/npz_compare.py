"""Exact comparison of two run records (reference vs fresh), by name.

    python -m tests.npz_compare ref.npz new.npz

Compares every history column present in both files (by ``hist_columns`` name when both carry
one, else positionally), every ledger channel present in both, the per-voxel arrays and the
step count.  Equality is exact (``np.array_equal`` with ``equal_nan=True``: ``j_bare_over_plane``
is NaN once no bare opening is left).  Nothing is rounded and no tolerance exists: the cleanup
of 2026-09-14 promises the production path bit for bit, and this is the instrument.
"""
from __future__ import annotations

import sys

import numpy as np

VOXEL_ARRAYS = ("state", "eps_s", "m_precipitated", "m_skin", "t_birth", "t_sealed", "d_sub_col",
                "t_nuc")


def _ledger(d):
    if "ledger_names" in d.files:
        return dict(zip(d["ledger_names"].tolist(), d["ledger_values"].tolist()))
    if "ledger" in d.files:                       # driver npz: repr of the report dict
        return eval(str(d["ledger"]), {"nan": np.nan, "inf": np.inf})
    return {}


def compare_run(ref, new):
    """Return a list of human-readable mismatch lines; empty means identical."""
    bad = []
    # history
    h0, h1 = np.asarray(ref["hist"], float), np.asarray(new["hist"], float)
    if "hist_columns" in ref.files and "hist_columns" in new.files:
        c0, c1 = [str(c) for c in ref["hist_columns"]], [str(c) for c in new["hist_columns"]]
        common = [c for c in c0 if c in c1]
        if h0.shape[0] != h1.shape[0]:
            bad.append(f"hist rows {h0.shape[0]} vs {h1.shape[0]}")
        else:
            for c in common:
                a, b = h0[:, c0.index(c)], h1[:, c1.index(c)]
                if not np.array_equal(a, b, equal_nan=True):
                    k = int(np.argmax(~((a == b) | (np.isnan(a) & np.isnan(b)))))
                    bad.append(f"hist[{c}] differs first at row {k}: {a[k]!r} vs {b[k]!r}")
        missing = [c for c in c0 if c not in c1]
        if missing:
            bad.append(f"note: columns absent in new (not compared): {missing}")
    else:
        if h0.shape != h1.shape:
            bad.append(f"hist shape {h0.shape} vs {h1.shape}")
        elif not np.array_equal(h0, h1, equal_nan=True):
            i, j = np.argwhere(~((h0 == h1) | (np.isnan(h0) & np.isnan(h1))))[0]
            bad.append(f"hist[{i},{j}] {h0[i, j]!r} vs {h1[i, j]!r}")
    # ledger
    l0, l1 = _ledger(ref), _ledger(new)
    for k in sorted(set(l0) & set(l1)):
        if not (l0[k] == l1[k] or (np.isnan(l0[k]) and np.isnan(l1[k]))):
            bad.append(f"ledger[{k}] {l0[k]!r} vs {l1[k]!r}")
    # voxel arrays and counts
    for k in VOXEL_ARRAYS + ("n_steps",):
        if k in ref.files and k in new.files:
            a, b = np.asarray(ref[k]), np.asarray(new[k])
            if a.shape != b.shape:
                bad.append(f"{k} shape {a.shape} vs {b.shape}")
            elif not np.array_equal(a, b, equal_nan=True):
                n = int(np.sum(~((a == b) | (np.isnan(a) & np.isnan(b)))) if a.dtype.kind == "f"
                        else np.sum(a != b))
                bad.append(f"{k} differs in {n} entries")
    return bad


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        sys.exit(__doc__)
    ref, new = (np.load(p, allow_pickle=True) for p in argv)
    bad = compare_run(ref, new)
    real = [b for b in bad if not b.startswith("note:")]
    for line in bad:
        print(line)
    print("IDENTICAL" if not real else f"{len(real)} MISMATCH(ES)")
    sys.exit(0 if not real else 1)


if __name__ == "__main__":
    main()
