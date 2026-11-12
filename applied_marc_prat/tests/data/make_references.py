"""Capture the reference records for tests/test_regression.py.

    OMP_NUM_THREADS=1 python tests/data/make_references.py            # all three
    OMP_NUM_THREADS=1 python tests/data/make_references.py ref_2d_cv  # one

Run this ONLY from a tree whose behaviour is the one to be pinned, on the platform the suite runs
on (frodo since 2026-09-15: the 2-D references from a clone at a0db925 = checkpoint 95c1acb plus
this harness, ref_3d_small from a clone at f0ff577, square pore openings; see reference_cases.py).  Re-running it on a later tree replaces the
reference with whatever that tree does, so it must never be part of a fix for a failing
regression test.  Each file records numpy/scipy/platform and the git hash it came from.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tests.reference_cases import CASES, DATA, call, capture, meta       # noqa: E402


def main(names):
    DATA.mkdir(exist_ok=True)
    for name in names:
        kind, kwargs = CASES[name]
        res = call(kind, kwargs)
        rec = capture(kind, res)
        rec.update({f"meta_{k}": np.array(v) for k, v in meta().items()})
        out = DATA / f"{name}.npz"
        np.savez_compressed(out, **rec)
        print(f"{name}: {rec['hist'].shape[0]} rows, {out.stat().st_size / 1e3:.0f} kB, git {meta()['git'][:7]}")


if __name__ == "__main__":
    main(sys.argv[1:] or list(CASES))
