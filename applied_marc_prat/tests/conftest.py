"""Test configuration: the import path and the shared small production run.

Everything the tests import is in ``saltv2`` (and ``tests`` itself), so the repository root is
the only path entry; ``saltv2.growth.prat`` finds prat2026/ on its own.  The small boxes are
``saltv2.cases.small``, re-exported here because the test modules import them from
``tests.conftest``; ``run_3d_small`` is the very case pinned bit for bit in
tests/data/ref_3d_small.npz (test_regression.py) and is shared as a session fixture because it
takes ~10 s.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from saltv2.cases.small import CASE_3D, SMALL_2D, SMALL_3D                    # noqa: E402,F401


@pytest.fixture(scope="session")
def run_3d_small():
    from saltv2.coupling import driver3d as v
    return v.run(**CASE_3D)
