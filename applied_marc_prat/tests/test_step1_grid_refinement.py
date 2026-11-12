"""Step 1B acceptance: the 10 um grid represents the Step 1 geometry exactly (no solve), and the
committed refinement record satisfies the conservation and consistency requirements."""
import json
from pathlib import Path

import numpy as np
import pytest

from saltv2.cases import step1_frozen_salt as s1

REC = Path(__file__).resolve().parents[1] / "runs" / "S_surface" / "step1_grid_refinement" / "results" / "step1_refinement.json"


@pytest.mark.parametrize("a, n_side, n_h", ((20e-6, 27, 15), (10e-6, 54, 30)))
def test_grid_represents_the_physical_pillar_exactly(a, n_side, n_h):
    grid, foot, L_ext = s1.build(s1.H_BOX, a)
    assert int(foot.sum()) == n_side ** 2
    pillar, faces = s1.wet_faces(grid, foot, 0.3e-3, sides=True)
    assert int(pillar.sum()) == n_side ** 2 * n_h
    assert int(faces[(2, -1)].sum()) == n_side ** 2
    n_side_faces = sum(int(m.sum()) for k, m in faces.items() if k[0] != 2)
    assert n_side_faces == 4 * n_side * n_h
    assert n_side ** 2 * a * a == pytest.approx(s1.PORE_SIDE ** 2, rel=1e-12)
    assert n_side_faces * a * a == pytest.approx(4 * s1.PORE_SIDE * 0.3e-3, rel=1e-12)
    # each side face belongs to a gas cell, and no face is counted twice
    for k, m in faces.items():
        assert not np.any(m & pillar)
    assert L_ext == pytest.approx(18.6e-3)


@pytest.mark.skipif(not REC.exists(), reason="run runs/S_surface/step1_grid_refinement/run.py first")
def test_refinement_record_is_conservative_and_consistent():
    r = json.loads(REC.read_text())
    for lab in ("coarse", "fine"):
        a1, a2 = r[lab]["A1"], r[lab]["A2"]
        assert a1["E_side"] == 0.0 and a1["E_total"] == a1["E_top"]
        assert a2["E_total"] == pytest.approx(a2["E_top"] + a2["E_side"], rel=1e-15)
        for o in (a1, a2):
            assert abs(o["imbalance"]) <= 1e-8 and o["residual"] <= 1e-10
            assert o["A_top_wet"] == pytest.approx(s1.PORE_SIDE ** 2, rel=1e-12)
        assert a2["A_side_wet"] == pytest.approx(4 * s1.PORE_SIDE * 0.3e-3, rel=1e-12)
    assert r["fine"]["A1"]["n_unknown"] == 70 * 70 * 200 - 54 * 54 * 30
    p = r["paired"]
    for lab in ("coarse", "fine"):
        assert p[lab]["dE"] == pytest.approx(r[lab]["A2"]["E_total"] - r[lab]["A1"]["E_total"])
        assert 0.0 < p[lab]["f_side"] < 1.0
