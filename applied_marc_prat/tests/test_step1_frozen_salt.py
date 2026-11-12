"""Step 1 acceptance (verify/step1): evaporation from a prescribed wet salt pillar,
saltv2.cases.step1_frozen_salt.  Gas only, 50 um grid for speed (seconds); the 20 um
primary results are the execution record under runs/S_surface/step1_frozen_salt."""
import numpy as np
import pytest

from saltv2.cases import step0_gas as s0
from saltv2.cases import step1_frozen_salt as s1

A = s1.A_VOXEL_COARSE


def test_face_selective_assembly_reproduces_the_production_kernel():
    cc = s1.cross_check(h=0.3e-3, a=A)
    assert abs(cc["A2_E_rel"]) < 1e-11 and cc["A2_X_rel_max"] < 1e-13
    assert abs(cc["A2_E_top_rel"]) < 1e-11
    assert cc["A2_E_base_kernel"] == 0.0                         # the covered opening never evaporates
    assert abs(cc["A0_E_rel"]) < 1e-11 and cc["A0_X_rel_max"] < 1e-13


def test_A0_reproduces_step0_on_the_same_grid():
    o, *_ = s1.run_case("A0", 0.0, a=A)
    ref, _ = s0.solve_pore(s1.H_BOX)
    assert o["E_total"] == pytest.approx(ref["E"], rel=1e-12)
    assert o["E_top"] == o["E_total"] and o["E_side"] == 0.0


def test_zero_height_continuity_A1_A2_equal_A0():
    o0, X0, _, m0 = s1.run_case("A0", 0.0, a=A)
    for cond in ("A1", "A2"):
        o, X, _, m = s1.run_case(cond, 0.0, a=A)
        assert set(m["faces"]) == set(m0["faces"]) == {(2, -1)}
        assert np.array_equal(m["faces"][(2, -1)], m0["faces"][(2, -1)])
        assert o["E_total"] == o0["E_total"] and np.array_equal(X, X0)


def test_A1_side_faces_are_inert_and_A2_faces_are_counted_once():
    h = 0.3e-3; n_h = int(round(h / A))
    o1, X1, j1, m1 = s1.run_case("A1", h, a=A)
    assert set(m1["faces"]) == {(2, -1)} and o1["E_side"] == 0.0 and o1["n_side_faces"] == 0
    # inert sides: the gas cells next to the pillar carry no Dirichlet term, so the lateral
    # gradient INTO the pillar is not a flux; the field next to the sides is below X_wet
    pillar = m1["pillar"]; grid, foot, _ = s1.build(s1.H_BOX, A)
    side_cells = np.zeros(grid.shape, bool)
    for ax in (0, 1):
        for st in (+1, -1):
            side_cells |= grid.shift(pillar, ax, st, False) & ~pillar
    assert np.all(X1[side_cells] < o1["X_wet"])
    o2, X2, j2, m2 = s1.run_case("A2", h, a=A)
    n_side = int(np.sqrt(foot.sum()))
    assert o2["n_side_faces"] == 4 * n_side * n_h
    assert o2["A_side_wet"] == pytest.approx(4 * n_side * A * h)
    assert np.all(X2[side_cells] < o2["X_wet"]) and o2["E_side"] > 0.0
    assert o2["E_total"] == o2["E_top"] + o2["E_side"]
    # the opening under the pillar is not a gas cell
    assert np.all(np.isnan(X1[foot, 0])) and np.all(np.isnan(X2[foot, 0]))


@pytest.mark.parametrize("cond", ("A1", "A2"))
def test_vapour_mass_conservation(cond):
    o, *_ = s1.run_case(cond, 0.3e-3, a=A)
    assert abs(o["imbalance"]) < 1e-10 and o["residual"] < 1e-12


def test_20um_grid_resolves_the_opening_exactly():
    grid, foot, L_ext = s1.build(s1.H_BOX, s1.A_VOXEL)
    assert grid.shape == (35, 35, 100) and int(foot.sum()) == 27 * 27
    assert foot.sum() * s1.A_VOXEL ** 2 == pytest.approx(s1.PORE_SIDE ** 2, rel=1e-12)
    assert L_ext == pytest.approx(18.6e-3)
