"""Supply-limit root finding of the pre-experiment (runs/E1_pre): the first failing height must be the
lowest crossing of need(h) = Pc(h), not a later one after a layer-boundary recovery."""
import numpy as np
import pytest

from saltv2.growth.crust import LEGACY_MIXED, SINGLE_SCALE
from saltv2.growth.supply import first_failure_height, supply_need

STACK = LEGACY_MIXED        # the numbers below were reviewed on the halving stack


def scan_first_crossing(U, d_sub, s_foot, c_k, h_top, n=40000):
    hs = np.linspace(1e-7, h_top, n)
    f = np.array([supply_need(h, U, d_sub, s_foot, c_k, STACK)[0] - supply_need(h, U, d_sub, s_foot, c_k, STACK)[1] for h in hs])
    i = np.argmax(f > 0)
    return hs[i] if f[i] > 0 else np.inf


def test_first_crossing_with_recovery_at_layer_boundary():
    # reviewed case: fails at ~4.09 mm, recovers at the 4.28 mm layer boundary, fails again at ~4.81 mm
    U, d_sub, s_foot, c_k = 3.2e-7, 68e-6, 0.0, 3e-3
    h1 = first_failure_height(U, d_sub, s_foot, c_k, 6e-3, STACK)
    ref = scan_first_crossing(U, d_sub, s_foot, c_k, 6e-3)
    assert h1 == pytest.approx(ref, abs=2e-6)
    assert 4.0e-3 < h1 < 4.2e-3
    # a single bisection over [0, 6 mm] lands on the later crossing -- the behaviour being replaced
    lo, hi = 0.0, 6e-3
    for _ in range(40):
        mid = 0.5 * (lo + hi); nd, pc = supply_need(mid, U, d_sub, s_foot, c_k, STACK)
        lo, hi = (lo, mid) if nd > pc else (mid, hi)
    assert lo > 4.5e-3


@pytest.mark.parametrize("U,d_sub,c_k", [(1e-7, 34e-6, 3e-3), (1e-6, 135e-6, 3e-3), (3e-8, 1080e-6, 1e-6)])
def test_first_crossing_matches_fine_scan(U, d_sub, c_k):
    h_top = 8e-3
    h1 = first_failure_height(U, d_sub, 200.0, c_k, h_top, STACK)
    ref = scan_first_crossing(U, d_sub, 200.0, c_k, h_top)
    if np.isinf(ref):
        assert np.isinf(h1)
    else:
        assert h1 == pytest.approx(ref, abs=2e-6)


def test_no_failure_below_top_returns_inf():
    assert np.isinf(first_failure_height(1e-9, 1080e-6, 0.0, 3e-3, 1e-3, STACK))


def test_negative_or_zero_demand_never_returns_negative_height():
    # a cluster whose faces net-condense (U < 0) must not produce a negative failure height
    for U in (-1e-6, -1e-9, 0.0):
        h = first_failure_height(U, 34e-6, 0.0, 3e-3, 4e-3, STACK)
        assert np.isinf(h) or h >= 0.0
    # gravity alone still limits a very tall column: with U = 0 the crossing exists only far above 4 mm
    assert np.isinf(first_failure_height(0.0, 34e-6, 0.0, 3e-3, 4e-3, STACK))


def test_a_single_scale_structure_fails_where_its_one_layer_does():
    """Without a hierarchy need(h) is ONE straight line and the threshold one number, so the
    first failure is the closed form h = (Pc - s_foot) / (rho g + mu U / k)."""
    from saltv2.materials.brine import GAMMA, MU, RHO_L
    from saltv2.materials.constants import G
    U, d_sub, s_foot, c_k = 3.2e-7, 68e-6, 50.0, 1e-7
    a1 = 0.14 * d_sub / 4
    k = c_k * 2 * a1 ** 3 / (12 * 0.14 * d_sub)
    h = (2 * GAMMA / a1 - s_foot) / (RHO_L * G + MU * U / k)
    assert first_failure_height(U, d_sub, s_foot, c_k, 10 * h, SINGLE_SCALE) == pytest.approx(h, rel=1e-12)
    assert np.isinf(first_failure_height(U, d_sub, s_foot, c_k, 0.5 * h, SINGLE_SCALE))
    need, Pc = supply_need(h, U, d_sub, s_foot, c_k, SINGLE_SCALE)
    assert need == pytest.approx(Pc, rel=1e-12)


def test_the_gate_has_no_default_structure():
    with pytest.raises(TypeError):
        first_failure_height(3.2e-7, 68e-6, 0.0, 3e-3, 6e-3)
