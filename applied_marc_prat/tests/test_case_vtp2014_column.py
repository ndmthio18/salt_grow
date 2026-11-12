"""The resolved bead column of the Veran-Tissoires & Prat (2014) series (``cases.vtp2014_column``).

Physical assertions on the network and on one short run; the long runs are execution records
(runs/E1_pre/vtp2014_column.py), not tests.
"""
import numpy as np
import pytest

from saltv2.cases import vtp2014_column as case
from saltv2.coupling import driver3d
from saltv2.coupling.driver3d import col
from saltv2.coupling.interface import PoreOpenings
from saltv2.growth.grid import Grid
from saltv2.materials.brine import C_SAT
from saltv2.pnm import PoreLayer


@pytest.fixture(scope="module")
def columns(tmp_path_factory):
    directory = tmp_path_factory.mktemp("columns")
    return {name: case.write_dgf(name, directory) for name in case.BEADS}


@pytest.mark.parametrize("name", list(case.BEADS))
def test_the_column_is_the_papers_network_at_that_bead_size(columns, name):
    d_b = case.BEADS[name]
    layer = PoreLayer.from_dgf(columns[name])
    assert layer.N * d_b == pytest.approx(case.L_COLUMN, abs=0.5 * d_b)      # 35 mm of beads
    assert layer.dz == pytest.approx(d_b, rel=1e-9)                          # node distance d_b
    assert np.allclose(layer.V, 0.363 * d_b ** 3, rtol=1e-12)                # eps d_b^3
    assert np.allclose(layer.V, (2 * layer.r) ** 3, rtol=1e-12)              # a cube, whole
    assert layer.A_open[layer.surface][0] / d_b ** 2 == pytest.approx(0.363 ** (2 / 3), rel=1e-12)
    for _, _, rt, L in layer.lat:
        assert rt == pytest.approx(0.15 * d_b, rel=1e-12)
        assert L + 2 * layer.r[0] == pytest.approx(d_b, rel=1e-12)           # pore + throat = spacing
    assert layer.rt_res[layer.bottom][0] == pytest.approx(0.15 * d_b, rel=1e-12)


def test_the_bead_size_does_not_change_the_supply(columns):
    """Flow ~ d^2, throat conductance ~ d, number of throats ~ 1/d: the Peclet number of the bare
    column is the same for every class, and far above the ln(C_sat / C0) below which a crystal on
    the surface pore could not be fed at all -- so the classes differ in the salt structure only."""
    pe = np.array([case.bare_column_peclet(name, columns[name]) for name in case.BEADS])
    assert np.ptp(pe) / pe.mean() < 0.02
    assert pe.min() > 10 * np.log(C_SAT / case.C0)
    assert pe.mean() == pytest.approx(1.0, abs=0.05)             # theirs, on calibrated throats: 1.48


@pytest.mark.parametrize("vpp", [case.VPP, case.VPP_FINE])
def test_the_voxels_resolve_the_opening_they_are_meant_to(columns, vpp):
    kw = case.case("1mm", columns["1mm"], vpp=vpp)
    layer = PoreLayer.from_dgf(kw["substrate_dgf"], L_box=kw["substrate_L_box"])
    a = layer.L_box / vpp
    cells = 5 * vpp // 7
    assert cells * a == pytest.approx(2 * layer.r[layer.surface][0], rel=2e-3)
    surf = layer.surface
    openings = PoreOpenings(Grid((vpp, vpp, 4), a), layer.x[surf], layer.y[surf], layer.r[surf],
                            layer.L_box)
    assert openings.voxel_area()[0] == pytest.approx((cells * a) ** 2, rel=1e-12)


def test_the_coarse_column_nucleates_and_its_books_close(columns):
    res = driver3d.run(**dict(case.case("1mm", columns["1mm"]), t_end=6 * 3600.0, verbose=False))
    assert res["params"]["d_sub"] == case.BEADS["1mm"]
    assert np.isfinite(res["t_nuc"]).all() and res["t_nuc"][0] < 3 * 3600.0
    top_down = np.argsort(-res["pore_z"])
    C = res["pore_C"][-1][top_down]
    assert C[0] == pytest.approx(C_SAT, rel=1e-12)               # under salt: held at saturation
    assert np.all(np.diff(C) <= 1e-15) and C[-1] >= case.C0      # falling to the reservoir
    for name in ("salt_network", "salt_interface", "salt_precipitation"):
        assert abs(res["ledger"][f"relative_residual_{name}"]) < 1e-12
    assert np.all(np.diff(col(res["hist"], "t_h")) > 0)
