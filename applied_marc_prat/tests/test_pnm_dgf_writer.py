"""``pnm.write_column_dgf`` writes what ``PoreLayer.from_dgf`` reads: no physics, a round trip."""
import numpy as np
import pytest

from saltv2.pnm import DGF_DUMUX_1D, PoreLayer, write_column_dgf

FIELDS = ("N", "n_layers", "n_side", "n_surface", "dz", "L_box", "reservoir_pore")
ARRAYS = ("r", "V", "A_open", "z", "height", "surface", "bottom", "rt_res", "Lt_res")


def test_it_writes_the_dumux_column_back(tmp_path):
    """The shipped 1-D file, rebuilt from its six numbers, must read back as the same network."""
    theirs = PoreLayer.from_dgf(DGF_DUMUX_1D)
    path = write_column_dgf(tmp_path / "column.dgf", n_pores=5, spacing=0.0006, pore_radius=0.0002,
                            pore_volume=6.4e-11, throat_radius=1.1e-4, throat_length=0.0002)
    ours = PoreLayer.from_dgf(path)
    for name in FIELDS:
        assert getattr(ours, name) == pytest.approx(getattr(theirs, name), rel=1e-12), name
    for name in ARRAYS:
        assert np.allclose(getattr(ours, name), getattr(theirs, name), rtol=1e-12, atol=0.0), name
    assert len(ours.lat) == len(theirs.lat)
    for (a, b, rt, L), (a0, b0, rt0, L0) in zip(ours.lat, theirs.lat):
        assert (a, b) == (a0, b0) and rt == rt0 and L == L0


def test_the_radius_and_the_volume_of_a_pore_are_independent(tmp_path):
    path = write_column_dgf(tmp_path / "c.dgf", n_pores=40, spacing=1e-4, pore_radius=3.5e-5,
                            pore_volume=2.0e-13, throat_radius=1.5e-5, throat_length=2.9e-5)
    layer = PoreLayer.from_dgf(path)
    assert layer.N == 40 and layer.surface.size == 1 and layer.bottom.size == 1
    assert np.all(layer.V == 2.0e-13) and np.all(layer.A_open == (2 * 3.5e-5) ** 2)
    assert layer.z.max() == 0.0 and layer.z.min() == pytest.approx(-39e-4, rel=1e-12)
    assert len(layer.lat) == 39
    assert layer.rt_res[layer.bottom][0] == 1.5e-5 and layer.Lt_res[layer.bottom][0] == 2.9e-5
    assert np.count_nonzero(layer.rt_res) == 1          # only the bottom pore touches the reservoir


@pytest.mark.parametrize("bad", [dict(n_pores=0), dict(spacing=0.0), dict(pore_volume=-1.0)])
def test_it_refuses_a_column_that_is_not_one(tmp_path, bad):
    good = dict(n_pores=3, spacing=1e-4, pore_radius=3e-5, pore_volume=1e-13,
                throat_radius=1e-5, throat_length=3e-5)
    with pytest.raises(ValueError):
        write_column_dgf(tmp_path / "bad.dgf", **{**good, **bad})
