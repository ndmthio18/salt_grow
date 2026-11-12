"""B' (2026-09-22): the kinetic share follows the crystals of each layer, and a voxel born at the
clamp F_min is dense and sealed at birth.

Regime: evaporation-wicking, the small 3-D box.  The 27.5 um class of the hierarchical structure
has F = 0.850 on layer 1 and 0.800 on layer 2, below F_min = (1 - 0.20)/(1 - 0.02) = 0.816, so
the structure grows porous over the first bead diameter and turns dense there.  Off, the same
box grows porous throughout (layer-1 F on every layer), as every record before this date did.
"""
import numpy as np
import pytest

from saltv2.cases.small import SMALL_3D
from saltv2.coupling.driver3d import col, run
from saltv2.growth.crust import HIERARCHICAL, SINGLE_SCALE, BirthClosure
from saltv2.growth.grid import SALT_SEALED
from saltv2.growth.prat import prat

KW = dict(SMALL_3D, c_k=1.0, t_end=6 * 3600.0, eps_min=0.02, salt_structure=HIERARCHICAL,
          d_sub=27.5e-6, verbose=False)


def test_layer_factor_is_the_pore_ratio_power_and_one_when_off():
    on = BirthClosure(d_sub=27.5e-6, structure=HIERARCHICAL, layer_av=True)
    off = BirthClosure(d_sub=27.5e-6, structure=HIERARCHICAL)
    n = np.arange(1, 6)
    assert np.allclose(on.layer_da_factor(n), 0.5 ** (n - 1))
    assert np.all(off.layer_da_factor(n) == 1.0)
    single = BirthClosure(d_sub=27.5e-6, structure=SINGLE_SCALE, layer_av=True)
    assert np.all(single.layer_da_factor(n) == 1.0)         # r = 1: nothing to follow
    assert on.describe()["layer_av"] is True and off.describe()["layer_av"] is False
    # the numbers of runs/E4_crust1d and README section 19: F by layer at 27.5 um
    Da1 = prat.Da_star(eps_s=0.20, d_sub=27.5e-6, kr=prat.kr_at_temperature(295.15), Ds=1.54e-9)
    F = prat.F(Da1 * on.layer_da_factor(np.array([1, 2, 3])))
    assert F == pytest.approx([0.850, 0.800, 0.739], abs=2e-3)


def test_off_keeps_growing_and_on_seals_the_second_layer():
    off = run(**KW)                                       # layer-1 F everywhere: porous throughout
    assert col(off["hist"], "n_sealed_cells").max() == 0
    assert col(off["hist"], "h_max")[-1] > 2 * 27.5e-6
    # a class whose FIRST layer is already below the limit still raises when the switch is off
    with pytest.raises(ValueError, match="crust limit"):
        run(**dict(KW, d_sub=12e-6, t_end=600.0))
    res = run(**KW, layer_av=True)
    p, c, hh = res["params"], res["crust"], res["hist"]
    assert p["layer_av"] is True and p["F_min"] == pytest.approx(0.8 / 0.98)
    sealed = c.state == SALT_SEALED
    assert p["n_dense_born"] > 0 and sealed.sum() >= p["n_dense_born"]
    # dense birth happens on layer 2 and above, never on layer 1
    dense_born = sealed & (c.eps_birth <= 0.02 + 1e-3)
    assert dense_born.sum() == p["n_dense_born"] and c.generation[dense_born].min() >= 2
    # a dense-born voxel carries the shrunk radius a_2 R_P sqrt(eps_min / eps_s) = 76 nm ...
    a2 = float(BirthClosure(d_sub=27.5e-6, structure=HIERARCHICAL).pore_radius(2))
    assert np.allclose(c.r_birth[dense_born], 0.5 * a2 * np.sqrt(c.eps_birth[dense_born] / 0.20), rtol=1e-12)
    assert c.r_birth[dense_born].max() < 80e-9
    # ... and the structure stops: no wet salt is left, the evaporation has fallen, and the
    # height froze at layer 2 (36.5 um: five voxels of 7.29 um)
    assert col(hh, "V_wet")[-1] == 0.0
    assert col(hh, "E_total").min() < 0.75 and col(hh, "E_total")[-1] < 0.9 < col(hh, "E_total")[0]
    assert col(hh, "h_max")[-1] == pytest.approx(5 * res["a"], rel=1e-12)
    # the books still close
    for k, v in res["ledger"].items():
        if k.startswith("relative_residual"):
            assert abs(v) < 1e-10, (k, v)
