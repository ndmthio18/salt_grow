"""The DuMuX 1-D network as the substrate of the 3-D closed loop (``saltv2.cases.dumux1d``).

Physical assertions only: the network must be the file's, the ledgers must close, and at the
reference case's own far field the surface pore must nucleate.  Nothing is pinned bit for bit.
"""
import numpy as np
import pytest

from saltv2.cases import dumux1d
from saltv2.pnm import PoreLayer


def test_the_network_is_the_dgf_file_and_not_a_template():
    """pnm_1d_homogeneous.dgf: cubic pores of r = 200 um (V = (2 r)^3 = 6.4e-11 m^3) 0.6 mm apart,
    throats of r = 110 um and length 0.6 mm - 2 r = 200 um; the pore labelled 1 is the reservoir
    and is therefore not a degree of freedom."""
    layer = PoreLayer.from_dgf(dumux1d.DGF)
    assert layer.N == 5 and layer.n_surface == 1 and len(layer.lat) == 4
    assert np.allclose(layer.r, 200e-6) and np.allclose(layer.V, 6.4e-11)
    assert layer.dz == pytest.approx(0.6e-3)
    assert all(rt == pytest.approx(110e-6) and L == pytest.approx(200e-6) for _, _, rt, L in layer.lat)
    assert np.count_nonzero(layer.rt_res) == 1                   # one reservoir throat, at the bottom
    assert layer.rt_res[layer.bottom] == pytest.approx(110e-6)
    assert layer.Lt_res[layer.bottom] == pytest.approx(200e-6)
    assert layer.z[layer.surface] == 0.0 and layer.z[layer.bottom] == pytest.approx(-2.4e-3)


def test_the_case_values_are_the_dumux_inputs():
    assert dumux1d.C0_SATURATED == pytest.approx(0.25402, abs=1e-5)     # x_NaCl = 0.095
    assert dumux1d.C0_DILUTE == pytest.approx(0.14583, abs=1e-5)        # x_NaCl = 0.05
    assert dumux1d.RH_ENV == pytest.approx(0.19959, abs=1e-5)


def test_the_closed_loop_runs_on_it_and_both_ledgers_close():
    from saltv2.coupling import driver3d
    res = driver3d.run(**dumux1d.RUN_3D, verbose=False)
    h = res["hist"]
    assert np.isfinite(res["t_nuc"]).all() and res["t_nuc"][0] < 3600.0     # nucleates within the hour
    assert driver3d.col(h, "h_max")[-1] > 0.0                               # and grows
    assert np.abs(driver3d.col(h, "res_water")).max() < 1e-12
    assert np.abs(driver3d.col(h, "res_salt")).max() < 1e-12
    assert abs(res["ledger"]["relative_residual_salt_network"]) < 1e-10


# ---- the column coupled to the voxel domain in its one-dimensional limit (RUN_3D_CHECK) ---------

@pytest.fixture(scope="module")
def check_run():
    from saltv2.coupling import driver3d
    return driver3d.run(**dumux1d.RUN_3D_CHECK, verbose=False)


def test_a_box_as_wide_as_the_opening_is_the_dumux_boundary_layer(check_run):
    """The opening fills the box floor, so the voxel vapour solve is one-dimensional and has to
    return rho_g D_v A (a_w(C0) - RH) W_sat / delta -- their top boundary condition -- exactly."""
    from saltv2.pnm import boundary_layer_evaporation
    expected = boundary_layer_evaporation(dumux1d.C0_DILUTE, dumux1d.OPENING ** 2,
                                          dumux1d.DELTA_BL, dumux1d.RH_ENV)
    assert check_run["M_ref"] == pytest.approx(float(expected), rel=1e-10)
    assert np.isnan(check_run["t_nuc"]).all()                  # the dilute column never nucleates


def test_the_pore_histories_belong_to_the_hist_rows(check_run):
    res = check_run
    layer = res["layer"]
    assert res["pore_C"].shape == (len(res["hist"]), layer.N)
    assert res["M_liquid_in"].shape == (len(res["hist"]),)
    assert np.array_equal(res["pore_z"], layer.z)
    first = res["pore_C"][0]                                  # one step after C0, not C0 itself
    assert (first >= dumux1d.C0_DILUTE).all() and (first <= res["pore_C"][-1]).all()
    top_down = np.argsort(-res["pore_z"])
    assert (np.diff(res["pore_C"][-1][top_down]) < 0).all()    # enriched towards the surface
    assert (res["pore_C"][-1] > dumux1d.C0_DILUTE).all()


def test_the_coupled_column_follows_the_dumux_dilute_run(check_run):
    """The whole chain -- voxel vapour solve -> pore sink -> pressure -> salt transport -> a_w --
    against their dilute 20 C run.  Until 2026-09-19 this failed at 1.5 %: the pore's liquid-mass
    sink carried a factor rho_pore / RHO_L, 4.6 % short at this column's C = 0.20."""
    top_down = np.argsort(-check_run["pore_z"])
    ours = check_run["pore_C"][-1][top_down]
    theirs = np.array(dumux1d.DUMUX_DILUTE_STEADY)
    # measured 2.72e-3 (2026-09-21, conserving network step, the case's T = 293.15 K declared);
    # the pore network alone under their boundary condition: 3.02e-3.  The bound is the old one.
    assert np.abs(ours / theirs - 1.0).max() < 3.5e-3


def test_a_box_size_without_a_dgf_is_refused():
    from saltv2.coupling import driver3d
    with pytest.raises(ValueError, match="substrate_L_box"):
        driver3d.run(substrate_L_box=1e-3, rh=0.07, delta_eff_mm=18.0, verbose=False)


def test_the_order_of_the_pores_in_the_dgf_does_not_matter(tmp_path):
    """The same column written top pore FIRST.  In the shipped file the surface pore is the last
    one, and until 2026-09-19 the driver used a cluster's OPENING numbers as network indices, so
    on this substrate the structure drew its water and salt from the bottom pore; on a lattice,
    whose surface pores come first, the two numberings coincide and nothing showed."""
    from saltv2.coupling import driver3d
    from pathlib import Path
    lines = Path(dumux1d.DGF).read_text().splitlines()
    first = next(k for k, line in enumerate(lines) if line[:1].isdigit())
    n = lines.index("#") - first
    vertices = lines[first:first + n][::-1]
    start = next(k for k in range(first + n, len(lines)) if lines[k][:1].isdigit())
    throats = []
    for line in lines[start:lines.index("#", start)]:
        a, b, *rest = line.split()
        throats.append(" ".join([str(n - 1 - int(a)), str(n - 1 - int(b)), *rest]))
    flipped = tmp_path / "top_first.dgf"
    flipped.write_text("\n".join(lines[:first] + vertices + lines[first + n:start] + throats + ["#"]) + "\n")

    res = [driver3d.run(**dict(dumux1d.RUN_3D, substrate_dgf=path), verbose=False)
           for path in (dumux1d.DGF, flipped)]
    assert res[1]["layer"].surface[0] == 0 and res[0]["layer"].surface[0] == res[0]["layer"].N - 1
    by_height = [np.argsort(r["pore_z"]) for r in res]
    assert res[0]["t_nuc"][0] == pytest.approx(res[1]["t_nuc"][0], rel=1e-9)
    assert len(res[0]["hist"]) == len(res[1]["hist"])
    assert np.allclose(res[0]["pore_C"][:, by_height[0]], res[1]["pore_C"][:, by_height[1]],
                       rtol=1e-9, atol=0.0)
    for name in ("h_max", "V_wet", "E_total"):
        assert np.allclose(driver3d.col(res[0]["hist"], name), driver3d.col(res[1]["hist"], name),
                           rtol=1e-9, atol=1e-30), name


# ---- a pore under salt holds C_sat; what arrives beyond it is the structure's salt (2026-09-19) ---

@pytest.fixture(scope="module")
def undersaturated_growth():
    """RUN_3D: the reservoir (0.254) is BELOW saturation, so the structure can only be fed by what
    the column carries up against its own back-diffusion."""
    from saltv2.coupling import driver3d
    return driver3d.run(**dumux1d.RUN_3D, verbose=False)


def test_a_pore_under_salt_stays_at_saturation_and_the_column_above_its_reservoir(undersaturated_growth):
    """The 2014 JFM / DuMuX picture: brine in contact with crystal sits at C_sat.  Two earlier
    accountings failed here -- every face charging beta(C_sat) drained the column 7.6 % BELOW
    its reservoir, and an upwind hand-over let the covered pore relax to the reservoir's 0.254."""
    from saltv2.coupling import driver3d
    from saltv2.materials.brine import C_SAT
    res = undersaturated_growth
    t = driver3d.col(res["hist"], "t_h") * 3600.0
    after = t > res["t_nuc"][0]
    assert after.sum() > 100
    top_down = np.argsort(-res["pore_z"])
    C = res["pore_C"][after][:, top_down]
    assert np.allclose(C[:, 0], C_SAT, rtol=1e-12)                    # held, not approached
    assert C.min() > dumux1d.C0_SATURATED                             # nothing under the reservoir
    late = C[len(C) // 2:]                                            # the pre-nucleation hump has decayed
    assert (np.diff(late, axis=1) < 0).all()                          # C_sat at the top, falling to C0


def test_the_structure_gets_what_crystallises_and_not_a_gram_more(undersaturated_growth):
    """Every kilogram the structure holds left a pore's brine by crystallising: crystallised =
    supplied + seed + buffer.  Per kg of water the faces hand over at most beta(C_sat) -- Prat
    (2026) eq (29), the upper bound of this closure -- and here LESS, because a column fed from
    an undersaturated reservoir cannot deliver it."""
    from saltv2.materials.brine import C_SAT
    res = undersaturated_growth; led = res["ledger"]
    assert led["salt_crystallised"] > 0.0
    assert led["relative_residual_salt_interface"] == pytest.approx(0.0, abs=1e-10)
    assert led["salt_crystallised"] == pytest.approx(
        led["salt_supplied"] + led["salt_nucleated"] + float(res["interface_buffer"].sum()), rel=1e-12)
    assert led["salt_pnm_cap"] == 0.0 and led["salt_unsourced"] == 0.0     # nothing vanishes
    per_water = led["salt_released_by_evaporation"] / led["water_evaporated_on_salt"]
    assert 0.0 < per_water < C_SAT / (1 - C_SAT)
    assert led["relative_residual_salt_precipitation"] == pytest.approx(0.0, abs=1e-9)
    assert led["relative_residual_salt_network"] == pytest.approx(0.0, abs=1e-9)
