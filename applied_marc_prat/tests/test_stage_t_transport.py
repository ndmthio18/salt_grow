"""Stage-T: the opt-in conservative transport step and the paired B2 / C networks.

Regime: evaporation-wicking, before crystallisation (saturated rigid network on a reservoir,
prescribed surface water loss, no salt leaves the network).

What is pinned here:
1. ``pnm.conservative_step`` closes the liquid, salt and water books to round-off, on a network
   with lateral redistribution (hotspot demand on a 3-D lattice).
2. The production pair ``solve_pressure`` + ``salt_step`` is audited with the same yardstick and
   its residual IS the missing storage V w d rho (salt) / V d rho (liquid) -- the documented
   limitation, measured exactly, not hidden and not fixed.
3. B2 and C share every pore and every throat except the subsurface lateral ones.
4. Null control: under uniform demand on a uniform network B2 = C = one isolated column.
5. With a hotspot only C carries subsurface lateral flow.
6. At steady state (no storage) the conservative and the production paths give the same profile.
"""
import numpy as np
import pytest

from saltv2.cases import stage_t_wicking as st
from saltv2.materials import brine
from saltv2.materials.constants import G
from saltv2.pnm import (SUBSURFACE_LATERAL, SURFACE_LATERAL, VERTICAL, PoreLayer,
                        conservative_step, hydrostatic_head, salt_step, solve_pressure,
                        step_residuals, throat_kind)

C0 = st.C0


def _small(variant, n_side=3, n_layers=4):
    return st.build_network(variant, n_side=n_side, n_layers=n_layers)


def _march(layer, q, n_steps, dt, mode):
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    w = np.full(layer.N, C0); m = brine.salt_mass(w, layer.V)
    tot = dict(salt=0.0, water=0.0, liquid=0.0, salt_in=0.0, out=0.0, storage_salt=0.0,
               storage_liquid=0.0)
    for _ in range(n_steps):
        w_old = w
        if mode == "conservative":
            r = conservative_step(layer, w, C0, q, g_lat, g_res, d_lat, d_res, dt)
            w, rho, Ql, Qr = r["w"], r["rho_flux"], r["Q_lat"], r["Q_res"]
        else:
            rho = brine.density(w)
            _, Ql, Qr = solve_pressure(layer, q, g_lat, g_res, rho, float(brine.density(C0)))
            m = salt_step(layer, m, C0, Ql, Qr, d_lat, d_res, dt)
            w = brine.mass_fraction_from_salt_mass(m, layer.V)
        res = step_residuals(layer, w_old, w, rho, Ql, Qr, d_lat, d_res, C0, q, dt)
        for k in ("salt", "water", "liquid", "out"):
            tot[k] += res[k]
        tot["salt_in"] += res["salt_in"]
        drho = brine.density(w) - brine.density(w_old)
        tot["storage_salt"] += float(np.sum(layer.V * w * drho))
        tot["storage_liquid"] += float(np.sum(layer.V * drho))
    return w, tot


def test_conservative_step_closes_liquid_salt_and_water_to_round_off():
    layer = _small("C")
    q = st.demand(layer, "hotspot")
    w, tot = _march(layer, q, 120, 60.0, "conservative")
    assert w[layer.surface].max() > C0 + 5e-3                  # it did enrich
    assert abs(tot["salt"]) <= 1e-10 * tot["salt_in"]
    assert abs(tot["water"]) <= 1e-10 * tot["out"]
    assert abs(tot["liquid"]) <= 1e-10 * tot["out"]


def test_production_path_residual_is_exactly_the_missing_storage():
    """Documents the known limitation of the production path; it must NOT be rebaselined away."""
    layer = _small("C")
    q = st.demand(layer, "hotspot")
    _, tot = _march(layer, q, 120, 60.0, "production")
    assert tot["salt"] == pytest.approx(tot["storage_salt"], rel=1e-6)
    assert tot["liquid"] == pytest.approx(tot["storage_liquid"], rel=1e-6)
    assert tot["liquid"] > 1e-3 * tot["out"]                  # not a round-off effect


def test_B2_and_C_differ_only_in_the_subsurface_lateral_throats():
    C, B2, B1 = _small("C"), _small("B2"), _small("B1")
    for attr in ("x", "y", "z", "r", "V", "rt_res", "Lt_res", "surface", "bottom"):
        assert np.array_equal(getattr(C, attr), getattr(B2, attr)), attr
    kC, kB2 = throat_kind(C), throat_kind(B2)
    for kind in (SURFACE_LATERAL, VERTICAL):
        assert [t for t, k in zip(C.lat, kC) if k == kind] == \
               [t for t, k in zip(B2.lat, kB2) if k == kind]
    assert np.sum(kC == SUBSURFACE_LATERAL) == 3 * 2 * 9 and np.sum(kB2 == SUBSURFACE_LATERAL) == 0
    assert np.all(throat_kind(B1) == VERTICAL)
    assert np.array_equal(st.demand(C, "hotspot"), st.demand(B2, "hotspot"))


@pytest.mark.parametrize("pattern", ["uniform", "stripe", "hotspot"])
def test_every_pattern_imposes_the_same_total_demand(pattern):
    layer = _small("C", n_side=5)
    assert np.sum(st.demand(layer, pattern)) == pytest.approx(25 * st.J_MEAN * st.PITCH ** 2,
                                                              rel=1e-13)


def test_null_control_uniform_B2_equals_C_equals_one_column():
    kw = dict(pattern="uniform", n_layers=6, t_max_h=0.5)
    sC, aC = st.run("C", n_side=3, **kw)
    sB2, aB2 = st.run("B2", n_side=3, **kw)
    s1, a1 = st.run("C", n_side=1, **kw)
    wC, wB2, w1 = (a["snap_end_w"] for a in (aC, aB2, a1))
    assert np.max(np.abs(wC - wB2)) < 1e-13
    assert np.max(np.abs(wC.reshape(6, 9) - w1[:, None])) < 1e-13
    assert sC["snapshots"]["end"]["subsurface_lateral_liquid_over_E"] < 1e-9


def test_hotspot_drives_subsurface_lateral_flow_in_C_only():
    kw = dict(pattern="hotspot", n_side=3, n_layers=6, t_max_h=0.25)
    sC, _ = st.run("C", **kw)
    sB2, _ = st.run("B2", **kw)
    assert sC["snapshots"]["end"]["subsurface_lateral_liquid_over_E"] > 1e-3
    assert sB2["snapshots"]["end"]["subsurface_lateral_liquid_over_E"] == 0.0
    for s in (sC, sB2):
        assert abs(s["rel_res_salt"]) < 1e-9 and abs(s["rel_res_water"]) < 1e-9


def test_conservative_and_production_paths_share_the_steady_profile():
    """Storage vanishes at steady state, so the two paths must agree there."""
    layer = PoreLayer(1, n_side=1, pore_radius=10e-6, n_layers=8, dz=200e-6, Lt_res=200e-6)
    g_lat, g_res = layer.conductances(); d_lat, d_res = layer.diffusive()
    q = np.zeros(layer.N); q[layer.surface] = 1.5e-16 * 1190.0
    w = np.full(layer.N, C0); m = brine.salt_mass(w, layer.V)
    for _ in range(80):
        w = conservative_step(layer, w, C0, q, g_lat, g_res, d_lat, d_res, 1e5)["w"]
        rho = brine.density(brine.mass_fraction_from_salt_mass(m, layer.V))
        _, Ql, Qr = solve_pressure(layer, q, g_lat, g_res, rho, float(brine.density(C0)))
        m = salt_step(layer, m, C0, Ql, Qr, d_lat, d_res, 1e5)
    w_prod = brine.mass_fraction_from_salt_mass(m, layer.V)
    assert w[0] > C0 * 1.05
    assert np.max(np.abs(w - w_prod)) < 1e-10


def test_hydrostatic_head_integrates_the_column_density():
    C = _small("C", n_layers=4); A = _small("A", n_layers=4)
    w = np.full(C.N, 0.2)
    rho = float(brine.density(0.2))
    H = 4 * st.DZ
    assert np.allclose(hydrostatic_head(C, w, H), G * rho * H, rtol=1e-13)
    assert np.allclose(hydrostatic_head(C, w, 2.5 * st.DZ), G * rho * 2.5 * st.DZ, rtol=1e-13)
    assert np.allclose(hydrostatic_head(A, np.full(A.N, 0.2), H), G * rho * H, rtol=1e-13)


def test_fixed_background_hotspot_is_the_hotspot_on_5x5_and_keeps_the_local_field():
    five = _small("C", n_side=5)
    assert np.allclose(st.demand(five, "hotspot_fixedbg"), st.demand(five, "hotspot"), rtol=1e-14)
    ten = _small("C", n_side=10)
    q5 = st.demand(five, "hotspot_fixedbg")[five.surface]
    q10 = st.demand(ten, "hotspot_fixedbg")[ten.surface]
    assert q10.max() == pytest.approx(q5.max(), rel=1e-14)
    assert np.min(q10) == pytest.approx(np.min(q5), rel=1e-14)
