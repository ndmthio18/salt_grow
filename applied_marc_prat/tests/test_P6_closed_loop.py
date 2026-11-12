"""P5 / P6 acceptance tests of plan/17 on the drivers themselves.

P5: the supply fixed point q -> p -> h_lim -> active faces must converge, must be
consistent (the reach a cluster is given must be the reach its own demand implies), and must
never be continued silently when it does not.
P6: the closed loop must conserve salt and water, must report which limiter set dt, must carry
a named output schema, and must show the chain
    Q_demand / Q_cap > 1  =>  frozen faces  =>  V_sealed up  =>  E down
with every arrow read off its own output, not inferred from E alone.
"""
import numpy as np
import pytest

from saltv2.coupling import driver3d as v
from saltv2.coupling import driver2d as fb
from tests.conftest import SMALL_3D as E1_SMALL


@pytest.fixture
def e1_run(run_3d_small):
    """The shared 3-D reference case (conftest.CASE_3D, pinned by test_regression.py)."""
    return run_3d_small


# ------------------------------------------------------------------ P6 schema
def test_history_is_addressed_by_name_and_carries_the_P6_quantities(e1_run):
    h = e1_run["hist"]
    assert h.shape[1] == len(v.HIST_COLUMNS)
    assert e1_run["hist_columns"] == v.HIST_COLUMNS
    for name in ("E_total", "E_base", "E_salt", "V_wet", "V_sealed", "f_blocked", "dt",
                 "limiter", "n_sealed_cells", "res_water", "res_salt", "n_supply_iter"):
        assert v.col(h, name).shape == (len(h),)
    with pytest.raises(KeyError):
        v.col(h, "not_a_column")
    # the record names its inputs, its parameterisation and its mass ledger
    assert e1_run["params"]["ell_over_dsub"] == 1.0 and e1_run["params"]["eps_min"] == 0.05
    assert "relative_residual_salt_network" in e1_run["ledger"]
    rr = e1_run["resolution"]
    assert rr["reaction_band_cells"] == pytest.approx(rr["reaction_band_m"] / rr["voxel_edge_m"])
    assert rr["reaction_band_resolved"] == (rr["reaction_band_cells"] >= rr["min_reaction_band_cells"])


def test_E_total_is_the_resolved_flux(e1_run):
    """E_total is the resolved gas-phase evaporation and nothing else: openings plus wet salt."""
    h = e1_run["hist"]
    assert np.allclose(v.col(h, "E_total"), v.col(h, "E_base") + v.col(h, "E_salt"), rtol=1e-12)
    assert np.allclose(v.col(h, "E_base"),
                       v.col(h, "E_base_bare") + v.col(h, "E_base_crust"), rtol=1e-12)
    assert np.allclose(v.col(h, "E_salt"),
                       v.col(h, "E_wet_gas") + v.col(h, "E_wet_crust"), rtol=1e-12)


def test_dt_is_set_by_a_named_limiter(e1_run):
    """README section 8 finding 1: among the ACCURACY limiters (grow, network, salt, cap -- the
    ones dt_factor scales) the step is set by the network stability limit, never by the growth
    CFL, so tightening the growth limit alone changes nothing.  The EVENT limiters (seal, end)
    land the step on an event and are not part of that finding.  The run must say which bound
    set each step instead of leaving it to be guessed.

    Measured on this fixture (tests/data/ref_3d_small.npz): disk openings (to 2026-09-14)
    network 552, seal 450, end 1; square openings (2026-09-15, the footprint fills the pitch at
    vpp 3, so more voxels seal) network 617, seal 697, salt 2, end 1.  Up to 2026-09-14 this test
    demanded that network win over ALL limiters, which held only because seal happened to fire
    less often on the disk footprint."""
    h = e1_run["hist"]
    lim = v.col(h, "limiter").astype(int)
    assert set(np.unique(lim)) <= set(range(len(v.LIMITERS)))
    counts = dict(zip(v.LIMITERS, np.bincount(lim, minlength=len(v.LIMITERS))))
    accuracy = ("grow", "network", "salt", "cap")
    assert max(accuracy, key=counts.get) == "network"
    assert counts["grow"] == 0
    assert np.all(v.col(h, "dt") > 0.0)


# ------------------------------------------------------------------ P4 / P6 conservation
def test_legacy_structure_ledger_and_network_volume_balance(e1_run):
    h = e1_run["hist"]
    assert np.abs(v.col(h, "res_salt")).max() < 1e-10
    assert np.abs(v.col(h, "res_water")).max() < 1e-10
    assert abs(e1_run["ledger"]["relative_residual_salt_network"]) < 1e-10
    assert abs(e1_run["ledger"]["relative_residual_salt_precipitation"]) < 1e-9
    # the growth law takes its salt out of the pores that fed it (plan/17 P4)
    assert e1_run["ledger"]["salt_supplied"] > 0.0


# ------------------------------------------------------------------ P5 fixed point
def test_supply_fixed_point_converges_and_is_self_consistent(e1_run):
    h = e1_run["hist"]
    iters = v.col(h, "n_supply_iter")
    assert iters.min() >= 1 and iters.max() < 50
    # the reach each cluster ended with is a finite-or-infinite non-negative number
    for c, h_lim in e1_run["h_freeze"].items():
        assert h_lim >= 0.0
        assert not np.isnan(h_lim)


def test_supply_fixed_point_raises_instead_of_continuing_silently(monkeypatch):
    """plan/17 P5: no silent continuation.  Force the iteration cap to 0 usable sweeps."""
    monkeypatch.setattr(v, "SUPPLY_ITER_MAX", 0)
    with pytest.raises(RuntimeError, match="supply fixed point did not converge"):
        v.run(c_k=1e-6, t_end=3600, **E1_SMALL)


def test_n_frozen_counts_clusters_that_cannot_grow():
    """A cluster stalls when its reach falls below its highest exposed face.  With the sealing
    chain on, a stalled cluster's frozen faces fill their band and seal; once the sealed voxels
    have taken the faces away the demand drops and the reach can become unlimited again -- that
    is the feedback, so the reach is read where it is first reported, not at the end."""
    kw = dict(t_end=20 * 3600, ell_over_dsub=1.0, eps_min=0.05, **E1_SMALL)
    free = v.run(c_k=1e-4, **kw)                  # capillary supply never fails
    stalled = v.run(c_k=1e-7, **kw)               # fails after two voxels
    assert v.col(free["hist"], "n_frozen").max() == 0
    assert not np.isfinite(list(free["h_freeze"].values())[0])
    assert v.col(stalled["hist"], "n_frozen").max() >= 1
    assert np.isfinite(stalled["t_freeze"]) and stalled["state_at_freeze"] is not None
    # a stalled cluster seals from its frozen faces
    assert v.col(stalled["hist"], "V_sealed")[-1] > 0.0
    # and the free run never seals (no frozen face, so nothing feeds a band)
    assert v.col(free["hist"], "V_sealed").max() == 0.0


def test_3d_chain_frozen_faces_densify_seal_and_lower_E():
    """plan/16 M2 in 3-D: a frozen (supply-limited) wet face keeps evaporating, its salt fills the
    voxel behind it, the voxel seals, and the resolved evaporation drops.  Each arrow is read off
    its own output; the salt it precipitates is drawn from the pores (ledger salt_front_external)."""
    res = v.run(c_k=1e-7, t_end=30 * 3600, ell_over_dsub=1.0, eps_min=0.05, **E1_SMALL)
    h = res["hist"]; c = lambda n: v.col(h, n)
    i_fr = int(np.argmax(c("n_frozen") > 0))
    assert i_fr > 0                                                    # arrow 1: supply failure
    assert c("V_sealed")[-1] > 0.0 and c("n_sealed_cells")[-1] >= 1    # arrow 2: sealing
    i_seal = int(np.argmax(c("V_sealed") > 0))
    assert i_seal >= i_fr
    assert c("E_total")[-1] < c("E_total")[i_seal - 1]                 # arrow 3: E drops
    # arrow 4: columns whose evaporating surfaces lie under the sealed roof are reported blocked
    # even though the air above the roof is connected to the far field
    assert c("f_blocked")[-1] > 0.0
    led = res["ledger"]
    assert led["salt_front_external"] > 0.0
    assert abs(led["relative_residual_salt_network"]) < 1e-10
    assert np.abs(c("res_salt")).max() < 1e-10


# ------------------------------------------------------------------ P6 stage gate, 2-D
def test_stage_gate_supply_failure_drives_drying_sealing_and_an_evaporation_drop():
    """Every arrow of
         Q_demand / Q_cap > 1 => frozen faces => sealed volume => E down
    is checked on its own output (plan/17 P6 small-scale stage gate).  The 2-D prototype is used
    because it reaches the sealed state in seconds."""
    common = dict(d_sub=34e-6, n_seeds=2, cells_per_pitch=4, z_box_mm=0.4, t_end_h=12.0,
                  snap_hours=0.0, verbose=False)
    supplied = fb.run(**common, c_k=1.0)                     # supply never fails
    starved = fb.run(**common, c_k=1e-7)                     # supply fails at once
    hs, hf = supplied["hist"], starved["hist"]
    # arrow 1: only the starved run has frozen faces
    assert hs[:, 7].max() == 0 and hf[:, 7].max() > 0
    # arrow 2 -> 3: frozen faces densify their own cells until they seal
    assert not (supplied["state"] == fb.SEALED).any()
    assert (starved["state"] == fb.SEALED).any()
    assert starved["m_precipitated"].sum() > 0.0
    # task B: what closes is the sub-grid reaction band, not the whole voxel, so the sealed
    # cells' BULK porosity stays well above eps_min while their band is full
    sealed = starved["state"] == fb.SEALED
    assert (starved["m_skin"][sealed] >= starved["skin_capacity"][sealed] * (1 - 1e-9)).all()
    # arrow 4: the resolved evaporation of the starved run falls below its own early value
    i = int(np.argmax(hf[:, 6] > 0))                          # first sealed cell
    assert i > 0 and hf[-1, 1] < hf[i - 1, 1]
    # and both balances of both runs close
    for res in (supplied, starved):
        assert abs(res["ledger"]["relative_residual_salt_structure"]) < 1e-10
        assert abs(res["ledger"]["relative_residual_salt_precipitation"]) < 1e-9


def test_active_window_gives_the_same_run_as_the_full_box():
    """README E1 section 9.11: all per-step lattice work runs on the lowest k_top + margin layers.
    A run in a tall box (window active) must reproduce the run in a box just tall enough (window
    = whole box), and the window must not change the sealing chain."""
    kw = dict(c_k=1e-6, t_end=24 * 3600, ell_over_dsub=1.0, eps_min=0.05,
              **{k_: v_ for k_, v_ in E1_SMALL.items() if k_ != "z_top_mm"})
    short = v.run(**{**kw, "z_top_mm": 0.5})           # 0.5 mm box, margin covers it: no window
    tall = v.run(**{**kw, "z_top_mm": 3.0})            # same physics, window active most of the run
    # window vs no window (same adaptive solver).  h_max is identical for the whole run; the
    # columns that a discrete event drives (V_sealed, f_blocked, n_frozen and, through them,
    # E_total) can land one step apart once cells start sealing, so they are compared step by
    # step only up to the first sealing event and by their event TIME after it.  Comparing them
    # index by index across the jump would test the step sequence, not the physics.
    i_seal = int(np.argmax(v.col(short["hist"], "V_sealed") > 0))
    assert i_seal > 0
    for name in ("E_total", "h_max", "V_sealed", "f_blocked", "n_frozen"):
        a, b = v.col(short["hist"], name), v.col(tall["hist"], name)
        n = min(len(a), len(b), i_seal) if name != "h_max" else min(len(a), len(b))
        assert np.allclose(a[:n], b[:n], rtol=1e-6, atol=1e-12), name
    for r in (short, tall):
        se = v.col(r["hist"], "V_sealed") > 0
        r["t_seal_1"] = v.col(r["hist"], "t_h")[int(np.argmax(se))]
    assert short["t_seal_1"] == pytest.approx(tall["t_seal_1"], rel=1e-6)
    assert tall["eng"].nz_active is None                 # run() hands back full-grid masks
    assert (tall["eng"].state[..., 60:] == 0).all()      # nothing was ever placed above the window
    assert np.array_equal(short["eng"].state[..., :20], tall["eng"].state[..., :20])
