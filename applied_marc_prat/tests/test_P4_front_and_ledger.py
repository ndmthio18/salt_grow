"""P4 acceptance tests of plan/17: front densification, the analytic sealing
time, and the water/salt mass ledger.

The forbidden practices of P4 are tested as well: front deposition may not be
accumulated into one global scalar (it must land in named voxels), and a
non-finite or unnamed contribution must not be able to enter the ledger
silently.
"""
import numpy as np
import pytest

from saltv2 import (SALT_SEALED, SALT_WET, BirthClosure, CrustField,
                    Grid, GeometryEngine, MassLedger)

A = 1e-4
RHO_CR = 2170.0            # NaCl crystal density [kg/m^3]
C_SAT = 0.264              # NaCl saturation mass fraction
D_SUB = 68e-6


def field_with_wet_layer(shape=(4, 4, 10), h=3, **kw):
    g = Grid(shape, A)
    eng = GeometryEngine(g)
    eng.state[..., :h] = SALT_WET
    cf = CrustField(g, eng, BirthClosure(d_sub=D_SUB, **kw), rho_cr=RHO_CR)
    return g, eng, cf


def tau_seal(eps_s0, eps_min, a, F, j):
    """plan/17 eq. of section 4.4:
    tau = (eps_s0 - eps_min) a rho_cr (1 - C_eq) / (F j C_eq)."""
    return (eps_s0 - eps_min) * a * RHO_CR * (1.0 - C_SAT) / (F * j * C_SAT)


# ------------------------------------------------------- front densification
def test_deposition_lands_in_named_voxels_not_in_a_global_scalar():
    g, eng, cf = field_with_wet_layer()
    dm = np.zeros(g.shape)
    dm[1, 1, 2] = 1e-10
    dm[2, 3, 2] = 3e-10
    cf.deposit(dm)
    assert cf.m_precipitated[1, 1, 2] == pytest.approx(1e-10)
    assert cf.m_precipitated[2, 3, 2] == pytest.approx(3e-10)
    assert cf.m_precipitated.sum() == pytest.approx(4e-10)
    # porosity of exactly those two salt voxels dropped, of no other
    salt = eng.solid_mask()
    dropped = salt & (cf.eps_s < cf.closure.eps_s0 - 1e-15)
    assert dropped.sum() == 2 and dropped[1, 1, 2] and dropped[2, 3, 2]
    assert cf.eps_s[1, 1, 2] == pytest.approx(
        cf.closure.eps_s0 - 1e-10 / (RHO_CR * g.cell_volume))


def test_deposition_outside_salt_or_negative_is_rejected():
    g, eng, cf = field_with_wet_layer()
    bad = np.zeros(g.shape); bad[0, 0, 8] = 1e-9        # a gas cell
    with pytest.raises(ValueError, match="salt voxels"):
        cf.deposit(bad)
    neg = np.zeros(g.shape); neg[0, 0, 2] = -1e-9
    with pytest.raises(ValueError, match="non-negative"):
        cf.deposit(neg)
    too_much = np.zeros(g.shape); too_much[0, 0, 2] = 1e-8   # more than the pore volume
    with pytest.raises(ValueError, match="overfills"):
        cf.deposit(too_much)


def test_numerical_sealing_time_matches_the_analytic_tau_seal():
    """Constant front flux j on one face of one voxel: the voxel must seal at
    exactly tau_seal, and dt_to_seal must never let it overshoot eps_min."""
    eps_min = 0.05
    F, j = 0.6, 2.0e-5                                  # [-], [kg m^-2 s^-1]
    g, eng, cf = field_with_wet_layer(shape=(3, 3, 6), h=2)
    eps0 = cf.closure.eps_s0
    # salt mass rate precipitated in the front voxel: F j C/(1-C) per unit area
    rate = np.zeros(g.shape)
    rate[1, 1, 1] = F * j * C_SAT / (1.0 - C_SAT) * g.face_area
    expected = tau_seal(eps0, eps_min, g.a, F, j)
    t = 0.0
    for _ in range(500):
        dt = min(cf.dt_to_seal(rate, eps_min), 0.1 * expected)
        cf.deposit(rate * dt)
        t += dt
        assert cf.eps_s[1, 1, 1] >= eps_min - 1e-14     # never steps past eps_min
        if cf.seal(eps_min, t).any():
            break
    assert eng.state[1, 1, 1] == SALT_SEALED
    assert t == pytest.approx(expected, rel=1e-10)
    assert cf.t_sealed[1, 1, 1] == pytest.approx(expected, rel=1e-10)
    # and the salt that did the sealing is exactly the pore volume it filled
    assert cf.m_precipitated[1, 1, 1] == pytest.approx(
        (eps0 - eps_min) * RHO_CR * g.cell_volume, rel=1e-10)


def test_dt_to_seal_is_infinite_when_nothing_densifies():
    g, eng, cf = field_with_wet_layer()
    assert not np.isfinite(cf.dt_to_seal(np.zeros(g.shape), 0.05))


def test_salt_mass_is_continuous_across_sealing():
    """The state transition moves solid between states; it must not create or
    destroy salt (plan/17 P4 analytic test 4)."""
    g, eng, cf = field_with_wet_layer(shape=(4, 4, 8), h=3)
    m0 = cf.salt_solid_mass()
    dm = np.zeros(g.shape); dm[:, :, 2] = 2e-10
    cf.deposit(dm)
    m1 = cf.salt_solid_mass()
    assert m1 - m0 == pytest.approx(dm.sum(), rel=1e-12)
    cf.eps_s[:, :, 2] = 0.02
    cf.seal(0.02, 20.0)
    assert eng.state[0, 0, 2] == SALT_SEALED
    # sealing itself changed no mass: only eps_s (set by hand here) did
    assert cf.salt_solid_mass() == pytest.approx(
        float(np.sum(1.0 - cf.eps_s[eng.solid_mask()])) * g.cell_volume * RHO_CR, rel=1e-14)


def test_growth_and_front_salt_add_up_to_the_structure_mass():
    """Growth births voxels at eps_s0, the front fills their pores: the two
    channels together must equal what the structure holds (plan/17 P4)."""
    g = Grid((5, 5, 12), A)
    eng = GeometryEngine(g)
    eng.state[:, :, 0] = SALT_WET
    cf = CrustField(g, eng, BirthClosure(d_sub=D_SUB), rho_cr=RHO_CR)
    led = MassLedger()
    # the balance the driver actually checks: what the structure holds must
    # equal what the two source channels put into it
    led.register("salt", plus=["salt_supplied"],
                 minus=["salt_grown", "salt_front"], rtol=1e-12)
    led.add("salt_grown", cf.salt_solid_mass())            # the seed layer
    rng = np.random.default_rng(5)
    for _ in range(25):
        Vn = 1e-8 * (0.3 + rng.random(g.shape))
        rate = eng.volume_rate(Vn)
        if rate.max() <= 0.0:
            break
        dt = eng.dt_occupancy(rate, f_max=0.5)
        rep = eng.advance(rate, dt)
        if rep.flipped.any():
            cf.birth(rep.flipped, eng.t)
            led.add("salt_grown", rep.flipped.sum() * g.cell_volume
                    * (1.0 - cf.closure.eps_s0) * RHO_CR)
        wet = eng.solid_mask()
        dm = np.where(wet, 1e-12, 0.0)
        cf.deposit(dm)
        led.add("salt_front", dm.sum())
    # occupancy that has not yet flipped is not part of the structure mass, so
    # compare against the voxels that exist
    assert led["salt_grown"] + led["salt_front"] == pytest.approx(
        cf.salt_solid_mass(), rel=1e-12)
    led.add("salt_supplied", cf.salt_solid_mass())
    led.check("salt")


# ---------------------------------------------------------------- the ledger
def test_ledger_rejects_unknown_channels_and_non_finite_values():
    led = MassLedger()
    with pytest.raises(KeyError, match="unknown ledger channel"):
        led.add("salt_frnot", 1.0)
    with pytest.raises(ValueError, match="non-finite"):
        led.add("salt_front", np.nan)
    with pytest.raises(KeyError):
        led.register("bad", plus=["nope"], minus=[])


def test_ledger_balance_detects_a_leak():
    led = MassLedger()
    led.register("water", plus=["salt_supplied"],
                 minus=["water_evaporated", "water_evaporated_base"], rtol=1e-10)
    led.add("salt_supplied", 1.0); led.add("water_evaporated", 0.6); led.add("water_evaporated_base", 0.4)
    assert led.relative_residual("water") == pytest.approx(0.0, abs=1e-15)
    led.check("water")
    led.add("water_evaporated", 1e-6)                   # a leak of 1 ppm
    assert led.relative_residual("water") == pytest.approx(-1e-6, rel=1e-6)
    with pytest.raises(RuntimeError, match="mass balance 'water' violated"):
        led.check("water")


def test_ledger_report_carries_every_channel_and_residual():
    led = MassLedger()
    led.register("salt", plus=["salt_supplied"],
                 minus=["salt_grown", "salt_front", "salt_pnm_cap"])
    for c, v in (("salt_supplied", 2.0), ("salt_grown", 1.0), ("salt_front", 0.5), ("salt_pnm_cap", 0.5)):
        led.add(c, v)
    rep = led.report()
    for c in ("salt_grown", "salt_front", "salt_pnm_cap", "salt_supplied",
              "water_evaporated", "water_evaporated_base", "volume_unplaced"):
        assert c in rep
    assert rep["residual_salt"] == pytest.approx(0.0, abs=1e-15)
    assert "relative_residual_salt" in rep


def test_a_voxel_seals_in_one_dt_to_seal_step_for_any_rate():
    """Regression (2026-09-05): dt_to_seal cuts dt so that the limiting voxel lands on eps_min,
    but deposit() recomputes the drop from rate * dt and round-off can leave eps one ulp above
    eps_min.  Without a tolerance seal() misses it and the next dt_to_seal is a zero-length step,
    repeated forever (the 2-D verification case stalled at 26 000 zero steps).  Every rate must
    seal the voxel in exactly one step."""
    from saltv2 import Grid, GeometryEngine, BirthClosure, CrustField, SALT_WET, SALT_SEALED
    g = Grid((3, 4), 4.41e-6); eng = GeometryEngine(g); eng.seed(1, height_cells=2)
    crust = CrustField(g, eng, BirthClosure(d_sub=34e-6, c_k=3e-3), rho_cr=2170.0)
    eps_min = 0.02
    rng = np.random.default_rng(0)
    for _ in range(500):
        crust.eps_s[1, 1] = 0.43; eng.state[1, 1] = SALT_WET
        crust.t_sealed[1, 1] = np.nan
        dm_rate = np.zeros(g.shape); dm_rate[1, 1] = 10 ** rng.uniform(-16, -12)
        dt = crust.dt_to_seal(dm_rate, eps_min)
        assert np.isfinite(dt) and dt > 0.0
        crust.deposit(dm_rate * dt); crust.seal(eps_min, 0.0)
        assert eng.state[1, 1] == SALT_SEALED
        assert crust.dt_to_seal(dm_rate, eps_min) == np.inf
