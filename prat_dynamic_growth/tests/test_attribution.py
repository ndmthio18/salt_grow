import json
import numpy as np
import pytest
from prat_growth.parameters import PhysicalParameters, NumericalParameters
from prat_growth.dynamic import DynamicALEModel, snapshot_diagnostics
from prat_growth.frozen import solve_frozen
from prat_growth.analytic import general_growth_velocity
from prat_growth.verification import local_audit
from prat_growth.provenance import identity, cache_matches
from prat_growth.campaign import save_result


def test_q0_converges_to_paper():
    p = PhysicalParameters()
    errors = []
    for count in (160, 320, 640):
        n = NumericalParameters(n_cells=count, grid_power=3)
        q = solve_frozen(p, n, p.h0, np.full(count, p.epsilon_ref), "Q0")
        errors.append(abs(q.velocity/general_growth_velocity(p, p.h0)-1))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < 6e-5


def test_frozen_retains_pore_rate_and_matches_dynamic_concentration_operator():
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=100, grid_power=3)
    m = DynamicALEModel(p, n, "B2")
    y = m.initial_state()
    y[:100] = .3+.13*np.exp(-(1-m.centers)/.001)
    q = solve_frozen(p, n, p.h0, y[:100], "Q6")
    y[100:200] = q.concentration
    dy = m.rhs(0, y)
    assert np.max(abs(q.pore_rate)) > 0
    # Cancellation in absolute mass fractions limits the thinnest cell.
    assert np.max(abs(dy[100:200])) < 2e-7
    audit = local_audit(m, y)
    assert max(audit.values()) < 1e-9
    damaged = dy.copy()
    damaged[130] += .01
    assert local_audit(m, y, damaged)["salt_linf"] > 1e-4


def test_bad_top_state_fails_closed():
    m = DynamicALEModel(PhysicalParameters(), NumericalParameters(n_cells=16), "B2")
    with pytest.raises(RuntimeError, match="no physical top root"):
        m.top_state(-1000., .43, m.p.D_eff_ref, m.p.h0)


def test_cache_checks_configuration_code_and_data(tmp_path):
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=16, t_end=1, output_interval=1)
    r = DynamicALEModel(p, n, "B2").run("quasi_steady")
    save_result(tmp_path, "test", r)
    meta, data = tmp_path/"runs/test.json", tmp_path/"runs/test.npz"
    expected = identity(p, n, "B2", "quasi_steady")
    assert cache_matches(meta, data, expected)
    assert not cache_matches(meta, data, identity(p.with_updates(beta=3), n, "B2", "quasi_steady"))
    assert not cache_matches(meta, data, {**expected, "code_sha256": "changed"})
    data.write_bytes(data.read_bytes()+b"damage")
    assert not cache_matches(meta, data, expected)


@pytest.mark.parametrize("closure", ["B1", "B2"])
def test_quasisteady_initial_salt_inventory(closure):
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=80, grid_power=3, t_end=20, output_interval=10, max_step=1)
    r = DynamicALEModel(p, n, closure).run("quasi_steady")
    assert abs(snapshot_diagnostics(r, 0)["salt_residual_normalized"]) < 1e-14
    assert abs(snapshot_diagnostics(r, -1)["salt_residual_normalized"]) < 1e-8


def test_uniform_structure_controls_collapse():
    p = PhysicalParameters()
    n = NumericalParameters(n_cells=80, grid_power=3)
    eps = np.full(n.n_cells, p.epsilon_ref)
    controls = [solve_frozen(p, n, p.h0, eps, q) for q in ("Q2", "Q3", "Q4", "Q5", "Q6")]
    np.testing.assert_allclose([q.velocity for q in controls], controls[0].velocity, rtol=1e-12)


def test_archived_m3_convergence_and_attribution():
    import csv
    from hashlib import sha256
    from pathlib import Path
    out = Path(__file__).resolve().parents[1]/"outputs/m3"
    if not (out/"manifest.json").exists():
        pytest.skip("run python -m prat_growth.attribution to generate evidence")
    def rows(name):
        with (out/name).open() as stream:
            return list(csv.DictReader(stream))
    convergence = {row["case"]: row for row in rows("convergence.csv")}
    for key in ("velocity", "h", "epsilon_min", "L95_interpolated", "internal_fraction"):
        assert abs(float(convergence["grid_640"][key])/float(convergence["grid_480"][key])-1) < .005
    assert int(convergence["grid_640"]["cells_delta"]) >= 8
    for case in ("step_30", "step_120"):
        assert abs(float(convergence[case]["velocity"])/float(convergence["grid_640"]["velocity"])-1) < 1e-6
    for row in rows("local_audits.csv"):
        for key in ("liquid_linf", "salt_linf", "solid_linf"):
            assert float(row[key]) < 1e-8
        assert float(row["interface_salt_relative"]) < 1e-10
    contributions = rows("contributions.csv")
    matched = rows("matched.csv")
    for t in (300., 3600., 10800., 21600.):
        parts = [r for r in contributions if r["n_cells"] == "640" and float(r["time_s"]) == t]
        transient = next(r for r in matched if r["case"] == "T2" and float(r["time_s"]) == t)
        delta = float(transient["velocity"])-general_growth_velocity(PhysicalParameters(), .001)
        assert abs(sum(float(r["delta_velocity"]) for r in parts)-delta) < 1e-20
    manifest = json.loads((out/"manifest.json").read_text())
    for path in (out/"runs").glob("*.json"):
        meta = json.loads(path.read_text())
        expected = identity(PhysicalParameters(**meta["physical"]), NumericalParameters(**meta["numerical"]),
                            meta["closure"], meta["initial_profile"])
        # M4 adds source files, so M3 is an intact archive rather than a current
        # reusable cache. Do not rewrite its metadata to force a hash match.
        assert meta["identity"]["config_sha256"] == expected["config_sha256"]
        assert meta["identity"]["code_sha256"] == manifest["code_sha256"]
        assert meta["data_sha256"] == sha256(path.with_suffix(".npz").read_bytes()).hexdigest()
