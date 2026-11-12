"""The production path is bit for bit what it was before the 2026-09-14 cleanup.

Regime: evaporation-wicking.  Each test re-runs one of the cases in ``reference_cases.CASES``
and compares its history, ledger, per-voxel state and salt masses with the record captured from
the checkpoint tree 95c1acb (``tests/data/ref_*.npz``).  Exact equality, no tolerance: the
cleanup deletes code, it does not touch a single number.  If one of these fails after an edit,
the edit changed behaviour -- fix the edit, never regenerate the reference.  The only exception
is a behaviour change the user decided on, re-captured on purpose and listed in
``reference_cases.py`` (so far: ref_3d_small, square pore openings, 2026-09-15).

References are platform-pinned (numpy/scipy/platform in the file) to FRODO since 2026-09-15:
ref_2d_small and ref_2d_cv were re-captured there from a0db925 (the checkpoint 95c1acb plus this
harness, no model change), ref_3d_small from f0ff577 (square openings).  Run the suite on frodo.
On another build (the local machine) the small cases differ at the ulp level, which then grows
into a few shifted events over ~1000 steps; that is a platform question, settled by
regenerating on that platform from those same commits, never by loosening the assert.

The medium case (``scancrust 7``: 68 um, 2 x 2 pores, vpp 6, 60 h, ~3000 steps) lives on frodo
in ``runs/E1_pre/crust_out/`` and is compared there with ``python -m tests.npz_compare``; it is
gated behind ``SALTV2_MEDIUM=1`` because it takes tens of minutes.
"""
import os
from pathlib import Path

import numpy as np
import pytest

from tests.npz_compare import compare_run
from tests.reference_cases import CASES, DATA, call, capture

SMALL = [n for n in CASES]


def _fresh_as_npz(kind, kwargs, tmp_path):
    rec = capture(kind, call(kind, kwargs))
    p = tmp_path / "fresh.npz"
    np.savez(p, **rec)
    return np.load(p, allow_pickle=True)


@pytest.mark.parametrize("name", SMALL)
def test_small_case_matches_reference(name, tmp_path, run_3d_small):
    kind, kwargs = CASES[name]
    ref = np.load(DATA / f"{name}.npz", allow_pickle=True)
    if name == "ref_3d_small":
        # the shared session fixture IS this case; make sure the two definitions never drift
        from tests.conftest import CASE_3D
        import inspect
        from saltv2.coupling import driver3d as v
        accepted = inspect.signature(v.run).parameters
        assert {k: val for k, val in kwargs.items() if k in accepted} == CASE_3D
        rec = capture(kind, run_3d_small)
        np.savez(tmp_path / "fresh.npz", **rec)
        new = np.load(tmp_path / "fresh.npz", allow_pickle=True)
    else:
        new = _fresh_as_npz(kind, kwargs, tmp_path)
    bad = [b for b in compare_run(ref, new) if not b.startswith("note:")]
    assert not bad, "\n".join(bad)
    # every quantity the reference carries must still be produced (a silently dropped column
    # or channel is a change too); columns removed ON PURPOSE are listed here by name
    removed_columns = {"E_vent", "f_vent_open", "m_front", "V_dry_open", "h_max_wet", "n_rewet",
                       "h_front"}
    removed_channels = {"salt_front_kinetic", "salt_base_unplaced", "water_supplied",
                        "water_in_pores"}
    if "hist_columns" in ref.files:
        gone = set(map(str, ref["hist_columns"])) - set(map(str, new["hist_columns"]))
        assert gone <= removed_columns, f"unexpected columns dropped: {gone - removed_columns}"
    gone = set(ref["ledger_names"].tolist()) - set(new["ledger_names"].tolist())
    assert gone <= removed_channels, f"unexpected ledger channels dropped: {gone - removed_channels}"


@pytest.mark.skipif(not os.environ.get("SALTV2_MEDIUM"), reason="set SALTV2_MEDIUM=1 on frodo")
def test_scancrust_7_matches_reference():
    """Reference: crust_out/ref_scancrust_7.npz (the checkpoint tree's own re-run of case 7);
    fresh: crust_out/vapor_crust_7.npz written by ``pre_experiment_vapor.py scancrust 7 snap``
    from the current tree."""
    here = Path(__file__).resolve().parents[1] / "runs" / "E1_pre" / "crust_out"
    ref = np.load(here / "ref_scancrust_7.npz", allow_pickle=True)
    new = np.load(here / "vapor_crust_7.npz", allow_pickle=True)
    bad = [b for b in compare_run(ref, new) if not b.startswith("note:")]
    assert not bad, "\n".join(bad)
