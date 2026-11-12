"""Stage-T launcher: pre-crystallisation substrate transport, paired B2 / C networks.

Cases live in ``saltv2.cases.stage_t_wicking.CASES``; this script only runs them and writes
records.  Regime: evaporation-wicking, saturated network, prescribed surface water loss.

    python runs/T1_substrate_wicking/stage_t.py list
    python runs/T1_substrate_wicking/stage_t.py manifest --out DIR      # provenance (host, code)
    python runs/T1_substrate_wicking/stage_t.py run C_hotspot B2_hotspot --out DIR
    python runs/T1_substrate_wicking/stage_t.py run-all --jobs 8 --out DIR
    python runs/T1_substrate_wicking/stage_t.py report --out DIR        # DIR/summary.md

Each case writes DIR/<name>.json (summary, tracked) and DIR/<name>.npz (fields, git-ignored).
Compute runs on frodo (run_frodo.sh); locally only ``list`` / ``report`` / tiny checks.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from saltv2.cases import stage_t_wicking as st  # noqa: E402

HASHED = ["saltv2/pnm/conservative.py", "saltv2/pnm/flow.py", "saltv2/pnm/transport.py",
          "saltv2/pnm/network.py", "saltv2/pnm/grids/N5x5_R270um_SD10um_L3.5mmx3.5mm_n1_top.dgf",
          "saltv2/materials/brine.py", "saltv2/materials/constants.py",
          "saltv2/cases/stage_t_wicking.py", "saltv2/cases/vtp2014_column.py",
          "runs/T1_substrate_wicking/stage_t.py"]


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance():
    import scipy
    prov = dict(host=socket.gethostname(), platform=platform.platform(),
                python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                executable=sys.executable, cwd=str(ROOT),
                time_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                omp_num_threads=os.environ.get("OMP_NUM_THREADS"),
                files={f: _sha256(ROOT / f) for f in HASHED})
    p = ROOT / "PROVENANCE.json"                 # written by run_frodo.sh into the remote tree
    if p.exists():
        prov["code"] = json.loads(p.read_text())
    else:
        try:
            prov["code"] = dict(head=subprocess.check_output(
                ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
                dirty=bool(subprocess.check_output(
                    ["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip()))
        except Exception as e:  # noqa: BLE001
            prov["code"] = dict(error=str(e))
    return prov


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    return x


def run_one(name, out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    try:
        summary, arrays = st.run(**st.CASES[name])
    except Exception:  # noqa: BLE001  -- a failed case is recorded, never silently dropped
        (out / f"{name}.FAILED.txt").write_text(traceback.format_exc())
        return name, "FAILED"
    summary["name"] = name
    summary["case"] = st.CASES[name]
    summary["provenance"] = provenance()
    np.savez_compressed(out / f"{name}.npz", **arrays)
    (out / f"{name}.json").write_text(json.dumps(_jsonable(summary), indent=1))
    return name, f"t_sat {summary['t_sat_first_h']} h, wall {summary['wall_s']:.0f} s"


def cmd_report(out):
    out = Path(out)
    rows = {p.stem: json.loads(p.read_text()) for p in sorted(out.glob("*.json"))
            if p.stem not in ("manifest", "PROVENANCE")}
    lines = ["| case | mode | dt [s] | layers | n | t_sat first [h] | first pore | "
             "var w_surf @1h | var @t_sat | sub-lat liquid/E @1h | hot salt from below/lateral @1h | "
             "rel res salt | rel res water | p_req G1 max [Pa] | margin open G1 | wall [s] |",
             "|" + "---|" * 16]
    for name, s in rows.items():
        s1 = s["snapshots"].get("t1h", {}); sf = s["snapshots"].get("first_saturation", {})
        t = "censored" if s["censored"] else f"{s['t_sat_first_h']:.4f}"
        hot = (f"{s1.get('hot_salt_from_below', float('nan')):.3e} / "
               f"{s1.get('hot_salt_lateral', float('nan')):.3e}") if s1 else "-"
        lines.append(
            f"| {name} | {s['mode']} | {s['dt']:g} | {s['n_layers']} | {s['n_side']} | {t} | "
            f"{s['first_pore']} | {s1.get('w_surf_var', float('nan')):.3e} | "
            f"{sf.get('w_surf_var', float('nan')):.3e} | "
            f"{s1.get('subsurface_lateral_liquid_over_E', float('nan')):.3e} | {hot} | "
            f"{s['rel_res_salt']:+.2e} | {s['rel_res_water']:+.2e} | "
            f"{sf.get('p_req_G1_max', float('nan')):.1f} | {sf.get('margin_open_G1', float('nan')):.3f} | "
            f"{s['wall_s']:.0f} |")
    failed = sorted(p.name for p in out.glob("*.FAILED.txt"))
    if failed:
        lines += ["", "FAILED: " + ", ".join(failed)]
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["list", "manifest", "run", "run-all", "report"])
    ap.add_argument("names", nargs="*")
    ap.add_argument("--out", type=Path, default=ROOT / "runs/T1_substrate_wicking/results/local")
    ap.add_argument("--jobs", type=int, default=1)
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for k, v in st.CASES.items():
            print(f"{k:28s} {v}")
        return
    if a.cmd == "manifest":
        a.out.mkdir(parents=True, exist_ok=True)
        m = dict(provenance=provenance(), cases=st.CASES)
        (a.out / "manifest.json").write_text(json.dumps(_jsonable(m), indent=1))
        print(json.dumps(m["provenance"], indent=1))
        return
    if a.cmd == "report":
        cmd_report(a.out)
        return
    names = list(st.CASES) if a.cmd == "run-all" and not a.names else a.names
    unknown = [n for n in names if n not in st.CASES]
    if unknown:
        raise SystemExit(f"unknown cases: {unknown}")
    # largest first, so the long lateral-domain runs do not start last
    names.sort(key=lambda n: -st.CASES[n]["n_side"] ** 2 * st.CASES[n]["n_layers"] / st.CASES[n]["dt"])
    if a.jobs <= 1:
        for n in names:
            print(*run_one(n, a.out), flush=True)
    else:
        with ProcessPoolExecutor(a.jobs) as ex:
            for name, msg in ex.map(run_one, names, [a.out] * len(names)):
                print(name, msg, flush=True)
    cmd_report(a.out)


if __name__ == "__main__":
    main()
