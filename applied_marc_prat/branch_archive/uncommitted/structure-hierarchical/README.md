# Uncommitted `structure-hierarchical` worktree overlay

- Source worktree:
  `/home/zhixin/Documents/Humboldt/salt/salt_grow_v2-hierarchical`
- Base branch tip: `6faf96f41d78d0cfbf97e7a40dc234de3dc7a38a`
- Captured: 2026-11-12
- Status: four modified tracked files and one untracked test.

Files captured:

- `runs/E1_pre/baseline.py`
- `runs/E1_pre/run_baseline.sh`
- `saltv2/coupling/driver3d.py`
- `saltv2/growth/crust.py`
- `tests/test_layer_av.py` (untracked on liqin)

`files/` contains the complete captured versions; `uncommitted.patch` contains
the tracked diff; `worktree-status.txt` records the original status; and
`file-hashes.csv` records local hashes.

The overlay implements a layer-dependent `Da*`/`F` route with clamping and a
dense-at-birth/sealing transition. These choices repair one inconsistency in
the committed hierarchical route but introduce additional unvalidated
closures. They are development material, not a registered branch result.

Restoring the bundle, checking out `structure-hierarchical`, and overlaying
these five files gave `2 passed in 3.60 s` for `tests/test_layer_av.py` on the
current Windows/Python environment. This verifies capture and the dedicated
test only; it is not physical validation or a full-suite qualification.

