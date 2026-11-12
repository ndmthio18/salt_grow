# Source and verification provenance

## Primary snapshot

- Host alias: `liqin`
- Remote path:
  `/home/zhixin/Documents/Humboldt/salt/salt_grow_v2-step2a`
- Branch: `verify/step2a-prat-single-scale`
- Commit: `f7dbf3876df9a563c3572844d48afd886395614f`
- Commit date: 2026-10-26 22:46:10 +0100
- Commit subject: `Step 2A: Prat (2026) single-scale growth law against the independent 1-D reaction-transport reference (verify/step2a)`

The copied paths are `prat2026/`, `saltv2/`, `tests/`,
`runs/S_surface/`, `knowledge/`, and `ideas/`. They were exported with
`git archive` from the exact commit above. The archive SHA-256 was:

`9CE240721567AA7C204BBD27BD9602FDE3A818922D5B0A0543CB3F5813E99976`

After extraction, the local relative-path inventory matched the remote Git
tree exactly: 264 files locally, 264 files remotely, zero missing and zero
extra files.

## Later applicability review

- Remote path:
  `/home/zhixin/Documents/Humboldt/salt/salt_grow_v2/ideas/Prat生长律适用性审查_2026-10-25.md`
- Introducing commit:
  `b7419a24e2151467a3f4757e06cceb22df939237`
- Git blob: `dbad5882ccb3946e2b01415b69273d33e4699d76`
- Local and remote SHA-256:
  `EFC4DA5FC69B8D062F644350DFA1C96508D9DECC2330E8F835079E93BFE05282`

This review was copied separately because the Step 2A verification branch
diverged before it was added, while the review directly evaluates the same
`saltv2` production path.

## Environments inspected on 2026-11-12

Liqin's shared historical virtual environment reported:

- Python 3.13.12
- NumPy 2.5.1
- SciPy 1.18.0
- Matplotlib 3.11.0
- Pillow 12.3.0
- pytest 9.1.1

The local verification environment reported:

- Python 3.12.4 on Windows
- NumPy 2.0.0
- SciPy 1.14.0
- Matplotlib 3.9.2
- Pillow 10.4.0
- pytest 9.1.1

The historical repository does not include a complete lock file, so the liqin
package versions above should not be assumed to be the exact versions that
created every stored reference array.

## Local verification

### Core migrated material

With `PYTHONUTF8=1`, the analytical package, Step 2A test, and architecture
tests gave:

`48 passed in 2.58 s`

Command:

```powershell
python -m pytest -q tests/test_layering.py tests/test_step2a_prat.py prat2026/tests
```

### Complete suite

The unmodified complete suite gave:

`269 passed, 1 skipped, 6 failed in 790.60 s`

Three failures in `tests/test_layering.py` were caused by Windows' default GBK
decoding of UTF-8 source files because the test uses `Path.read_text()` without
an explicit encoding. Those tests pass when `PYTHONUTF8=1` is set.

Three failures in `tests/test_regression.py` were caused by its intentional
bit-for-bit equality rule. The discrete state arrays remained exactly equal.
The largest observed absolute differences were approximately:

- `ref_3d_small`: `1.73e-9` in history/time-like fields, `1.53e-14` in
  porosity, and `1.29e-26` in precipitated mass;
- `ref_2d_small`: `3.33e-16` in history and `8.27e-25` in precipitated mass;
- `ref_2d_cv`: `8.67e-13` in history, `6.19e-14` in porosity, and
  `1.21e-21` in precipitated mass.

These results show that the imported source runs in the current environment
and that the core Prat verification passes. They do **not** establish portable
bitwise reproducibility of the historical full-driver baselines, and they are
not evidence of experimental or physical validation.

## Complete branch archive

All 14 local branch refs in the liqin repository were exported to
`branch_archive/salt_grow_v2_all_branches.bundle`. `git bundle verify` passed
and reported complete history. Bundle SHA-256:

`0F8F44BC80C921269112920793B515D0011145333B8947640B42C87E98C3CA88`

Seven branch tips that are separate from the Step 2A ancestry were also
exported as browsable deltas. Their path inventories match Git exactly, and
every Git-normalised local blob matches the corresponding branch-tip blob.
The only dirty worktree found was `structure-hierarchical`; its four tracked
modifications and one untracked test are preserved separately and are not
represented as a commit.

