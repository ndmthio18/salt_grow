# liqin branch archive

This directory preserves and organises all 14 local Git branches found in the
liqin `salt_grow_v2` repository on 2026-11-12.

The archive has two complementary layers:

1. `salt_grow_v2_all_branches.bundle` is the authoritative, byte-preserving
   Git history. It contains all branch refs and their complete reachable
   history.
2. `deltas/` contains browsable exports of the seven branch tips that are not
   ancestors of `verify/step2a-prat-single-scale`. Each export contains only
   the files changed on that route relative to its merge base with Step 2A.

`uncommitted/structure-hierarchical/` separately preserves the dirty worktree
found on liqin. It is not part of any Git branch or registered result.

Read `BRANCH_INDEX.md` for the scientific role and limitations of each route.

## Bundle integrity

- File: `salt_grow_v2_all_branches.bundle`
- Size: 48,921,084 bytes
- SHA-256:
  `0F8F44BC80C921269112920793B515D0011145333B8947640B42C87E98C3CA88`
- `git bundle verify`: passed; Git reported a complete history with 14 heads.

To restore every branch exactly:

```powershell
git clone --mirror salt_grow_v2_all_branches.bundle salt_grow_v2.git
git --git-dir=salt_grow_v2.git branch --list
```

For an ordinary working checkout:

```powershell
git clone salt_grow_v2_all_branches.bundle salt_grow_v2_restored
cd salt_grow_v2_restored
git branch -a
git switch --track origin/review/prat-beta-hierarchy
```

## Delta layout

Every `deltas/<route>/` directory contains:

- `files/`: files added or modified on that route;
- `commits.tsv`: commits unique relative to Step 2A;
- `name-status.tsv`: Git path/status inventory;
- `diffstat.txt`: compact change summary;
- `changes.patch`: text diff, with binary changes identified but not embedded;
- `file-hashes.csv`: branch blob IDs, local SHA-256 values, and sizes.

Windows Git exported browsable text files with CRLF line endings. The manifest
therefore records both raw local blob IDs and Git-normalised blob IDs. All
normalised blob IDs match their branch-tip blobs; binary files are unchanged.
The bundle remains the authoritative source for exact original bytes.

Historical `CLAUDE.md` files inside deltas are archived project records, not
active instructions for this repository.

