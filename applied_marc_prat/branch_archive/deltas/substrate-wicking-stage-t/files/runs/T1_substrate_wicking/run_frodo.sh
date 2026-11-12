#!/bin/bash
# Stage-T on frodo, in an ISOLATED directory (never the production tree on frodo).
# Run from the LOCAL repo root.  The code shipped is the working tree as a content-addressed git
# tree (git write-tree through a temporary index: tracked + untracked-not-ignored files, no commit
# made), so the record carries the baseline HEAD, the tree hash and whether the tree was dirty.
#
#   bash runs/T1_substrate_wicking/run_frodo.sh launch [TAG] [JOBS] [CASES...]  # default: all cases
#   bash runs/T1_substrate_wicking/run_frodo.sh ship-test [TAG]                  # ship the tree, full pytest, no cases
#   bash runs/T1_substrate_wicking/run_frodo.sh test TAG                         # full pytest there
#   bash runs/T1_substrate_wicking/run_frodo.sh status TAG
#   bash runs/T1_substrate_wicking/run_frodo.sh fetch TAG                        # -> results/TAG/
set -euo pipefail
LOCAL=$(git rev-parse --show-toplevel)
REMOTE=frodo
BASE=/home/zhixin/saltv2_stage_t
PY=/home/zhixin/Documents/Humboldt/salt/salt_grow_v2/.venv/bin/python   # interpreter only (read)
cmd=${1:-launch}; TAG=${2:-}

if [ "$cmd" = "launch" ] || [ "$cmd" = "ship-test" ]; then
  HEAD=$(git -C "$LOCAL" rev-parse HEAD)
  IDX=$(mktemp); trap 'rm -f "$IDX"' EXIT
  GIT_INDEX_FILE=$IDX git -C "$LOCAL" read-tree HEAD
  GIT_INDEX_FILE=$IDX git -C "$LOCAL" add -A
  TREE=$(GIT_INDEX_FILE=$IDX git -C "$LOCAL" write-tree)
  DIRTY=$([ "$TREE" = "$(git -C "$LOCAL" rev-parse HEAD^{tree})" ] && echo false || echo true)
  TAG=${TAG:-$(date -u +%Y%m%dT%H%M%SZ)_${TREE:0:10}}
  JOBS=${3:-8}; shift 3 2>/dev/null || shift $#
  CASES="$*"
  RD=$BASE/$TAG
  ssh $REMOTE "test ! -e $RD" || { echo "$RD exists on $REMOTE; choose another TAG" >&2; exit 1; }
  git -C "$LOCAL" archive "$TREE" | ssh $REMOTE "mkdir -p $RD/code $RD/results $RD/logs && tar -x -C $RD/code"
  ssh $REMOTE "cat > $RD/code/PROVENANCE.json" <<EOF
{"baseline_head": "$HEAD", "tree": "$TREE", "dirty_vs_head": $DIRTY,
 "branch": "$(git -C "$LOCAL" branch --show-current)", "shipped_utc": "$(date -u +%FT%TZ)",
 "local_host": "$(hostname)", "remote_dir": "$RD"}
EOF
  if [ "$cmd" = "ship-test" ]; then
    ssh $REMOTE "cd $RD/code && OMP_NUM_THREADS=1 timeout 3600 $PY -m pytest tests/ -q -p no:cacheprovider > $RD/logs/pytest.log 2>&1; tail -n 15 $RD/logs/pytest.log"
    echo "tested tree $TREE (head $HEAD, dirty $DIRTY) in $REMOTE:$RD"
    exit 0
  fi
  unit="saltv2_staget_${TAG//[^A-Za-z0-9]/_}"
  ssh $REMOTE "cd $RD/code && OMP_NUM_THREADS=1 $PY runs/T1_substrate_wicking/stage_t.py manifest --out $RD/results > $RD/logs/manifest.log 2>&1 && \
      systemd-run --user --unit=$unit --collect -p WorkingDirectory=$RD/code \
      /usr/bin/bash -c \"OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        exec $PY -u runs/T1_substrate_wicking/stage_t.py run-all $CASES --jobs $JOBS --out $RD/results > $RD/logs/run.log 2>&1\""
  echo "launched $unit in $REMOTE:$RD (tree $TREE, head $HEAD, dirty $DIRTY)"
  exit 0
fi
[ -n "$TAG" ] || { echo "TAG required" >&2; exit 1; }
RD=$BASE/$TAG
case "$cmd" in
  test)
    ssh $REMOTE "cd $RD/code && OMP_NUM_THREADS=1 timeout 3600 $PY -m pytest tests/ -q -p no:cacheprovider > $RD/logs/pytest.log 2>&1; tail -n 15 $RD/logs/pytest.log" ;;
  status)
    ssh $REMOTE "systemctl --user list-units 'saltv2_staget_*' --no-pager --all | head; ls $RD/results | wc -l; tail -n 5 $RD/logs/run.log" ;;
  fetch)
    mkdir -p "$LOCAL/runs/T1_substrate_wicking/results/$TAG"
    rsync -az "$REMOTE:$RD/results/" "$REMOTE:$RD/logs" "$LOCAL/runs/T1_substrate_wicking/results/$TAG/"
    rsync -az "$REMOTE:$RD/code/PROVENANCE.json" "$LOCAL/runs/T1_substrate_wicking/results/$TAG/" ;;
  *) echo "unknown command $cmd" >&2; exit 1 ;;
esac
