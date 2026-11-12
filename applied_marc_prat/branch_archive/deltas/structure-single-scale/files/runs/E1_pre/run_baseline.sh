#!/bin/bash
# The self-consistent baseline on frodo: one systemd user unit per (case, salt structure, vpp).
# Run from the LOCAL machine, from the worktree this script lives in: syncs THAT tree (never the
# logs or results) to its OWN remote directory -- not the one of the main checkout -- and starts
# the units.
#   bash runs/E1_pre/run_baseline.sh                                  # the default table below: THIS branch's structure, single_scale
#   RUNS="1mm:single_scale:7 dumux1d:hierarchical:12" bash runs/E1_pre/run_baseline.sh
#   bash runs/E1_pre/run_baseline.sh status
#   bash runs/E1_pre/run_baseline.sh fetch                            # records and figures back
# The tree is synced without .git, so the commit hash travels in CODE_VERSION (as in run_crust.sh).
set -e
LOCAL=$(cd "$(dirname "$0")/../.." && pwd)
REMOTE=frodo
D=/home/zhixin/Documents/Humboldt/salt/$(basename "$LOCAL")
PY=/home/zhixin/Documents/Humboldt/salt/salt_grow_v2/.venv/bin/python
LOGS=$D/runs/E1_pre/logs_baseline
RUNS="${RUNS:-5-50um:single_scale:7 50-62um:single_scale:7 100-125um:single_scale:7 1mm:single_scale:7 dumux1d:single_scale:12 1mm:single_scale:14 100-125um:single_scale:14}"

if [ "$1" = "status" ]; then
  ssh $REMOTE "systemctl --user list-units 'saltv2_base_*' --no-pager --all | head -20; \
               for f in $LOGS/*.log; do echo \"--- \$f\"; tail -n 2 \$f; done"
  exit 0
fi
if [ "$1" = "fetch" ]; then
  rsync -az --include 'baseline_*/' --include 'baseline_*/*' --exclude '*' \
        "$REMOTE:$D/runs/E1_pre/simulation_outputs/" "$LOCAL/runs/E1_pre/simulation_outputs/"
  exit 0
fi

git -C "$LOCAL" describe --always --dirty > "$LOCAL/runs/E1_pre/CODE_VERSION"
rsync -az --exclude 'logs*' --exclude '*.npz' --exclude '.git' --exclude '.venv' \
      --exclude simulation_outputs "$LOCAL/" "$REMOTE:$D/"
for run in $RUNS; do
  IFS=: read -r cls st vpp <<< "$run"
  unit="saltv2_base_${cls//-/_}_${st}_vpp${vpp}"
  ssh $REMOTE "mkdir -p $LOGS && systemd-run --user --unit=$unit --collect \
      -p WorkingDirectory=$D \
      /usr/bin/bash -c \"OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        exec $PY -u runs/E1_pre/baseline.py run $cls $st --vpp $vpp > $LOGS/${cls}_${st}_vpp${vpp}.log 2>&1\""
done
ssh $REMOTE "systemctl --user list-units 'saltv2_base_*' --no-pager"
