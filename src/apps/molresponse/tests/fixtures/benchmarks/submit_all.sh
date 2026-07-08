#!/bin/bash
# Submit all benchmark SLURM jobs. By default submits h2o, c2h4, c6h6 in order.
# Pick a subset by argument: ./submit_all.sh h2o c2h4
#
# Records submitted job ids in jobs.txt next to this script.

set -euo pipefail

DIR=$(cd "$(dirname "$0")" && pwd)
SYSTEMS=("$@")
[[ ${#SYSTEMS[@]} -eq 0 ]] && SYSTEMS=(h2o c2h4 c6h6)

JOBS_FILE=$DIR/jobs.txt
: > "$JOBS_FILE"

for sys in "${SYSTEMS[@]}"; do
  script=$DIR/$sys/benchmark.slurm
  if [[ ! -f $script ]]; then
    echo "SKIP $sys — no benchmark.slurm at $script" >&2
    continue
  fi
  cd "$DIR/$sys"
  out=$(sbatch benchmark.slurm)
  jobid=$(awk '{print $NF}' <<<"$out")
  echo "$sys $jobid" | tee -a "$JOBS_FILE"
done

echo
echo "submitted jobs recorded in $JOBS_FILE"
echo "check status:  squeue -u \$USER"
echo "tail logs:     tail -f \$WORK_BASE/<sys>_<jobid>/run_legacy_tda/run.log"
