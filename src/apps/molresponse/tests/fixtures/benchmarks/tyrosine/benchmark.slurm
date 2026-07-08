#!/bin/bash
#SBATCH --job-name=es_tyrosine_bench
#SBATCH --output=es_tyrosine_bench-%j.out
#SBATCH --error=es_tyrosine_bench-%j.err
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --partition=hbm-medium-96core

# Excited-state benchmark: tyrosine (auto-generated from h2o template).
# ntasks-per-node held at 4 to fit Xeon Max memory budget for large
# n_occ. See CLAUDE.md "Empirical data" table.

set -euo pipefail

REPO=${REPO:-$HOME/Projects/madness}
BUILD=${BUILD:-$REPO/build}
WORK_BASE=${WORK_BASE:-/gpfs/scratch/$USER/madness_es_bench}

SYSTEM=tyrosine
FIXTURE=$REPO/src/apps/molresponse/tests/fixtures/systems/${SYSTEM}_hf
BENCH=$REPO/src/apps/molresponse/tests/fixtures/benchmarks/$SYSTEM

BIN_MOLDFT=$BUILD/src/apps/moldft/moldft
BIN_LEGACY=$BUILD/src/apps/molresponse_legacy/molresponse_legacy
BIN_V1=$BUILD/src/apps/molresponse/molresponse

for f in "$BIN_MOLDFT" "$BIN_LEGACY" "$BIN_V1"; do
  [[ -x $f ]] || { echo "ERROR: missing binary $f" >&2; exit 1; }
done
[[ -d $FIXTURE ]] || { echo "ERROR: fixture dir not found: $FIXTURE" >&2; exit 1; }

export MAD_NUM_THREADS=20
source ~/load_xeonmax.sh

LAUNCHER=${LAUNCHER:-mpirun --map-by numa numactl --preferred-many=8-15}

WORK=$WORK_BASE/${SYSTEM}_${SLURM_JOB_ID:-local}
mkdir -p "$WORK" && cd "$WORK"
echo "WORK=$WORK"

cp "$FIXTURE/moldft.in"  ./moldft.in
cp "$FIXTURE/geometry.xyz" ./geometry.xyz 2>/dev/null || true

echo "=== ground state (moldft) ==="
time $LAUNCHER "$BIN_MOLDFT" moldft.in > moldft.log 2>&1
ls mad.* 2>&1 | head

TIMING_CSV=$WORK/timing.csv
echo "system,code,mode,wall_seconds,iters,omega1,omega2,omega3,omega4,exit" > "$TIMING_CSV"

run_one() {
  local code=$1 mode=$2 bin=$3 input_name=$4
  local rundir=$WORK/run_${code}_${mode}
  mkdir -p "$rundir" && cd "$rundir"
  ln -sf "$WORK/mad.restartdata.00000" .
  ln -sf "$WORK/mad.restartaodata"     . 2>/dev/null || true
  ln -sf "$WORK/mad.calc_info.json"    . 2>/dev/null || true
  cp "$BENCH/response_${mode}.in" "./$input_name"
  local log=$rundir/run.log
  local t0=$(date +%s.%N)
  set +e; $LAUNCHER "$bin" > "$log" 2>&1; local rc=$?; set -e
  local t1=$(date +%s.%N)
  local wall=$(awk "BEGIN{printf \"%.2f\", $t1 - $t0}")
  local iters=$(grep -cE "^[[:space:]]*Iteration [0-9]+" "$log" 2>/dev/null || echo 0)
  local final=$(awk '/Final excitation energies:/{getline; print}' "$log" | tail -1 | tr -s ' ' | sed 's/^ *//')
  local w1 w2 w3 w4
  read -r _ w1 w2 w3 w4 _ <<<"$final"
  : "${w1:=NA}" "${w2:=NA}" "${w3:=NA}" "${w4:=NA}"
  echo "$SYSTEM,$code,$mode,$wall,$iters,$w1,$w2,$w3,$w4,$rc" >> "$TIMING_CSV"
  echo "  -> $code/$mode: ${wall}s, $iters iters, exit $rc"
  cd "$WORK"
}

echo; echo "=== legacy TDA ==="; run_one legacy tda "$BIN_LEGACY" input
echo; echo "=== legacy RPA ==="; run_one legacy rpa "$BIN_LEGACY" input
echo; echo "=== v1 TDA     ==="; run_one v1     tda "$BIN_V1"     response.in
echo; echo "=== v1 RPA     ==="; run_one v1     rpa "$BIN_V1"     response.in

echo; echo "=== summary ==="
column -t -s, "$TIMING_CSV"
