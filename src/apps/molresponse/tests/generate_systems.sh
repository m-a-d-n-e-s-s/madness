#!/usr/bin/env bash
#
# generate_systems.sh — Generate ground-state checkpoints for all fixture systems.
#
# Usage:
#   ./generate_systems.sh [--force]
#
# Environment:
#   MOLDFT_BIN    Path to moldft binary (default: moldft)
#   FIXTURE_DIR   Path to fixtures directory (default: ./fixtures)
#   MPI_NP        Number of MPI processes (default: 1)
#
# Each system directory under fixtures/systems/ must contain a moldft.in.
# This script runs moldft on each input and stores the checkpoint in the
# system directory. Systems with an existing moldft.restartdata are skipped
# unless --force is passed.

set -euo pipefail

MOLDFT_BIN="${MOLDFT_BIN:-moldft}"
FIXTURE_DIR="${FIXTURE_DIR:-$(dirname "$0")/fixtures}"
MPI_NP="${MPI_NP:-1}"
FORCE=false

if [[ "${1:-}" == "--force" ]]; then
    FORCE=true
fi

SYSTEMS_DIR="${FIXTURE_DIR}/systems"

if [[ ! -d "$SYSTEMS_DIR" ]]; then
    echo "ERROR: Systems directory not found: $SYSTEMS_DIR"
    exit 1
fi

echo "============================================"
echo "  molresponse_v3 fixture generator"
echo "============================================"
echo "MOLDFT_BIN:  $MOLDFT_BIN"
echo "FIXTURE_DIR: $FIXTURE_DIR"
echo "MPI_NP:      $MPI_NP"
echo "FORCE:       $FORCE"
echo ""

generated=0
skipped=0
failed=0

for system_dir in "$SYSTEMS_DIR"/*/; do
    system_name=$(basename "$system_dir")
    input_file="${system_dir}/moldft.in"
    restart_file="${system_dir}/moldft.restartdata.00000"

    if [[ ! -f "$input_file" ]]; then
        echo "SKIP  $system_name  (no moldft.in)"
        skipped=$((skipped + 1))
        continue
    fi

    if [[ -f "$restart_file" ]] && [[ "$FORCE" == "false" ]]; then
        echo "SKIP  $system_name  (checkpoint exists, use --force to regenerate)"
        skipped=$((skipped + 1))
        continue
    fi

    echo "RUN   $system_name ..."

    # Run moldft in the system directory so output files land there
    pushd "$system_dir" > /dev/null

    if [[ "$MPI_NP" -gt 1 ]]; then
        run_cmd="mpirun -np $MPI_NP $MOLDFT_BIN --input=moldft.in"
    else
        run_cmd="$MOLDFT_BIN --input=moldft.in"
    fi

    if $run_cmd > moldft.stdout 2> moldft.stderr; then
        echo "DONE  $system_name"
        generated=$((generated + 1))
    else
        echo "FAIL  $system_name  (exit code $?, see ${system_dir}moldft.stderr)"
        failed=$((failed + 1))
    fi

    popd > /dev/null
done

echo ""
echo "============================================"
echo "  Summary"
echo "============================================"
echo "Generated: $generated"
echo "Skipped:   $skipped"
echo "Failed:    $failed"
echo ""

if [[ "$failed" -gt 0 ]]; then
    exit 1
fi
