#!/usr/bin/env bash
# baseline/verify.sh — re-run the Phase 0.2 baseline cases and diff against
# the recorded outputs in baseline/{poisson2d_serial,poisson2d_mpi2}/.
#
# Exit 0 if all outputs reproduce; non-zero otherwise. Used as the gate
# after each step of the library port.
#
# Assumes Exasim has already been built into ./build/ (text2code,
# cput2cEXASIM, cpumpit2cEXASIM).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
APP="$ROOT/apps/poisson/poisson2d"
BUILD="$ROOT/build"
BASE="$ROOT/baseline"

fail=0

run_case () {
    local label="$1"; shift
    rm -rf "$APP/dataout"
    "$BUILD/text2code" "$APP/pdeapp.txt" >/dev/null
    ( cd "$APP" && "$@" >/dev/null )

    diff -q "$APP/dataout/outqoi.txt" "$BASE/$label/outqoi.txt" \
        || { echo "[FAIL] $label outqoi.txt differs"; fail=1; }

    local got; got="$(cd "$APP/dataout" && md5 -r *.bin | sort)"
    local want; want="$(sort "$BASE/$label/md5.txt")"
    if [ "$got" != "$want" ]; then
        echo "[FAIL] $label dataout md5s differ"
        diff <(echo "$want") <(echo "$got") || true
        fail=1
    else
        echo "[ OK ] $label"
    fi
}

run_case poisson2d_serial "$BUILD/cput2cEXASIM" "$APP/pdeapp.txt"
run_case poisson2d_mpi2   mpirun -np 2 "$BUILD/cpumpit2cEXASIM" "$APP/pdeapp.txt"

exit "$fail"
