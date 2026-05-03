#!/usr/bin/env bash
# Run the four supported tutorial paths against poisson2d for one
# backend variant. Verify each produces a converged solution.
#
# Usage:  bash tutorial/run-all.sh [--variant cpu|gpu|mpi|mpi_gpu]
#                                  [--build <dir>] [--np N]
#
# Defaults:  --variant cpu   --build ./build   --np 2 (mpi variants)

set -uo pipefail

EXASIM="$(cd "$(dirname "$0")/.." && pwd)"
VARIANT="cpu"
BUILD=""
NP="2"

while [ $# -gt 0 ]; do
    case "$1" in
        --variant) VARIANT="$2"; shift 2 ;;
        --build)   BUILD="$2";   shift 2 ;;
        --np)      NP="$2";      shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

if [ -z "$BUILD" ]; then
    case "$VARIANT" in
        cpu)     BUILD="$EXASIM/build" ;;
        gpu)     BUILD="$EXASIM/build_gpu" ;;
        mpi)     BUILD="$EXASIM/build_mpi" ;;
        mpi_gpu) BUILD="$EXASIM/build_mpi_gpu" ;;
        *) echo "unknown variant: $VARIANT" >&2; exit 2 ;;
    esac
fi

case "$VARIANT" in
    cpu)     SUFFIX="";          RUNNER=();                 IS_MPI=0 ;;
    gpu)     SUFFIX="_gpu";      RUNNER=();                 IS_MPI=0 ;;
    mpi)     SUFFIX="_mpi";      RUNNER=(mpirun -np "$NP"); IS_MPI=1 ;;
    mpi_gpu) SUFFIX="_mpi_gpu";  RUNNER=(mpirun -np "$NP"); IS_MPI=1 ;;
esac

if [ -d "$BUILD" ]; then
    BUILD="$(cd "$BUILD" && pwd)"
fi

APP="$EXASIM/apps/library_example"

echo "[tutorial] VARIANT = $VARIANT (suffix=$SUFFIX)"
echo "[tutorial] BUILD   = $BUILD"
[ "$IS_MPI" = 1 ] && echo "[tutorial] NP      = $NP"

echo "[tutorial] regenerating poisson2d codegen artifacts ..."
bash "$APP/regenerate.sh" poisson2d > /dev/null

PASSED=()
FAILED=()
SKIPPED=()
record_pass() { PASSED+=("$1"); echo "[ OK ]   $1"; }
record_fail() { FAILED+=("$1"); echo "[FAIL]   $1"; }
record_skip() { SKIPPED+=("$1"); echo "[SKIP]   $1"; }

run_cmd() {
    # Wrapper: run a command, capture exit. Returns 0/non-zero.
    "$@"
}

# ---- Path 01 ---- hand-written + run.hpp
P1_NAME="01: hand-written + run.hpp"
P1_BIN="$BUILD/poisson2d_template${SUFFIX}"
if [ -x "$P1_BIN" ]; then
    cd "$APP/poisson2d"
    rm -rf datain dataout
    mkdir -p datain dataout
    if [ "$IS_MPI" = 1 ]; then
        sed "s/^mpiprocs *=.*/mpiprocs = $NP;/" pdeapp.txt > pdeapp_run.txt
    else
        sed 's/^mpiprocs *=.*/mpiprocs = 1;/' pdeapp.txt > pdeapp_run.txt
    fi
    if "${RUNNER[@]}" "$P1_BIN" ./pdeapp_run.txt > /tmp/tut_01.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$P1_NAME"
    else
        record_fail "$P1_NAME ($P1_BIN)"
    fi
else
    record_skip "$P1_NAME ($P1_BIN not built)"
fi

# ---- Path 02 ---- hand-written + facade
P2_NAME="02: hand-written + ExasimSolver<M>"
P2_BIN="$BUILD/solve_square_facade${SUFFIX}"
if [ -x "$P2_BIN" ]; then
    cd "$EXASIM"
    if EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$P2_BIN" > /tmp/tut_02.log 2>&1 \
        && grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_02.log; then
        record_pass "$P2_NAME"
    else
        record_fail "$P2_NAME ($P2_BIN)"
    fi
else
    record_skip "$P2_NAME ($P2_BIN not built)"
fi

# ---- Path 03 ---- codegen + run.hpp
P3_NAME="03: codegen + run.hpp"
ARGS=(--no-regen --build "$BUILD")
[ "$IS_MPI" = 1 ] && ARGS=(--variant "$VARIANT" --np "$NP" "${ARGS[@]}")
[ "$VARIANT" = "gpu" ] && ARGS=(--variant gpu "${ARGS[@]}")
if bash "$APP/validate_codegen.sh" "${ARGS[@]}" poisson2d > /tmp/tut_03.log 2>&1 \
    && grep -q "OK ]" /tmp/tut_03.log; then
    record_pass "$P3_NAME"
else
    record_fail "$P3_NAME"
fi

# ---- Path 04 ---- codegen + facade
P4_NAME="04: codegen + ExasimSolver<M>::load_pdeapp"
if bash "$APP/validate_codegen.sh" --facade "${ARGS[@]}" poisson2d > /tmp/tut_04.log 2>&1 \
    && grep -q "OK ]" /tmp/tut_04.log; then
    record_pass "$P4_NAME"
else
    record_fail "$P4_NAME"
fi

# ---- Path 05 ---- codegen + AbiAdapter
case "$VARIANT" in
    cpu)     P5_NAME_BIN="cput2cEXASIM" ;;
    mpi)     P5_NAME_BIN="cpumpit2cEXASIM" ;;
    gpu)     P5_NAME_BIN="gput2cEXASIM" ;;
    mpi_gpu) P5_NAME_BIN="gpumpit2cEXASIM" ;;
esac
P5_NAME="05: codegen + AbiAdapter ($P5_NAME_BIN)"
P5_BIN="$BUILD/$P5_NAME_BIN"
if [ -x "$P5_BIN" ]; then
    cd "$APP/poisson2d_codegen"
    if "$EXASIM/build/text2code" ./pdeapp.txt > /tmp/tut_05_t2c.log 2>&1 \
        && "${RUNNER[@]}" "$P5_BIN" ./pdeapp.txt > /tmp/tut_05.log 2>&1 \
        && grep -q "Updated Norm:" /tmp/tut_05.log; then
        record_pass "$P5_NAME"
    else
        record_fail "$P5_NAME"
    fi
else
    record_skip "$P5_NAME ($P5_BIN not built; needs WITH_TEXT2CODE=ON)"
fi

echo
echo "[tutorial] $VARIANT summary:"
echo "  passed:  ${#PASSED[@]}"
echo "  failed:  ${#FAILED[@]}"
echo "  skipped: ${#SKIPPED[@]}"
if [ ${#FAILED[@]} -gt 0 ]; then
    echo
    for f in "${FAILED[@]}"; do echo "  - $f"; done
    exit 1
fi
exit 0
