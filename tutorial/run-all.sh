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

# Default build dir per variant if not provided.
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
    cpu)     SUFFIX="";          RUNNER=();                       IS_MPI=0 ;;
    gpu)     SUFFIX="_gpu";      RUNNER=();                       IS_MPI=0 ;;
    mpi)     SUFFIX="_mpi";      RUNNER=(mpirun -np "$NP");       IS_MPI=1 ;;
    mpi_gpu) SUFFIX="_mpi_gpu";  RUNNER=(mpirun -np "$NP");       IS_MPI=1 ;;
esac

# Make BUILD absolute since the subshells `cd` into example dirs.
if [ -d "$BUILD" ]; then
    BUILD="$(cd "$BUILD" && pwd)"
fi

APP="$EXASIM/apps/library_example"

echo "[tutorial] VARIANT = $VARIANT (suffix=$SUFFIX)"
echo "[tutorial] BUILD   = $BUILD"
[ "$IS_MPI" = 1 ] && echo "[tutorial] NP      = $NP"

echo "[tutorial] regenerating poisson2d codegen artifacts ..."
bash "$APP/regenerate.sh" poisson2d > /dev/null

pass() { echo "[ OK ]  $1"; }
fail() { echo "[FAIL] $1"; exit 1; }

# ---- Path 01 ---- hand-written + run.hpp
# poisson2d_template doesn't have multi-backend variants in the
# library_example layout — only the codegen / facade families do.
# Skip path 01 for non-CPU variants and note it.
P1_BIN="$BUILD/poisson2d_template${SUFFIX}"
echo "[tutorial] path 01 — hand-written + legacy CLI"
if [ -x "$P1_BIN" ]; then
    cd "$APP/poisson2d"
    rm -rf datain dataout
    mkdir -p datain dataout
    if [ "$IS_MPI" = 1 ]; then
        "${RUNNER[@]}" "$P1_BIN" ./pdeapp.txt > /tmp/tut_01.log 2>&1 || fail "01: $P1_BIN"
    else
        "$P1_BIN" ./pdeapp.txt > /tmp/tut_01.log 2>&1 || fail "01: $P1_BIN"
    fi
    [ -f "dataout/outudg_np0.bin" ] || fail "01: no outudg"
    pass "01"
else
    echo "[SKIP] 01: $P1_BIN not built (poisson2d_template only built on CPU)"
fi

# ---- Path 02 ---- hand-written + facade
# solve_square_facade reads master-node tables relative to EXASIM_DIR
# at runtime; pass it via env var.
P2_BIN="$BUILD/solve_square_facade${SUFFIX}"
echo "[tutorial] path 02 — hand-written + ExasimSolver<M>"
if [ -x "$P2_BIN" ]; then
    cd "$EXASIM"
    if [ "$IS_MPI" = 1 ]; then
        EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$P2_BIN" > /tmp/tut_02.log 2>&1 || fail "02: $P2_BIN"
    else
        EXASIM_DIR="$EXASIM" "$P2_BIN" > /tmp/tut_02.log 2>&1 || fail "02: $P2_BIN"
    fi
    grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_02.log || fail "02: solve_square didn't converge"
    pass "02"
else
    echo "[SKIP] 02: $P2_BIN not built"
fi

# ---- Path 03 ---- codegen + run.hpp
echo "[tutorial] path 03 — codegen + legacy CLI"
ARGS=(--no-regen --build "$BUILD")
[ "$IS_MPI" = 1 ] && ARGS=(--variant "$VARIANT" --np "$NP" "${ARGS[@]}")
[ "$VARIANT" = "gpu" ]     && ARGS=(--variant gpu     "${ARGS[@]}")
[ "$VARIANT" = "mpi_gpu" ] && ARGS=("${ARGS[@]}")  # already set above
bash "$APP/validate_codegen.sh" "${ARGS[@]}" poisson2d > /tmp/tut_03.log 2>&1 \
    || fail "03: validate_codegen.sh poisson2d"
grep -q "OK ]" /tmp/tut_03.log || fail "03: no OK"
pass "03"

# ---- Path 04 ---- codegen + facade
echo "[tutorial] path 04 — codegen + ExasimSolver<M>::load_pdeapp"
bash "$APP/validate_codegen.sh" --facade "${ARGS[@]}" poisson2d > /tmp/tut_04.log 2>&1 \
    || fail "04: validate_codegen.sh --facade poisson2d"
grep -q "OK ]" /tmp/tut_04.log || fail "04: no OK"
pass "04"

# ---- Path 05 ---- codegen + AbiAdapter
# Only attempted on CPU single-rank since cput2cEXASIM is gated on
# WITH_TEXT2CODE=ON. The binary's RPATH is hardcoded to
# backend/Model/, so we re-run text2code without --out-dir to
# populate that location before invoking the binary.
P5_BIN="$BUILD/cput2cEXASIM"
if [ "$VARIANT" = "cpu" ] && [ -x "$P5_BIN" ]; then
    echo "[tutorial] path 05 — codegen + AbiAdapter"
    cd "$APP/poisson2d_codegen"
    "$EXASIM/build/text2code" ./pdeapp.txt > /tmp/tut_05_t2c.log 2>&1 \
        || fail "05: text2code failed"
    "$P5_BIN" ./pdeapp.txt > /tmp/tut_05.log 2>&1 \
        || fail "05: $P5_BIN"
    grep -q "Updated Norm:" /tmp/tut_05.log || fail "05: no Newton convergence"
    pass "05"
else
    echo "[SKIP] 05: cput2cEXASIM not available for variant $VARIANT (needs WITH_TEXT2CODE=ON, CPU)"
fi

echo "[tutorial] $VARIANT: all paths pass"
