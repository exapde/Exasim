#!/usr/bin/env bash
# Run the tutorial paths against poisson2d for one backend variant.
# Each path is a self-contained subdirectory with its own
# `main.cpp`, `pdeapp.txt`, `CMakeLists.txt`, and so on. This
# script generates the codegen artifacts where needed, builds each
# tutorial binary, runs it, and verifies the output.
#
# Usage:  bash tutorial/run-all.sh [--variant cpu|gpu|mpi|mpi_gpu]
#                                  [--build <dir>] [--np N]
#
# Defaults:  --variant cpu   --build $EXASIM/build   --np 2 (mpi variants)

set -uo pipefail

EXASIM="$(cd "$(dirname "$0")/.." && pwd)"
TUT="$EXASIM/tutorial"
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

echo "[tutorial] VARIANT = $VARIANT (suffix=$SUFFIX)"
echo "[tutorial] BUILD   = $BUILD"
[ "$IS_MPI" = 1 ] && echo "[tutorial] NP      = $NP"

PASSED=()
FAILED=()
SKIPPED=()
record_pass() { PASSED+=("$1"); echo "[ OK ]   $1"; }
record_fail() { FAILED+=("$1"); echo "[FAIL]   $1"; }
record_skip() { SKIPPED+=("$1"); echo "[SKIP]   $1"; }

# Run text2code on a tutorial section's pdeapp.txt to populate
# my_model.hpp and libpdemodel*.{so,dylib} in that section's
# directory. The codegen-CLI and codegen-facade tutorials need
# this to exist before `cmake --build` can link the binary.
generate_codegen() {
    local dir="$1"
    ( cd "$TUT/$dir" && "$EXASIM/build/text2code" ./pdeapp.txt --out-dir "$TUT/$dir" ) \
        > "/tmp/tut_${dir}_t2c.log" 2>&1
}

# Build a tutorial target via cmake. Returns 0 on success.
build_target() {
    local target="$1"
    cmake --build "$BUILD" --target "$target" \
        > "/tmp/tut_${target}_build.log" 2>&1
}

# ---- Path 01 ---- hand-written + run.hpp -----------------------
P1_NAME="01: hand-written + run.hpp"
P1_TARGET="tutorial_01_handwritten_cli${SUFFIX}"
P1_BIN="$BUILD/$P1_TARGET"
if build_target "$P1_TARGET" && [ -x "$P1_BIN" ]; then
    cd "$TUT/01-handwritten-cli"
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
    record_skip "$P1_NAME (target $P1_TARGET not buildable)"
fi

# ---- Path 02 ---- hand-written + facade ------------------------
P2_NAME="02: hand-written + ExasimSolver<M>"
P2_TARGET="tutorial_02_handwritten_facade${SUFFIX}"
P2_BIN="$BUILD/$P2_TARGET"
if build_target "$P2_TARGET" && [ -x "$P2_BIN" ]; then
    cd "$TUT/02-handwritten-facade"
    if EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$P2_BIN" > /tmp/tut_02.log 2>&1 \
        && grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_02.log; then
        record_pass "$P2_NAME"
    else
        record_fail "$P2_NAME ($P2_BIN)"
    fi
else
    record_skip "$P2_NAME (target $P2_TARGET not buildable)"
fi

# ---- Path 03 ---- codegen + run.hpp ----------------------------
P3_NAME="03: codegen + run.hpp"
P3_TARGET="tutorial_03_codegen_cli${SUFFIX}"
P3_BIN="$BUILD/$P3_TARGET"
if generate_codegen "03-codegen-cli" \
   && cmake "$BUILD" -B "$BUILD" > /dev/null 2>&1 \
   && build_target "$P3_TARGET" \
   && [ -x "$P3_BIN" ]; then
    cd "$TUT/03-codegen-cli"
    rm -rf datain dataout
    mkdir -p datain dataout
    if [ "$IS_MPI" = 1 ]; then
        sed "s/^mpiprocs *=.*/mpiprocs = $NP;/" pdeapp.txt > pdeapp_run.txt
    else
        sed 's/^mpiprocs *=.*/mpiprocs = 1;/' pdeapp.txt > pdeapp_run.txt
    fi
    if "${RUNNER[@]}" "$P3_BIN" ./pdeapp_run.txt > /tmp/tut_03.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$P3_NAME"
    else
        record_fail "$P3_NAME ($P3_BIN)"
    fi
else
    record_skip "$P3_NAME (target $P3_TARGET not buildable)"
fi

# ---- Path 04 ---- codegen + facade -----------------------------
P4_NAME="04: codegen + ExasimSolver<M>::load_pdeapp"
P4_TARGET="tutorial_04_codegen_facade${SUFFIX}"
P4_BIN="$BUILD/$P4_TARGET"
if generate_codegen "04-codegen-facade" \
   && cmake "$BUILD" -B "$BUILD" > /dev/null 2>&1 \
   && build_target "$P4_TARGET" \
   && [ -x "$P4_BIN" ]; then
    cd "$TUT/04-codegen-facade"
    rm -rf datain dataout
    mkdir -p datain dataout
    if [ "$IS_MPI" = 1 ]; then
        sed "s/^mpiprocs *=.*/mpiprocs = $NP;/" pdeapp.txt > pdeapp_run.txt
    else
        sed 's/^mpiprocs *=.*/mpiprocs = 1;/' pdeapp.txt > pdeapp_run.txt
    fi
    if "${RUNNER[@]}" "$P4_BIN" ./pdeapp_run.txt > /tmp/tut_04.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$P4_NAME"
    else
        record_fail "$P4_NAME ($P4_BIN)"
    fi
else
    record_skip "$P4_NAME (target $P4_TARGET not buildable)"
fi

# ---- Path 05 ---- codegen + AbiAdapter -------------------------
case "$VARIANT" in
    cpu)     P5_NAME_BIN="cput2cEXASIM" ;;
    mpi)     P5_NAME_BIN="cpumpit2cEXASIM" ;;
    gpu)     P5_NAME_BIN="gput2cEXASIM" ;;
    mpi_gpu) P5_NAME_BIN="gpumpit2cEXASIM" ;;
esac
P5_NAME="05: codegen + AbiAdapter ($P5_NAME_BIN)"
P5_BIN="$BUILD/$P5_NAME_BIN"
if [ -x "$P5_BIN" ]; then
    cd "$TUT/05-codegen-abi-adapter"
    rm -rf datain dataout
    mkdir -p datain dataout
    if [ "$IS_MPI" = 1 ]; then
        sed "s/^mpiprocs *=.*/mpiprocs = $NP;/" pdeapp.txt > pdeapp_run.txt
    else
        sed 's/^mpiprocs *=.*/mpiprocs = 1;/' pdeapp.txt > pdeapp_run.txt
    fi
    if "$EXASIM/build/text2code" ./pdeapp_run.txt > /tmp/tut_05_t2c.log 2>&1 \
        && "${RUNNER[@]}" "$P5_BIN" ./pdeapp_run.txt > /tmp/tut_05.log 2>&1 \
        && grep -q "Updated Norm:" /tmp/tut_05.log; then
        record_pass "$P5_NAME"
    else
        record_fail "$P5_NAME"
    fi
else
    record_skip "$P5_NAME ($P5_BIN not built; configure with WITH_TEXT2CODE=ON)"
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
