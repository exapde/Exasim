#!/usr/bin/env bash
# Run the tutorial sections against poisson2d for one backend
# variant. Each section is a self-contained subdirectory with its
# own `main.cpp`, `pdeapp.txt`, `CMakeLists.txt`, and so on. This
# script generates the codegen artifacts where needed, builds each
# tutorial binary, runs it, and verifies the output.
#
# Sections are ordered by decreasing automation:
#   01 — generated + prebuilt    (drop two text files; prebuilt binary)
#   02 — generated + CLI         (3-line main, header-only)
#   03 — generated + embedded    (load_pdeapp from C++)
#   04 — handwritten + CLI
#   05 — handwritten + embedded  (set_mesh, in-memory)
#   06 — handwritten + distributed (set_mesh_distributed; MPI only)
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
    cpu)     RUNNER=();                 IS_MPI=0 ;;
    gpu)     RUNNER=();                 IS_MPI=0 ;;
    mpi)     RUNNER=(mpirun -np "$NP"); IS_MPI=1 ;;
    mpi_gpu) RUNNER=(mpirun -np "$NP"); IS_MPI=1 ;;
esac

if [ -d "$BUILD" ]; then
    BUILD="$(cd "$BUILD" && pwd)"
fi

# Tutorial CMake targets carry no variant suffix; each backend
# variant has its own build dir, which already disambiguates the
# binary location.
echo "[tutorial] VARIANT = $VARIANT"
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
# directory. The generated-CLI and generated-embedded sections need
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

set_mpiprocs() {
    # Rewrites pdeapp.txt → pdeapp_run.txt with mpiprocs adjusted to
    # match the variant's actual rank count.
    local in="$1" out="$2"
    if [ "$IS_MPI" = 1 ]; then
        sed "s/^mpiprocs *=.*/mpiprocs = $NP;/" "$in" > "$out"
    else
        sed 's/^mpiprocs *=.*/mpiprocs = 1;/' "$in" > "$out"
    fi
}

# Generate grid.bin in a tutorial section that needs one.
# The script imports SquareMesh from frontends/Python/Mesh and
# writes a 16x16 Cartesian quad mesh on the unit square in the
# legacy writebin format. Pure NumPy; no gmsh dependency.
PY="${PYTHON:-python3}"
generate_grid_bin() {
    local dir="$1"
    "$PY" "$EXASIM/tutorial/tools/squaregrid.py" 16 \
        "$TUT/$dir/grid.bin" \
        > "/tmp/tut_${dir}_grid.log" 2>&1
}

# ---- Section 01 ---- generated + prebuilt ----------------------
case "$VARIANT" in
    cpu)     S1_BIN_NAME="cput2cEXASIM" ;;
    mpi)     S1_BIN_NAME="cpumpit2cEXASIM" ;;
    gpu)     S1_BIN_NAME="gput2cEXASIM" ;;
    mpi_gpu) S1_BIN_NAME="gpumpit2cEXASIM" ;;
esac
S1_NAME="01: generated + prebuilt ($S1_BIN_NAME)"
S1_BIN="$BUILD/$S1_BIN_NAME"
if [ -x "$S1_BIN" ]; then
    cd "$TUT/01-generated-prebuilt"
    rm -rf datain dataout
    mkdir -p datain dataout
    generate_grid_bin "01-generated-prebuilt"
    set_mpiprocs pdeapp.txt pdeapp_run.txt
    if "$EXASIM/build/text2code" ./pdeapp_run.txt > /tmp/tut_01_t2c.log 2>&1 \
        && "${RUNNER[@]}" "$S1_BIN" ./pdeapp_run.txt > /tmp/tut_01.log 2>&1 \
        && grep -q "Updated Norm:" /tmp/tut_01.log; then
        record_pass "$S1_NAME"
    else
        record_fail "$S1_NAME"
    fi
else
    record_skip "$S1_NAME ($S1_BIN not built; configure with WITH_TEXT2CODE=ON)"
fi

# ---- Section 02 ---- generated + CLI ---------------------------
S2_NAME="02: generated + CLI"
S2_TARGET="tutorial_02_generated_cli"
S2_BIN="$BUILD/$S2_TARGET"
if generate_codegen "02-generated-cli" \
   && cmake "$BUILD" -B "$BUILD" > /dev/null 2>&1 \
   && build_target "$S2_TARGET" \
   && [ -x "$S2_BIN" ]; then
    cd "$TUT/02-generated-cli"
    rm -rf datain dataout
    mkdir -p datain dataout
    generate_grid_bin "02-generated-cli"
    set_mpiprocs pdeapp.txt pdeapp_run.txt
    if "${RUNNER[@]}" "$S2_BIN" ./pdeapp_run.txt > /tmp/tut_02.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$S2_NAME"
    else
        record_fail "$S2_NAME ($S2_BIN)"
    fi
else
    record_skip "$S2_NAME (target $S2_TARGET not buildable)"
fi

# ---- Section 03 ---- generated + embedded ----------------------
S3_NAME="03: generated + embedded"
S3_TARGET="tutorial_03_generated_embedded"
S3_BIN="$BUILD/$S3_TARGET"
if generate_codegen "03-generated-embedded" \
   && cmake "$BUILD" -B "$BUILD" > /dev/null 2>&1 \
   && build_target "$S3_TARGET" \
   && [ -x "$S3_BIN" ]; then
    cd "$TUT/03-generated-embedded"
    rm -rf datain dataout
    mkdir -p datain dataout
    generate_grid_bin "03-generated-embedded"
    set_mpiprocs pdeapp.txt pdeapp_run.txt
    if "${RUNNER[@]}" "$S3_BIN" ./pdeapp_run.txt > /tmp/tut_03.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$S3_NAME"
    else
        record_fail "$S3_NAME ($S3_BIN)"
    fi
else
    record_skip "$S3_NAME (target $S3_TARGET not buildable)"
fi

# ---- Section 04 ---- handwritten + CLI -------------------------
S4_NAME="04: handwritten + CLI"
S4_TARGET="tutorial_04_handwritten_cli"
S4_BIN="$BUILD/$S4_TARGET"
if build_target "$S4_TARGET" && [ -x "$S4_BIN" ]; then
    cd "$TUT/04-handwritten-cli"
    rm -rf datain dataout
    mkdir -p datain dataout
    generate_grid_bin "04-handwritten-cli"
    set_mpiprocs pdeapp.txt pdeapp_run.txt
    if "${RUNNER[@]}" "$S4_BIN" ./pdeapp_run.txt > /tmp/tut_04.log 2>&1 \
        && [ -f "dataout/outudg_np0.bin" ]; then
        record_pass "$S4_NAME"
    else
        record_fail "$S4_NAME ($S4_BIN)"
    fi
else
    record_skip "$S4_NAME (target $S4_TARGET not buildable)"
fi

# ---- Section 05 ---- handwritten + embedded --------------------
S5_NAME="05: handwritten + embedded"
S5_TARGET="tutorial_05_handwritten_embedded"
S5_BIN="$BUILD/$S5_TARGET"
if build_target "$S5_TARGET" && [ -x "$S5_BIN" ]; then
    cd "$TUT/05-handwritten-embedded"
    if EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$S5_BIN" > /tmp/tut_05.log 2>&1 \
        && grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_05.log; then
        record_pass "$S5_NAME"
    else
        record_fail "$S5_NAME ($S5_BIN)"
    fi
else
    record_skip "$S5_NAME (target $S5_TARGET not buildable)"
fi

# ---- Section 06 ---- handwritten + distributed (MPI) -----------
S6_NAME="06: handwritten + distributed"
S6_TARGET="tutorial_06_handwritten_distributed"
S6_BIN="$BUILD/$S6_TARGET"
if [ "$IS_MPI" != 1 ]; then
    record_skip "$S6_NAME (variant $VARIANT is not MPI; section 06 is MPI-only)"
elif build_target "$S6_TARGET" && [ -x "$S6_BIN" ]; then
    cd "$TUT/06-handwritten-distributed"
    if EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$S6_BIN" > /tmp/tut_06.log 2>&1 \
        && grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_06.log; then
        record_pass "$S6_NAME"
    else
        record_fail "$S6_NAME ($S6_BIN)"
    fi
else
    record_skip "$S6_NAME (target $S6_TARGET not buildable)"
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
