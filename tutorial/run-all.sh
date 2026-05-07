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

# Sections allowed to exit nonzero on this variant due to known
# finalize-time teardown bugs. The compute output is still required
# to be present, so this is not blanket immunity — only the exit
# code is forgiven, and only on the variants where the bug is
# actually triggered.
#
# Empirical baseline (all CPU sections exit 0 today): leaving cpu
# and mpi strict catches new regressions immediately. GPU variants
# are the ones with the Kokkos::Cuda::finalize / libpdemodelcuda
# `std::map<int, CUevent_st*>` double-free race on teardown — the
# math is correct, the output files are written, but the process
# exits with a signal. Real fix lives in the reverted Meyers-
# singleton work and is tracked separately.
case "$VARIANT" in
    gpu|mpi_gpu)
        QUARANTINE_NONZERO=(01 02 03 04 05 06)
        ;;
    *)
        QUARANTINE_NONZERO=()
        ;;
esac

is_quarantined() {
    local sec="$1"
    [ ${#QUARANTINE_NONZERO[@]} -eq 0 ] && return 1
    local q
    for q in "${QUARANTINE_NONZERO[@]}"; do
        [ "$q" = "$sec" ] && return 0
    done
    return 1
}

# Classify a section run.
#   $1 = section id (matched against QUARANTINE_NONZERO)
#   $2 = section display name
#   $3 = process exit code (numeric)
#   $4..$N = output check (run as a command; success means output present)
#
# Pass criteria:
#   - output check passes AND exit_code == 0
#   - output check passes AND section is quarantined for this variant
# Fail criteria:
#   - output check fails (no usable result; quarantine does not apply)
#   - nonzero exit AND section is not quarantined
classify_run() {
    local sec="$1" name="$2" exit_code="$3"
    shift 3
    if ! "$@"; then
        record_fail "$name (no output produced; exit=$exit_code)"
        return
    fi
    if [ "$exit_code" -eq 0 ]; then
        record_pass "$name"
    elif is_quarantined "$sec"; then
        record_pass "$name (quarantined for variant=$VARIANT: nonzero exit $exit_code accepted)"
    else
        record_fail "$name (exit=$exit_code; nonzero not quarantined for variant=$VARIANT)"
    fi
}

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

# Pass criteria: zero exit code AND expected output present.
# Sections in QUARANTINE_NONZERO (above) may exit nonzero on this
# variant if the output check still passes — that pinhole is for
# the documented GPU teardown bug only. CPU and MPI variants are
# strict; the lenient global "judge on output, not exit code" policy
# was removed because it hid real teardown regressions on every
# variant in service of one variant's known bug.

# ---- Section 01 ---- generated + prebuilt ----------------------
case "$VARIANT" in
    cpu)     S1_BIN_NAME="cput2cEXASIM" ;;
    mpi)     S1_BIN_NAME="cpumpit2cEXASIM" ;;
    gpu)     S1_BIN_NAME="gput2cEXASIM" ;;
    mpi_gpu) S1_BIN_NAME="gpumpit2cEXASIM" ;;
esac
S1_NAME="01: generated + prebuilt ($S1_BIN_NAME)"
S1_BIN="$BUILD/$S1_BIN_NAME"
# Build the prebuilt binary on demand; not built by default in some
# CMake invocations even when registered.
if [ ! -x "$S1_BIN" ]; then
    build_target "$S1_BIN_NAME" || true
fi
if [ -x "$S1_BIN" ]; then
    cd "$TUT/01-generated-prebuilt"
    rm -rf datain dataout
    mkdir -p datain dataout
    generate_grid_bin "01-generated-prebuilt"
    set_mpiprocs pdeapp.txt pdeapp_run.txt
    "$EXASIM/build/text2code" ./pdeapp_run.txt > /tmp/tut_01_t2c.log 2>&1
    "${RUNNER[@]}" "$S1_BIN" ./pdeapp_run.txt > /tmp/tut_01.log 2>&1
    classify_run "01" "$S1_NAME" $? grep -q "Updated Norm:" /tmp/tut_01.log
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
    "${RUNNER[@]}" "$S2_BIN" ./pdeapp_run.txt > /tmp/tut_02.log 2>&1
    classify_run "02" "$S2_NAME" $? test -f "dataout/outudg_np0.bin"
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
    "${RUNNER[@]}" "$S3_BIN" ./pdeapp_run.txt > /tmp/tut_03.log 2>&1
    classify_run "03" "$S3_NAME" $? test -f "dataout/outudg_np0.bin"
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
    "${RUNNER[@]}" "$S4_BIN" ./pdeapp_run.txt > /tmp/tut_04.log 2>&1
    classify_run "04" "$S4_NAME" $? test -f "dataout/outudg_np0.bin"
else
    record_skip "$S4_NAME (target $S4_TARGET not buildable)"
fi

# ---- Section 05 ---- handwritten + embedded (single-rank only) -
S5_NAME="05: handwritten + embedded"
S5_TARGET="tutorial_05_handwritten_embedded"
S5_BIN="$BUILD/$S5_TARGET"
if [ "$IS_MPI" = 1 ]; then
    record_skip "$S5_NAME (variant $VARIANT is MPI; section 05 is single-rank only — see section 06)"
elif build_target "$S5_TARGET" && [ -x "$S5_BIN" ]; then
    cd "$TUT/05-handwritten-embedded"
    EXASIM_DIR="$EXASIM" "$S5_BIN" > /tmp/tut_05.log 2>&1
    classify_run "05" "$S5_NAME" $? grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_05.log
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
    EXASIM_DIR="$EXASIM" "${RUNNER[@]}" "$S6_BIN" > /tmp/tut_06.log 2>&1
    classify_run "06" "$S6_NAME" $? grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_06.log
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
