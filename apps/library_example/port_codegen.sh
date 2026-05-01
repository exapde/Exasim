#!/usr/bin/env bash
# port_codegen.sh — automate setting up a codegen example app from
# apps/poisson/<name>/ or apps/navierstokes/<name>/.
#
# Usage:  bash port_codegen.sh <name> [--build] [--validate]
#         e.g.  bash port_codegen.sh poisson3d --build --validate
#
# Steps (idempotent — safe to re-run):
#   1. Locate source app (apps/poisson/<name> or apps/navierstokes/<name>).
#   2. Create apps/library_example/<name>_codegen/ if missing.
#   3. Copy main.cpp from poisson2d_codegen (it's generic, only uses
#      `GeneratedModel`, no app-specific includes).
#   4. Copy pdeapp.txt from source, force serial (mpiprocs=1, gendatain=1).
#   5. Call `regenerate.sh <name>` to:
#        - copy grid.bin + pdemodel.txt from source
#        - run text2code to produce my_model.hpp + datain/
#   6. (--build)    Configure cmake (idempotent) and build <name>_codegen.
#                   Caller must register the target in install/CMakeLists.txt
#                   first (`add_codegen_example(<name>)`); script prints a
#                   reminder if the target is missing.
#   7. Record a baseline by running cput2cEXASIM on the same setup,
#      saving outqoi.txt (if present) or md5 of dataout/*.bin to
#      baseline/<name>_serial/.
#   8. (--validate) Run validate_codegen.sh <name> and exit with its status.
#
# Always idempotent: re-running is the natural debug loop. Edit text2code
# or kernels, rerun the script, see if validation passes.

set -uo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <name> [--build] [--validate]" >&2
    exit 1
fi
NAME="$1"
shift || true

DO_BUILD=0
DO_VALIDATE=0
for arg in "$@"; do
    case "$arg" in
        --build)    DO_BUILD=1    ;;
        --validate) DO_VALIDATE=1 ;;
        *) echo "unknown arg: $arg" >&2; exit 1 ;;
    esac
done

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LIB="$ROOT/apps/library_example"
BASE="$ROOT/baseline"
BUILD="$ROOT/build"

# 1. Locate source app.
SRC=""
for cand in "$ROOT/apps/poisson/$NAME" "$ROOT/apps/navierstokes/$NAME"; do
    if [ -d "$cand" ]; then SRC="$cand"; break; fi
done
if [ -z "$SRC" ]; then
    echo "error: no source app at apps/poisson/$NAME or apps/navierstokes/$NAME" >&2
    exit 1
fi
echo "[port] source: $SRC"

# 2. Set up codegen dir.
CG="$LIB/${NAME}_codegen"
mkdir -p "$CG"

# 3. main.cpp — copied from poisson2d_codegen (generic).
if [ ! -f "$CG/main.cpp" ]; then
    cp "$LIB/poisson2d_codegen/main.cpp" "$CG/main.cpp"
    echo "[port] copied main.cpp from poisson2d_codegen"
fi

# 4. pdeapp.txt — copy from source with mpiprocs=1, gendatain=1, and a
# truncated dt array. naca0012unsteady's `dt = [0.005, 0.01, 0.02,
# 0.04, 0.06, repeat(0.08, 155)]` runs for 160 timesteps (~16 min on
# this hardware) which makes the validation gate impractical. Replace
# any `repeat(<dt>, <n>)` clause with the same dt repeated once, so
# we end up with at most a handful of timesteps per app — enough to
# exercise tdfunc/sdg/DIRK without spending hours per validate.
cp "$SRC/pdeapp.txt" "$CG/pdeapp.txt"
sed -i.bak 's/^mpiprocs = .*/mpiprocs = 1;/; s/^gendatain = .*/gendatain = 1;/' \
    "$CG/pdeapp.txt"
# Replace `repeat(X, N)` -> `X` for any N. Works for the dt schedules
# in apps/navierstokes/naca0012unsteady (the only consumer in apps/).
sed -i.bak2 's/repeat(\([^,]*\), *[0-9]*)/\1/g' "$CG/pdeapp.txt"
rm -f "$CG/pdeapp.txt.bak" "$CG/pdeapp.txt.bak2"
echo "[port] wrote pdeapp.txt (serial + gendatain=1 + repeat() collapsed)"

# 5. Regenerate grid.bin + pdemodel.txt + datain/ + my_model.hpp.
if ! bash "$LIB/regenerate.sh" "$NAME" 2>&1 | tail -1; then
    echo "error: regenerate.sh failed" >&2
    exit 1
fi

# 6. Build (optional).
if [ $DO_BUILD -eq 1 ]; then
    if ! grep -q "add_codegen_example($NAME)" "$ROOT/install/CMakeLists.txt"; then
        echo "[port] reminder: add 'add_codegen_example($NAME)' to install/CMakeLists.txt" >&2
    fi
    ( cd "$BUILD" && cmake . >/dev/null && \
      cmake --build . --target ${NAME}_codegen -j4 2>&1 | tail -3 ) || {
        echo "error: cmake build of ${NAME}_codegen failed" >&2
        exit 1
    }
fi

# 7. Record baseline by running legacy cput2cEXASIM.
if [ -x "$BUILD/cput2cEXASIM" ]; then
    LEGACY="$(mktemp -d)"
    cp -r "$CG/datain" "$LEGACY/datain"
    cp "$CG/grid.bin" "$CG/pdemodel.txt" "$LEGACY/"
    cp "$CG/pdeapp.txt" "$LEGACY/pdeapp.txt"
    # Legacy needs an exasimpath in pdeapp.txt; we DON'T add one to the
    # committed codegen pdeapp.txt because the codegen run uses cwd-derived
    # exasimpath. Add it just for this temp legacy run.
    printf '\nexasimpath = "%s";\n' "$ROOT" >> "$LEGACY/pdeapp.txt"
    mkdir -p "$LEGACY/dataout"

    ( cd "$LEGACY" && "$BUILD/cput2cEXASIM" ./pdeapp.txt >/tmp/${NAME}_legacy.log 2>&1 ) || {
        echo "[port] warn: legacy cput2cEXASIM failed; check /tmp/${NAME}_legacy.log" >&2
        rm -rf "$LEGACY"
    }

    if [ -d "$LEGACY/dataout" ]; then
        BDIR="$BASE/${NAME}_serial"
        mkdir -p "$BDIR"
        if [ -f "$LEGACY/dataout/outqoi.txt" ]; then
            cp "$LEGACY/dataout/outqoi.txt" "$BDIR/outqoi.txt"
            echo "[port] recorded baseline: $BDIR/outqoi.txt"
        fi
        if compgen -G "$LEGACY/dataout/*.bin" > /dev/null; then
            cp "$LEGACY"/dataout/*.bin "$BDIR/" 2>/dev/null
            ( cd "$BDIR" && md5 -r *.bin 2>/dev/null | sort > md5.txt )
            echo "[port] recorded baseline: $BDIR/md5.txt + bins"
        fi
        rm -rf "$LEGACY"
    fi
fi

# 8. Validate (optional).
if [ $DO_VALIDATE -eq 1 ]; then
    bash "$LIB/validate_codegen.sh" "$NAME"
    exit $?
fi

echo "[port] done. Run \`bash apps/library_example/validate_codegen.sh $NAME\` to gate."
