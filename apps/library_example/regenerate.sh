#!/usr/bin/env bash
# regenerate.sh — populate a codegen example's `grid.bin`, `pdemodel.txt`,
# `datain/`, and `my_model.hpp` from the corresponding source app and
# the retargeted text2code emission.
#
# Usage:  bash regenerate.sh <app-name>
#         e.g.  bash regenerate.sh poisson3d
#
# Source layout assumed:
#   apps/poisson/<name>/             OR  apps/navierstokes/<name>/
#       grid.bin
#       pdemodel.txt
#   apps/library_example/<name>_codegen/
#       pdeapp.txt   (committed; serial-friendly: mpiprocs=1, gendatain=1,
#                                                exasimpath=...)
#       main.cpp     (committed; instantiates CSolution<GeneratedModel>)
#
# This script:
#   1. Copies grid.bin + pdemodel.txt from the source app.
#   2. Runs text2code on the codegen-app pdeapp.txt — which (a) preprocesses
#      the mesh into datain/, (b) writes my_model.hpp + libpdemodelserial.dylib
#      into backend/Model/. Step (b) leaves the legacy ABI's
#      libpdemodelserial.dylib pointing at THIS app's kernels — that's
#      fine because the codegen target only consumes my_model.hpp; the
#      legacy regression in baseline/verify.sh re-runs text2code on
#      apps/poisson/poisson2d before its check.
#   3. Copies the freshly-generated my_model.hpp into the codegen-app dir.

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <app-name>" >&2
    exit 1
fi
NAME="$1"

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CG_DIR="$ROOT/apps/library_example/${NAME}_codegen"

if [ ! -d "$CG_DIR" ]; then
    echo "error: $CG_DIR doesn't exist" >&2
    exit 1
fi

# Locate source app (poisson/* or navierstokes/*).
SRC_DIR=""
for cand in "$ROOT/apps/poisson/$NAME" "$ROOT/apps/navierstokes/$NAME"; do
    if [ -d "$cand" ]; then SRC_DIR="$cand"; break; fi
done
if [ -z "$SRC_DIR" ]; then
    echo "error: no source app found at apps/poisson/$NAME or apps/navierstokes/$NAME" >&2
    exit 1
fi

cp "$SRC_DIR/grid.bin"     "$CG_DIR/grid.bin"
cp "$SRC_DIR/pdemodel.txt" "$CG_DIR/pdemodel.txt"

# Copy any optional companion bins the app's pdeapp.txt expects.
# (isoq3d has xdg.bin; NS apps often have udg/vdg/wdg/xdg.bin for
# initial conditions; some have partition.bin from earlier runs.)
for opt in xdg.bin udg.bin vdg.bin wdg.bin partition.bin; do
    if [ -f "$SRC_DIR/$opt" ]; then
        cp "$SRC_DIR/$opt" "$CG_DIR/$opt"
    fi
done

cd "$CG_DIR"
# HOT.7.14 — write generated my_model.hpp + libpdemodel*.{so,dylib}
# into this example's own dir. The shared backend/Model/ output was
# the root cause of the test-matrix RESOURCE_LOCK serialization;
# per-target dirs let parallel ctest gates run without clobbering.
# Find a built text2code binary: prefer build/, fall back to any
# build_*/ (e.g. build_cpu/, build_mpi_local/). The codegen-only
# build dirs don't include text2code as a target.
TEXT2CODE=""
for cand in "$ROOT/build/text2code" "$ROOT"/build_*/text2code; do
    if [ -x "$cand" ]; then TEXT2CODE="$cand"; break; fi
done
if [ -z "$TEXT2CODE" ]; then
    echo "error: text2code binary not found (tried $ROOT/build/text2code and $ROOT/build_*/text2code)" >&2
    exit 1
fi
"$TEXT2CODE" ./pdeapp.txt --out-dir "$CG_DIR" >/dev/null

# Force the codegen / facade targets that link backend/Model/libpdemodel.cpp
# to rebuild on the next `cmake --build`: that TU `#include`s
# per-example Kokkos*.cpp / Hdg*.cpp / cpuInit*.cpp via the -I path
# order (per-example first, backend/Model fallback). The per-example
# files we just wrote aren't in the compiler-emitted .d (those record
# the fallback that existed at initial compile time), so without this
# touch make sees the .o as up-to-date and skips recompilation. The
# ctest RESOURCE_LOCK on backend_model_files serializes regens vs
# builds across all 12 examples, so `touch`-ing the shared TU is safe
# — each subsequent `cmake --build --target X_*` recompiles X's
# own libpdemodel.cpp.o against X's per-example dir.
touch "$ROOT/backend/Model/libpdemodel.cpp"

echo "[ OK ] regenerated ${NAME}_codegen (grid.bin, pdemodel.txt, datain/, my_model.hpp, libpdemodel*)"
