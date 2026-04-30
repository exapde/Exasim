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
"$ROOT/build/text2code" ./pdeapp.txt >/dev/null

cp "$ROOT/backend/Model/my_model.hpp" "$CG_DIR/my_model.hpp"

echo "[ OK ] regenerated ${NAME}_codegen (grid.bin, pdemodel.txt, datain/, my_model.hpp)"
