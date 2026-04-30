#!/usr/bin/env bash
# validate_codegen.sh — run each codegen example in-place from
# apps/library_example/<name>_codegen/ and gate against the
# baseline/<name>_serial/ recording.
#
# Two gate modes:
#   - QoI:   if baseline has outqoi.txt, diff `dataout/outqoi.txt`.
#            (Preferred — robust to floating-point reordering.)
#   - md5:   else compare md5s of `dataout/*.bin` against baseline md5s.
#
# Assumes:
#   - text2code has been built into ./build/
#   - <name>_codegen targets are built into ./build/
#   - bash regenerate.sh <name> has populated grid.bin / pdemodel.txt /
#     datain/ / my_model.hpp for each example.
#
# Exit 0 if all configured examples pass; non-zero otherwise.

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LIB="$ROOT/apps/library_example"
BASE="$ROOT/baseline"
BUILD="$ROOT/build"

EXAMPLES=("$@")
if [ ${#EXAMPLES[@]} -eq 0 ]; then
    EXAMPLES=(poisson2d poisson3d periodic)
fi

fail=0

for name in "${EXAMPLES[@]}"; do
    cg_dir="$LIB/${name}_codegen"
    base_dir="$BASE/${name}_serial"
    bin="$BUILD/${name}_codegen"

    if [ ! -x "$bin" ]; then
        echo "[SKIP] $name — $bin not built (cmake --build build --target ${name}_codegen)"
        continue
    fi
    if [ ! -d "$cg_dir/datain" ] || [ ! -f "$cg_dir/grid.bin" ]; then
        echo "[SKIP] $name — $cg_dir not populated (bash regenerate.sh $name)"
        continue
    fi
    if [ ! -d "$base_dir" ]; then
        echo "[SKIP] $name — no baseline at $base_dir"
        continue
    fi

    rm -rf "$cg_dir/dataout"
    mkdir -p "$cg_dir/dataout"
    ( cd "$cg_dir" && "$bin" ./pdeapp.txt >/dev/null )

    if [ -f "$base_dir/outqoi.txt" ]; then
        if diff -q "$cg_dir/dataout/outqoi.txt" "$base_dir/outqoi.txt" >/dev/null; then
            echo "[ OK ] $name (QoI byte-identical)"
        else
            echo "[FAIL] $name — outqoi.txt differs"
            diff "$cg_dir/dataout/outqoi.txt" "$base_dir/outqoi.txt" || true
            fail=1
        fi
    elif [ -f "$base_dir/md5.txt" ]; then
        got="$( cd "$cg_dir/dataout" && md5 -r *.bin 2>/dev/null | sort )"
        want="$(sort "$base_dir/md5.txt")"
        if [ "$got" = "$want" ]; then
            echo "[ OK ] $name (bin md5s match)"
        else
            echo "[FAIL] $name — bin md5s differ"
            diff <(echo "$want") <(echo "$got") || true
            fail=1
        fi
    else
        echo "[SKIP] $name — baseline has no outqoi.txt or md5.txt"
    fi
done

exit "$fail"
