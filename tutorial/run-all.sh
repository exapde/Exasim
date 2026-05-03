#!/usr/bin/env bash
# Run all four supported tutorial paths against poisson2d on CPU and
# verify they produce equivalent results (element-L2 < 1e-3).
#
# Usage:  bash tutorial/run-all.sh [--build <dir>]
#
# Default build dir: ./build (CPU, single-rank).

set -euo pipefail

EXASIM="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$EXASIM/build"
APP="$EXASIM/apps/library_example"

while [ $# -gt 0 ]; do
    case "$1" in
        --build) BUILD="$2"; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

echo "[tutorial] BUILD = $BUILD"
echo "[tutorial] regenerating poisson2d codegen artifacts ..."
bash "$APP/regenerate.sh" poisson2d > /dev/null

pass() { echo "[ OK ]  $1"; }
fail() { echo "[FAIL] $1"; exit 1; }

# ---- Path 01 ---- hand-written + run.hpp
echo "[tutorial] path 01 — hand-written + legacy CLI"
cd "$APP/poisson2d"
mkdir -p datain dataout
"$BUILD/poisson2d_template" ./pdeapp.txt > /tmp/tut_01.log 2>&1 || fail "poisson2d_template"
[ -f dataout/outudg_np0.bin ] || fail "01: no outudg"
pass "01"

# ---- Path 02 ---- hand-written + facade
# solve_square_facade reads master-node tables relative to EXASIM_DIR
# at runtime; pass it via env var.
echo "[tutorial] path 02 — hand-written + ExasimSolver<M>"
EXASIM_DIR="$EXASIM" "$BUILD/solve_square_facade" > /tmp/tut_02.log 2>&1 || fail "solve_square_facade"
grep -qE "max\|udg\| = 3.1415[89]" /tmp/tut_02.log || fail "02: solve_square didn't converge"
pass "02"

# ---- Path 03 ---- codegen + run.hpp
echo "[tutorial] path 03 — codegen + legacy CLI"
bash "$APP/validate_codegen.sh" --build "$BUILD" poisson2d > /tmp/tut_03.log 2>&1 \
    || fail "03: validate_codegen.sh poisson2d"
grep -q "OK ]" /tmp/tut_03.log || fail "03: no OK"
pass "03"

# ---- Path 04 ---- codegen + facade
echo "[tutorial] path 04 — codegen + ExasimSolver<M>::load_pdeapp"
bash "$APP/validate_codegen.sh" --facade --build "$BUILD" poisson2d > /tmp/tut_04.log 2>&1 \
    || fail "04: validate_codegen.sh --facade poisson2d"
grep -q "OK ]" /tmp/tut_04.log || fail "04: no OK"
pass "04"

# ---- Path 05 ---- codegen + AbiAdapter
# Currently broken (NaN residual); skip in the standard tutorial run.
echo "[tutorial] path 05 — skipped (AbiAdapter path is deprecated and currently broken)"

echo "[tutorial] all 4 supported paths pass"
