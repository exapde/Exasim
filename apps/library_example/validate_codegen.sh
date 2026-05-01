#!/usr/bin/env bash
# validate_codegen.sh — run each codegen example in-place from
# apps/library_example/<name>_codegen/ and gate against the
# baseline/<name>_serial/ recording.
#
# Usage:
#   validate_codegen.sh                              # CPU variant
#   validate_codegen.sh --variant gpu                # _codegen_gpu
#   validate_codegen.sh --variant mpi --np 2         # _codegen_mpi
#   validate_codegen.sh --variant mpi_gpu --np 2     # _codegen_mpi_gpu
#   validate_codegen.sh --build build_mpi/           # alt build dir
#   validate_codegen.sh --variant gpu poisson2d      # one example only
#
# Two gate modes:
#   - QoI:   if baseline has outqoi.txt, diff `dataout/outqoi.txt`.
#            (Preferred — robust to floating-point reordering.)
#   - md5:   else compare md5s of `dataout/*.bin` against baseline md5s.
#
# Assumes:
#   - text2code has been built into ./build/
#   - For each variant, the matching cmake build dir contains the binary
#     (default: build/ for cpu+gpu; pass --build for mpi variants).
#   - bash regenerate.sh <name> has populated grid.bin / pdemodel.txt /
#     datain/ / my_model.hpp for each example.
#
# Exit 0 if all configured examples pass; non-zero otherwise.

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
LIB="$ROOT/apps/library_example"
BASE="$ROOT/baseline"
BUILD="$ROOT/build"

VARIANT=""           # "" → _codegen, gpu → _codegen_gpu, etc.
NP=1                 # rank count for mpi variants
EXAMPLES=()

while [ $# -gt 0 ]; do
    case "$1" in
        --variant)  VARIANT="$2"; shift 2 ;;
        --np)       NP="$2"; shift 2 ;;
        --build)    BUILD="$2"; shift 2 ;;
        --help|-h)
            sed -n '2,30p' "$0"
            exit 0
            ;;
        *)          EXAMPLES+=("$1"); shift ;;
    esac
done

# Resolve BUILD to an absolute path. The harness `cd`s into each
# example's directory before running the binary; relative `--build`
# paths break under that cd.
if [ -d "$BUILD" ]; then
    BUILD="$(cd "$BUILD" && pwd)"
fi

if [ ${#EXAMPLES[@]} -eq 0 ]; then
    EXAMPLES=(poisson2d poisson3d periodic naca0012steady
              lshape isoq3d cone orion
              nsmach8 isoq sharpb2
              naca0012unsteady)
fi

# Map variant → binary suffix and runner prefix.
case "$VARIANT" in
    "")        SUFFIX=""             ; RUNNER=()             ;;
    gpu)       SUFFIX="_gpu"         ; RUNNER=()             ;;
    mpi)       SUFFIX="_mpi"         ; RUNNER=(mpirun -np "$NP") ;;
    mpi_gpu)   SUFFIX="_mpi_gpu"     ; RUNNER=(mpirun -np "$NP") ;;
    *)
        echo "Unknown variant: $VARIANT (expected: gpu | mpi | mpi_gpu)" >&2
        exit 2
        ;;
esac

fail=0

for name in "${EXAMPLES[@]}"; do
    cg_dir="$LIB/${name}_codegen"
    base_dir="$BASE/${name}_serial"
    bin="$BUILD/${name}_codegen${SUFFIX}"

    if [ ! -x "$bin" ]; then
        echo "[SKIP] $name — $bin not built (cmake --build $BUILD --target ${name}_codegen${SUFFIX})"
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

    # The codegen binary is linked against libpdemodelserial.dylib
    # (so the legacy AbiAdapter else-branch can resolve names at link
    # time, even though it's dead code at runtime under M=GeneratedModel).
    # `readbinaryfiles.hpp` calls `cpuInituDriver` non-templated in some
    # branches, which goes through whatever dylib is loaded at runtime.
    # Running text2code on a different app between port and validate
    # leaves the dylib pointing at the wrong PDE → spurious crash.
    # Re-run regenerate.sh before each validate to restore consistency.
    bash "$LIB/regenerate.sh" "$name" >/dev/null

    rm -rf "$cg_dir/dataout"
    mkdir -p "$cg_dir/dataout"
    ( cd "$cg_dir" && "${RUNNER[@]}" "$bin" ./pdeapp.txt >/dev/null )

    if [ -f "$base_dir/outqoi.txt" ]; then
        if diff -q "$cg_dir/dataout/outqoi.txt" "$base_dir/outqoi.txt" >/dev/null; then
            echo "[ OK ] $name (QoI byte-identical)"
        else
            # Codegen splits value+Jacobian into separate struct methods,
            # which prevents SymEngine from sharing CSE replacements
            # (e.g., `n*mu`) between them — legacy lumps both into one
            # `funcjac2cse` call. The value math is equivalent but the
            # FP operation order differs, which can accumulate to a few
            # parts in 1e-5 over a boundary integral. Accept relative
            # diff up to 1e-4 (well above the 1e-6 ~ 1e-8 NewtonTol used
            # in apps/, well below anything that would mask a real bug).
            close="$(python3 - <<PY
import re, sys
def parse(p):
    rows=[]
    for ln in open(p):
        toks=ln.split()
        try: rows.append([float(t) for t in toks])
        except: pass
    return rows
a=parse("$cg_dir/dataout/outqoi.txt")
b=parse("$base_dir/outqoi.txt")
if len(a)!=len(b) or any(len(r)!=len(s) for r,s in zip(a,b)):
    print("shape mismatch"); sys.exit(1)
worst=0.0
for r,s in zip(a,b):
    for x,y in zip(r,s):
        denom=max(abs(x),abs(y),1.0)
        worst=max(worst, abs(x-y)/denom)
if worst > 1e-4:
    print(f"max relative QoI diff = {worst:.3e} (>1e-4)")
    sys.exit(1)
print(f"max relative QoI diff = {worst:.3e}")
PY
)"
            if [ $? -eq 0 ]; then
                echo "[ OK ] $name (QoI numerically close — $close)"
            else
                echo "[FAIL] $name — outqoi.txt differs beyond solver tolerance"
                diff "$cg_dir/dataout/outqoi.txt" "$base_dir/outqoi.txt" || true
                fail=1
            fi
        fi
    elif [ -f "$base_dir/md5.txt" ]; then
        got="$( cd "$cg_dir/dataout" && md5 -r *.bin 2>/dev/null | sort )"
        want="$(sort "$base_dir/md5.txt")"
        if [ "$got" = "$want" ]; then
            echo "[ OK ] $name (bin md5s match)"
        else
            # md5 fail → fall back to a numerical-close check. SymEngine's
            # CSE can reorder floating-point sums (e.g., a+b+c vs b+c+a),
            # which produces last-bit differences that compound through
            # GMRES + Newton iterations. The legacy and codegen paths
            # solve the same equations to within solver tolerance; require
            # relative RMS < 1e-6 across every dataout/*.bin — that's
            # the typical NewtonTol used in apps/, beyond which the
            # solver itself doesn't claim convergence.
            close="$(python3 - <<PY
import os, struct, sys
fail = 0
def read(p):
    with open(p,'rb') as f: d=f.read()
    return struct.unpack(f'{len(d)//8}d', d)
for fn in sorted(os.listdir("$cg_dir/dataout")):
    if not fn.endswith('.bin'): continue
    bp = os.path.join("$base_dir", fn)
    cp = os.path.join("$cg_dir/dataout", fn)
    if not os.path.exists(bp): continue
    a, b = read(bp), read(cp)
    if len(a) != len(b):
        print(f"{fn}: size {len(a)} vs {len(b)}"); fail = 1; continue
    if not a:
        continue
    rms = (sum((x-y)*(x-y) for x,y in zip(a,b))/len(a))**0.5
    nrm = (sum(x*x for x in a)/len(a))**0.5
    rel = rms/max(nrm, 1e-30)
    if rel > 1e-6:
        print(f"{fn}: relrms={rel:.3e} (>1e-6)"); fail = 1
sys.exit(fail)
PY
)"
            if [ -z "$close" ]; then
                echo "[ OK ] $name (numerically close — relrms < 1e-6 across bins)"
            else
                echo "[FAIL] $name — bins differ beyond FP reorder"
                echo "$close"
                diff <(echo "$want") <(echo "$got") || true
                fail=1
            fi
        fi
    else
        echo "[SKIP] $name — baseline has no outqoi.txt or md5.txt"
    fi
done

exit "$fail"
