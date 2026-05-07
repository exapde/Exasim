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
STEM="codegen"       # codegen | facade — chooses binary family
SKIP_REGEN=0         # 1 → caller already populated grid.bin / datain
QUICK=0              # 1 → run a small fast subset (HOT.7.14)
EXAMPLES=()

while [ $# -gt 0 ]; do
    case "$1" in
        --variant)    VARIANT="$2"; shift 2 ;;
        --np)         NP="$2"; shift 2 ;;
        --build)      BUILD="$2"; shift 2 ;;
        --facade)     STEM="facade"; shift ;;
        --no-regen)   SKIP_REGEN=1; shift ;;
        --quick)      QUICK=1; shift ;;
        --help|-h)
            sed -n '2,30p' "$0"
            exit 0
            ;;
        *)            EXAMPLES+=("$1"); shift ;;
    esac
done

# Resolve BUILD to an absolute path. The harness `cd`s into each
# example's directory before running the binary; relative `--build`
# paths break under that cd.
if [ -d "$BUILD" ]; then
    BUILD="$(cd "$BUILD" && pwd)"
fi

if [ ${#EXAMPLES[@]} -eq 0 ]; then
    if [ "$QUICK" = "1" ]; then
        # --quick subset for inner-loop iteration: small fast Poissons
        # + a couple of heavier shapes that exercise NS + boundary
        # logic. Skips naca0012unsteady (5min outlier), the slow
        # lshape/cone/orion/isoq* runs, and sharpb2.
        EXAMPLES=(poisson2d periodic naca0012steady)
    else
        EXAMPLES=(poisson2d poisson3d periodic naca0012steady
                  lshape isoq3d cone orion
                  nsmach8 isoq sharpb2
                  naca0012unsteady)
    fi
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
    bin="$BUILD/${name}_${STEM}${SUFFIX}"

    if [ ! -x "$bin" ]; then
        echo "[SKIP] $name — $bin not built (cmake --build $BUILD --target ${name}_${STEM}${SUFFIX})"
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
    if [ "$SKIP_REGEN" = "0" ]; then
        bash "$LIB/regenerate.sh" "$name" >/dev/null
    fi

    # Per-gate dataout dir so codegen:X and facade:X (or different
    # variants) can run in parallel without colliding. The binary's
    # default dataoutpath is "$datapath/dataout"; override per-gate
    # via a sed'd pdeapp_run_<stem>.txt that names a unique subdir.
    out_subdir="dataout_${STEM}${SUFFIX}"
    out_path="$cg_dir/$out_subdir"
    rm -rf "$out_path"
    mkdir -p "$out_path"

    # Build the per-gate pdeapp_run file:
    #   - drop any existing dataoutpath line
    #   - for MPI variants: drop any existing gendatain line entirely,
    #     then append `gendatain = 0;` below. The struct default is 1,
    #     so a sed-replace alone misses pdeapp.txt files that omit the
    #     line entirely (the binary then skips ParallelPreprocessing
    #     and dies on missing meshN.bin). Drop+append is robust to
    #     either spelling.
    #   - append our per-gate dataoutpath
    pdeapp_run="$cg_dir/pdeapp_run_${STEM}${SUFFIX}.txt"
    sed_script='/^dataoutpath *=/d'
    if [ "$VARIANT" = "mpi" ] || [ "$VARIANT" = "mpi_gpu" ]; then
        sed_script="$sed_script
/^gendatain *=/d"
    fi
    sed "$sed_script" "$cg_dir/pdeapp.txt" > "$pdeapp_run"
    if [ "$VARIANT" = "mpi" ] || [ "$VARIANT" = "mpi_gpu" ]; then
        echo "gendatain = 0;" >> "$pdeapp_run"
    fi
    echo "dataoutpath = \"$out_subdir\";" >> "$pdeapp_run"
    pdeapp_in="./$(basename "$pdeapp_run")"

    # Run binary. Don't fail on non-zero exit alone — Kokkos / CUDA
    # cleanup code on lanka segfaults during process exit even after
    # the solver wrote correct dataout. Fall through to the dataout
    # check below: empty dataout = real failure, non-empty = pass.
    ( cd "$cg_dir" && "${RUNNER[@]}" "$bin" "$pdeapp_in" >/dev/null ) || true

    # Empty dataout = solver didn't run (e.g. mpirun -np 2 with no
    # datain/meshN.bin files, or programmatic crash before SaveSolutions).
    # Without this guard the bin-md5 / numerical fallback below
    # silently passes on zero files iterated and prints "[ OK ]".
    if [ -z "$(ls -A "$out_path" 2>/dev/null)" ]; then
        echo "[FAIL] $name — binary produced no dataout/ files"
        fail=1
        continue
    fi

    # Comparison strategy:
    #   1. element-L2 (HOT.7.18, preferred when both baseline + current
    #      have outelemid_np<r>.bin sidecars; partition- and FP-order-
    #      invariant — see apps/library_example/element_l2_diff.py).
    #   2. outqoi.txt diff (partition-invariant scalar QoI integrals).
    #   3. bin-md5 / numerical fallback (serial only; FP-fragile).
    #   4. "ran-cleanly" structural check (MPI without outqoi baseline).
    is_mpi=0
    if [ "$VARIANT" = "mpi" ] || [ "$VARIANT" = "mpi_gpu" ]; then
        is_mpi=1
    fi

    # Element-L2 takes priority when both sides have the elemid sidecar.
    if [ -f "$base_dir/outelemid_np0.bin" ] && [ -f "$out_path/outelemid_np0.bin" ]; then
        diff_out="$(python3 "$LIB/element_l2_diff.py" \
            --baseline-dir "$base_dir" \
            --current-dir  "$out_path" \
            --np "$NP" 2>&1)"
        diff_status=$?
        if [ $diff_status -eq 0 ]; then
            echo "[ OK ] $name (element-L2: $diff_out)"
        else
            echo "[FAIL] $name — element-L2: $diff_out"
            fail=1
        fi
        continue
    fi

    if [ "$is_mpi" = "1" ] && [ ! -f "$base_dir/outqoi.txt" ]; then
        # No partition-invariant baseline: just verify each rank wrote
        # its expected output bins. This is weaker than numerical
        # validation but catches the common failure modes (binary abort,
        # missing SaveSolutions, wrong rank count).
        missing=""
        for ((r=0; r<NP; r++)); do
            for stem in outudg outuhat; do
                f="$out_path/${stem}_np${r}.bin"
                if [ ! -s "$f" ]; then missing="$missing $stem _np${r}.bin"; fi
            done
        done
        if [ -z "$missing" ]; then
            echo "[ OK ] $name (ran cleanly, $NP ranks wrote bins; no QoI baseline for cross-rank check)"
        else
            echo "[FAIL] $name — missing per-rank bins:$missing"
            fail=1
        fi
        continue
    fi

    if [ -f "$base_dir/outqoi.txt" ]; then
        if diff -q "$out_path/outqoi.txt" "$base_dir/outqoi.txt" >/dev/null; then
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
a=parse("$out_path/outqoi.txt")
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
                diff "$out_path/outqoi.txt" "$base_dir/outqoi.txt" || true
                fail=1
            fi
        fi
    elif [ -f "$base_dir/md5.txt" ]; then
        got="$( cd "$out_path" && md5 -r *.bin 2>/dev/null | sort )"
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
for fn in sorted(os.listdir("$out_path")):
    if not fn.endswith('.bin'): continue
    bp = os.path.join("$base_dir", fn)
    cp = os.path.join("$out_path", fn)
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
