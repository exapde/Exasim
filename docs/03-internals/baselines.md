# Element-L2 baselines

The test harness compares each run's output to a recorded reference
using a **DG element-L2 relative norm** that is invariant to MPI
partition assignment and tolerant to platform numerics drift (e.g.
Apple Accelerate vs Intel MKL).

```
relative L2 = sqrt( sum_e ||u_run - u_baseline||² / sum_e ||u_baseline||² )
```

Pass threshold: `< 1e-3`. Well above platform-numerics drift (typical
~1e-7), well below anything that would mask a real solver bug
(typical >1e-2).

This file documents the sidecar format, the comparison script, and
how to record a new baseline.

## Properties

- **Partition-invariant.** Aggregate by global element ID and sort
  before comparing. Different `np`, different ParMETIS choices,
  different load balance — same metric.
- **FP-order-tolerant.** Relative L2 with a `1e-3` threshold absorbs
  last-bit drift between BLAS implementations / SIMD widths / FMA
  policies that accumulates through GMRES + Newton iterations.

## The sidecar

When `SaveSolutions` writes `outudg_np<r>.bin`, it also writes
`outelemid_np<r>.bin` — `ne_local` int64 values, one per local
element, holding the **global element ID** of each element in the
order `outudg` stored them.

```
outudg_np<r>.bin       (header[3 doubles]) + (npe × ncu × ne_local) doubles
outelemid_np<r>.bin    (ne_local int64 values)
```

The global IDs come from `mesh.elempart[0..ne_local)`, which is the
same array ParMETIS / DMD uses internally to map local element index
→ global element index.

The runtime function that writes the sidecar lives in
`backend/Solution/solution.hpp` `CSolution<M>::SaveSolutions`.

## The comparator

`apps/library_example/element_l2_diff.py` aggregates and compares.

```bash
python3 apps/library_example/element_l2_diff.py \
    --baseline-dir <path> \
    --current-dir  <path> \
    --np <N>          \
    [--baseline-np <N>]   # default 1 (serial baseline)
    [--rtol 1e-3]         # pass threshold
    [--stem outudg]       # default: outudg
```

Flow:

1. Read all `outudg_np<r>.bin` for `r in [0, np)` and the matching
   `outelemid_np<r>.bin`. Strip the 3-double header.
2. Slice the `outudg` payload into per-element blocks (size = total
   doubles / ne_local).
3. Build a dict `global_id → element_block` for both baseline and
   current.
4. Sort the union of keys; complain on missing pairs.
5. For each common key, accumulate `||a - b||²` and `||a||²` (sum of
   squares of nodal values inside that element block).
6. Return `sqrt(err_sq / ref_sq)` and pass if `< rtol`.

Exit code 0 on pass, non-zero on fail. The summary line is the
relative L2 magnitude, which the test harness echoes:

```
[ OK ] poisson2d (element-L2: outudg: relative element-L2 = 4.13e-12 over 256 elements)
```

The threshold of `1e-3` is conservative — typical lockstep-pass
magnitudes are 1e-12 to 1e-7, real bugs show up at >1e-2.

## Recording a new baseline

For a new example app:

1. Make sure the binary works in serial (run any output-producing
   variant with `np=1`):
   ```bash
   cmake --build build_cpu --target <name>_codegen
   bash apps/library_example/regenerate.sh <name>
   bash apps/library_example/validate_codegen.sh <name>   # FAIL — no baseline
   ```
2. Copy the freshly-produced output bins as the canonical baseline:
   ```bash
   mkdir -p baseline/<name>_serial
   cp apps/library_example/<name>_codegen/dataout_codegen/outudg_np0.bin \
      apps/library_example/<name>_codegen/dataout_codegen/outelemid_np0.bin \
      baseline/<name>_serial/
   # If the model defines QoIs, copy that too — outqoi.txt is
   # partition-invariant by nature and is a stronger gate when
   # available:
   cp apps/library_example/<name>_codegen/dataout_codegen/outqoi.txt \
      baseline/<name>_serial/   # optional
   ```
3. Re-run the validator — should pass with element-L2 ~ 0:
   ```bash
   bash apps/library_example/validate_codegen.sh <name>
   # → [ OK ] <name> (element-L2: outudg: relative element-L2 = 0.000e+00 …)
   ```
4. Commit `baseline/<name>_serial/` (the bin + sidecar files).

For an existing app whose recorded baseline is bit-fragile (drifts
across re-runs even on the same machine), record it once with
element-L2 sidecars and never look back.

## Cross-platform check

After recording on Mac and pushing, on lanka:

```bash
cmake --build build_cpu --target <name>_codegen
ctest --test-dir build_cpu -R "<name>:.*:cpu" --output-on-failure
```

Expected: passes with element-L2 in the 1e-12 to 1e-7 range. If it
gates at 1e-3 (the threshold), the magnitude tells you whether it's
platform-numerics or a real bug — the loose threshold catches real
divergence at 10-100× over the platform-drift baseline.

## What about the legacy bin-md5 path?

Still supported as a fallback. If a baseline directory has only
`outudg_np0.bin` + `md5.txt` (no sidecar), the harness drops to bin-
md5 comparison. New baselines should record the sidecar.

## What about MPI without QoI?

When the baseline lacks `outqoi.txt` AND lacks `outelemid_np0.bin`
(the truly-old baselines), MPI variants fall back to the structural
"each rank wrote its own bins" check — weaker than numerical
validation but catches the common failure modes (binary abort,
missing SaveSolutions, wrong rank count). Once you re-record with
the sidecar, the strong element-L2 check kicks in.
