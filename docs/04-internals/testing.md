# Testing

Three layers of tests live in this repo, in increasing scope:

1. **Unit-style:** `apps/library_example/validate_codegen.sh` invokes
   one binary with one example and gates against the matching
   `baseline/<example>_serial/`.
2. **Backend matrix:** four configure-time variants (cpu / gpu / mpi /
   mpi+gpu) of the same binary family, registered as CTest tests so
   `ctest -j` runs them all.
3. **CI smoke:** `.github/workflows/poisson2d-cpu.yml` runs a single
   hand-written-model example end-to-end on a fresh Ubuntu runner.
   Currently disabled — see [§ CI status](#3-ci-status).

This doc covers all three.

---

## 1. The validation script

```bash
bash apps/library_example/validate_codegen.sh \
    [--variant gpu|mpi|mpi_gpu] \
    [--np N]                  \
    [--build <build-dir>]     \
    [--facade]                \
    [--quick]                 \
    [--no-regen]              \
    [<example> ...]
```

| flag           | default                          | meaning                                                  |
| -------------- | -------------------------------- | -------------------------------------------------------- |
| `--variant`    | none → bin suffix `""` (cpu)     | `gpu`/`mpi`/`mpi_gpu` → adds matching suffix             |
| `--np`         | `1`                              | rank count for the `mpi*` variants                       |
| `--build`      | `./build`                        | directory holding the binary                             |
| `--facade`     | off (uses `_codegen` family)     | switches binary stem to `_facade` (the embedded library API path) |
| `--quick`      | off (full 12-example matrix)     | runs the small fast subset (poisson2d, periodic, naca0012steady) |
| `--no-regen`   | off (re-runs `regenerate.sh`)    | skip text2code regen — caller must have done it already  |
| `<example>...` | all 12 codegen examples          | pass one or more to validate a subset                    |

What it does, per example:

1. (Unless `--no-regen`) runs `regenerate.sh <example>` to make sure
   `grid.bin`, `pdemodel.txt`, `datain/`, and `my_model.hpp` are
   present and consistent with the source app.
2. Sets up a per-gate `dataout_<stem><variant>/` directory (so
   codegen and facade gates for the same example can run in parallel
   without clobbering each other's output).
3. For MPI variants: rewrites `pdeapp.txt` so `gendatain=0`, forcing
   the binary to call `ParallelPreprocessing` at runtime (the
   committed `pdeapp.txt` has `gendatain=1` for serial-fast paths).
4. Runs the binary; if it produces no output files, fail.
5. Compares against the baseline using the comparison priority below.

Output: one line per example, `[ OK ] <name>` or `[FAIL] <name>`.
Exit 0 if every gate passed; non-zero otherwise.

### Examples

```bash
# Single example, single variant:
bash apps/library_example/validate_codegen.sh poisson2d
bash apps/library_example/validate_codegen.sh --variant gpu --build build_gpu poisson2d

# Full matrix (legacy codegen path):
bash apps/library_example/validate_codegen.sh                           --build build_cpu
bash apps/library_example/validate_codegen.sh --variant gpu             --build build_gpu
bash apps/library_example/validate_codegen.sh --variant mpi --np 2      --build build_mpi
bash apps/library_example/validate_codegen.sh --variant mpi_gpu --np 2  --build build_mpi_gpu

# Same matrix, but through the ExasimSolver façade:
bash apps/library_example/validate_codegen.sh --facade                  --build build_cpu
bash apps/library_example/validate_codegen.sh --facade --variant gpu    --build build_gpu

# Inner-loop iteration — small fast subset:
bash apps/library_example/validate_codegen.sh --quick
```

### Comparison priority

The harness chooses one of four comparison strategies, in order of
preference, based on which artifacts the baseline contains:

1. **Element-L2** (preferred). If both baseline and current dataout
   contain `outelemid_np<r>.bin` sidecars, aggregate `outudg`
   element-by-element across ranks and compute a relative L2 norm.
   Partition-invariant and FP-order-tolerant. Threshold: `< 1e-3`.
   See [`baselines.md`](baselines.md) for details.
2. **outqoi.txt** (good fallback for variant invariance). Domain or
   boundary integrals are partition-invariant and FP-stable to ~1e-8;
   compare byte-identical first, then numerically with rel-diff
   threshold `< 1e-4` per cell.
3. **bin md5 / numerical RMS** (serial-only fallback). md5-compare
   `outudg_np0.bin` and friends; if the digest doesn't match, fall
   back to per-file rel-RMS `< 1e-6`. Sensitive to FP order — fragile
   across platforms.
4. **Ran-cleanly** (last resort, MPI variants without QoI baseline).
   Verify each rank in `[0, NP)` produced its `outudg/outuhat_np<r>.bin`.
   Catches binary aborts and missing SaveSolutions calls; doesn't
   gate numerics.

The element-L2 path is the architectural answer to bin-md5
fragility. See [`baselines.md`](baselines.md) for the sidecar format
and how to record a new baseline.

### Failure modes

| symptom                                    | usually means                                                          |
| ------------------------------------------ | ---------------------------------------------------------------------- |
| `[SKIP] <name> — <bin> not built`          | the corresponding cmake target wasn't built. run `cmake --build`.      |
| `[SKIP] <name> — <dir> not populated`     | `regenerate.sh` failed; check that text2code is built (`build/text2code`). |
| `[SKIP] <name> — no baseline`             | new example without a recorded baseline. run a recording pass first.   |
| `[FAIL] <name> — element-L2 = X (>1e-3)`   | numerical drift past threshold — bisect against last green commit. magnitude tells you whether it's platform numerics (e-3 to e-7) or a real bug (>1e-2). |
| `[FAIL] <name> — outqoi.txt differs`       | scalar QoI drift past 1e-4. usually a real boundary or source bug.     |
| `[FAIL] <name> — binary produced no dataout/ files` | the binary aborted before SaveSolutions ran.                  |
| segfault during run                        | likely a struct-init bug or library version mismatch — `gdb -batch -ex run -ex bt` it. |

When `codegen:<example>:<variant>` and `facade:<example>:<variant>`
both fail with the **same magnitude**, the divergence is not facade-
specific; both stems are agreeing on the same wrong answer. Most
common cause: platform-numerics drift vs the recorded baseline. If
only one stem fails, suspect a real regression in that path.

---

## 2. CTest integration

Every cmake build dir registers the same gates as CTest tests. The
test names are `<family>:<example>:<variant>`, e.g.
`facade:poisson2d:cpu` or `codegen:naca0012steady:gpu`.

```bash
# Run every test for the current build dir's variant:
ctest --test-dir build_cpu --output-on-failure -j16

# Filter by name (glob):
ctest --test-dir build_cpu -R poisson2d
ctest --test-dir build_cpu -R 'facade:.*:cpu'

# Show what would run, don't run:
ctest --test-dir build_cpu -N

# Re-run only the tests that failed last time:
ctest --test-dir build_cpu --rerun-failed
```

The tests parallelize cleanly with `-j8` because text2code's output
goes to per-example dirs, so different examples never race on
`backend/Model/libpdemodel*.{so,dylib}`. The codegen and facade
gates for the same example also run in parallel because each writes
to its own per-stem `dataout_<stem><variant>/` directory.

The TIMEOUT is 1800s per gate (large NS meshes like cone or sharpb2
under -j8 contention can take 10+ minutes per Newton+GMRES solve).

The set of test names a build dir registers depends on which backend
variant cmake was configured with:

| build dir         | configure flags                              | tests registered            |
| ----------------- | -------------------------------------------- | --------------------------- |
| `build_cpu/`      | `EXASIM_NOMPI=ON, EXASIM_CUDA=OFF`           | `*:*:cpu`                   |
| `build_gpu/`      | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON`            | `*:*:gpu`                   |
| `build_mpi/`      | `EXASIM_MPI=ON, EXASIM_NOMPI=OFF`            | `*:*:mpi` (np=2)            |
| `build_mpi_gpu/`  | `EXASIM_MPI=ON, EXASIM_CUDA=ON`              | `*:*:mpi_gpu` (np=2)        |

For each example × variant, two tests are registered:
`codegen:<ex>:<variant>` (legacy `<exasim/run.hpp>` driver) and
`facade:<ex>:<variant>` (embedded `ExasimSolver<M>` driver).
12 examples × 2 families = 24 gates per build dir, plus per-example
regen and per-target build fixtures.

### CTest fixture topology

```
_build_text2code              (FIXTURES_SETUP    text2code)
  ↓
_regen_<example>              (FIXTURES_SETUP    regen_<example>)
  ↓
_build_<example>_<binary>     (FIXTURES_SETUP    build_<example>_<binary>)
  ↓
<family>:<example>:<variant>  (FIXTURES_REQUIRED both above)
```

Each gate runs `validate_codegen.sh --no-regen <example>` so the
gate itself only runs the binary and diffs.

### Adding a new example

1. Drop a `<name>_codegen/` directory under `apps/library_example/`
   with `pdeapp.txt`, `main.cpp`, and (committed alongside the source
   app) `regenerate.sh`-compatible inputs.
2. Add `add_codegen_example(<name>)` to `install/CMakeLists.txt`.
3. Add `<name>` to the `_codegen_examples` list at the bottom of the
   same file (the CTest registration loop).
4. Record the baseline:
   ```bash
   cmake --build build_cpu --target <name>_codegen
   bash apps/library_example/validate_codegen.sh <name>   # FAIL — no baseline
   cp apps/library_example/<name>_codegen/dataout_codegen/outudg_np0.bin \
      apps/library_example/<name>_codegen/dataout_codegen/outelemid_np0.bin \
      baseline/<name>_serial/
   # If the model defines QoIs, also copy outqoi.txt:
   cp apps/library_example/<name>_codegen/dataout_codegen/outqoi.txt \
      baseline/<name>_serial/   # optional but partition-invariant
   bash apps/library_example/validate_codegen.sh <name>   # OK
   ```

See [`baselines.md`](baselines.md) for the full sidecar format.

---

## 3. CI status

`.github/workflows/poisson2d-cpu.yml` is the lone CI workflow. It
builds Kokkos serial + SymEngine + METIS + text2code + Exasim CPU
from source on a fresh Ubuntu 22.04 runner and runs `poisson2d_template`
as a smoke. Heavy deps are cached by source SHA.

**Status: disabled.** The GitHub-hosted runner has been flaky during
the SymEngine + Kokkos build (intermittent OOM and unexplained
hangs). The `workflow_dispatch` trigger is preserved so it can be
invoked manually for diagnosis, and the job is gated by `if: false`.
Steps to re-enable once the runner is sorted out:

1. Open `.github/workflows/poisson2d-cpu.yml`.
2. Remove the `if: false` line on the `build-and-run` job.
3. Restore the `push:` and `pull_request:` triggers (currently
   commented out).
4. Push — CI fires on the next commit.

A larger CI matrix (gpu / mpi / mpi+gpu via self-hosted runners or
container images) is queued but not yet built — the immediate
priority is stabilizing the existing CPU smoke first.

---

## 4. Local test loop

Inner-loop iteration (a few seconds per cycle):

```bash
# Edit a backend file …
# Build the affected target(s):
cmake --build build_cpu --target poisson2d_codegen poisson2d_facade

# Validate just that example:
bash apps/library_example/validate_codegen.sh poisson2d
bash apps/library_example/validate_codegen.sh --facade poisson2d

# … or via ctest:
ctest --test-dir build_cpu -R poisson2d --output-on-failure
```

Quick cross-section (~30s on Mac):

```bash
bash apps/library_example/validate_codegen.sh --quick
```

When a regression is suspected:

```bash
# Bisect candidate range:
git bisect start HEAD <good-sha>
git bisect run bash apps/library_example/validate_codegen.sh --facade poisson2d
```

Cross-arch validation (the four-build layout):

```bash
# Build all four variants in parallel:
PATH=/usr/local/cuda/bin:/path/to/openmpi/bin:$PATH \
    cmake --build build_cpu     -j16 &
PATH=/usr/local/cuda/bin:/path/to/openmpi/bin:$PATH \
    cmake --build build_gpu     -j8  &
PATH=/usr/local/cuda/bin:/path/to/openmpi/bin:$PATH \
    cmake --build build_mpi     -j16 &
PATH=/usr/local/cuda/bin:/path/to/openmpi/bin:$PATH \
    cmake --build build_mpi_gpu -j8  &
wait

# Validate each:
ctest --test-dir build_cpu     -j 8
ctest --test-dir build_gpu     -j 8
ctest --test-dir build_mpi     -j 8
ctest --test-dir build_mpi_gpu -j 8
```

The four ctests must run sequentially because they share
`apps/library_example/<name>_codegen/` directories for text2code
regen.
