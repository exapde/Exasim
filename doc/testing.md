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
   Currently disabled — see [§ CI status](#ci-status).

This doc covers all three.

---

## 1. The validation script

```bash
bash apps/library_example/validate_codegen.sh \
    [--variant gpu|mpi|mpi_gpu] \
    [--np N]                  \
    [--build <build-dir>]     \
    [--facade]                \
    [<example> ...]
```

| flag           | default                          | meaning                                                  |
| -------------- | -------------------------------- | -------------------------------------------------------- |
| `--variant`    | none → bin suffix `""` (cpu)     | `gpu`/`mpi`/`mpi_gpu` → adds matching suffix             |
| `--np`         | `1`                              | rank count for the `mpi*` variants                       |
| `--build`      | `./build`                        | directory holding the binary                             |
| `--facade`     | off (uses `_codegen` family)     | switches binary stem to `_facade` (HOT.7.9)              |
| `<example>...` | all 12 codegen examples          | pass one or more to validate a subset                    |

What it does, per example:

1. Runs `regenerate.sh <example>` to make sure `grid.bin`,
   `pdemodel.txt`, `datain/`, and `my_model.hpp` are present and
   consistent with the source app.
2. Clears `<example>_codegen/dataout/`, runs the binary with the
   committed `pdeapp.txt` from that directory.
3. Compares `dataout/outqoi.txt` to
   `baseline/<example>_serial/outqoi.txt` if it exists; otherwise
   compares md5s of all `*.bin`; otherwise falls back to "rel-RMS
   over all bins < 1e-6" (the typical solver tolerance).

Output: one line per example, `[ OK ] <name>` or `[FAIL] <name>`.
Exit 0 if every gate passed; non-zero otherwise.

### Examples

```bash
# Single example, single variant:
bash apps/library_example/validate_codegen.sh poisson2d
bash apps/library_example/validate_codegen.sh --variant gpu --build build_gpu poisson2d

# Full matrix on lanka-style four-build layout (legacy codegen path):
bash apps/library_example/validate_codegen.sh                        --build build_cpu
bash apps/library_example/validate_codegen.sh --variant gpu          --build build_gpu
bash apps/library_example/validate_codegen.sh --variant mpi --np 2   --build build_mpi
bash apps/library_example/validate_codegen.sh --variant mpi_gpu --np 2 --build build_mpi_gpu

# Same matrix, but through the ExasimSolver façade:
bash apps/library_example/validate_codegen.sh --facade                --build build_cpu
bash apps/library_example/validate_codegen.sh --facade --variant gpu  --build build_gpu
# … etc.
```

### Tolerances

The harness gates byte-identical first; failing that, two numerical
fallbacks kick in:

- For `outqoi.txt`: rel diff per cell ≤ 1e-4 (well above NewtonTol of
  1e-6 ~ 1e-8 used in apps).
- For `*.bin`: rel-RMS per file ≤ 1e-6.

These bounds are deliberately loose — SymEngine's CSE can reorder
floating-point sums (`a+b+c` vs `b+c+a`) which compounds through
GMRES + Newton. The bounds catch real bugs; they don't catch FP
ordering drift.

### Failure modes

| symptom                                  | usually means                                                          |
| ---------------------------------------- | ---------------------------------------------------------------------- |
| `[SKIP] <name> — <bin> not built`        | the corresponding cmake target wasn't built. run `cmake --build`.      |
| `[SKIP] <name> — <dir> not populated`    | `regenerate.sh` failed; check that text2code is built (`build/text2code`). |
| `[SKIP] <name> — no baseline`            | new example without a recorded baseline. run a recording pass first.   |
| `[FAIL] <name> — outqoi.txt differs`     | numerical drift past 1e-4 — bisect against last green commit.          |
| segfault during run                      | usually a struct-init bug in HOT.7.x; `gdb -batch -ex run -ex bt` it.  |

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

`ctest` honors the same dependency wiring as the validate script:
each test invokes `regenerate.sh` first, then the binary, then the
diff. The TIMEOUT is 600s per example (long enough for naca0012's
Newton+GMRES inner loop on slow runners).

The set of test names a build dir registers depends on which backend
variant cmake was configured with:

| build dir         | configure flags                              | tests registered            |
| ----------------- | -------------------------------------------- | --------------------------- |
| `build_cpu/`      | `EXASIM_NOMPI=ON, EXASIM_CUDA=OFF`           | `*:*:cpu`                   |
| `build_gpu/`      | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON`            | `*:*:gpu`                   |
| `build_mpi/`      | `EXASIM_MPI=ON, EXASIM_NOMPI=OFF`            | `*:*:mpi` (np=2)            |
| `build_mpi_gpu/`  | `EXASIM_MPI=ON, EXASIM_CUDA=ON`              | `*:*:mpi_gpu` (np=2)        |

For each example × variant, two tests are registered:
`codegen:<ex>:<variant>` (Family A from `doc/running_apps.md`) and
`facade:<ex>:<variant>` (Family B). 12 examples × 2 families = 24
tests per build dir.

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
   cp apps/library_example/<name>_codegen/dataout/outqoi.txt \
      baseline/<name>_serial/outqoi.txt
   bash apps/library_example/validate_codegen.sh <name>   # OK
   ```
5. (Optional but recommended) commit md5s of the bin files for
   md5-fallback gating:
   ```bash
   cd apps/library_example/<name>_codegen/dataout
   md5 -r *.bin > /Users/teoc/projects/exasim2/Exasim/baseline/<name>_serial/md5.txt
   ```

### Recording a new baseline

Use the legacy `cput2cEXASIM` runner on the same `pdeapp.txt`:

```bash
bash baseline/record.sh <name>
```

(See `baseline/README.md` for the full procedure, including the
"text2code-rerun-on-poisson2d" reset that prevents `libpdemodelserial`
from pointing at the wrong app's kernels.)

---

## 3. CI status

`/.github/workflows/poisson2d-cpu.yml` is the lone CI workflow. It
builds Kokkos serial + SymEngine + METIS + text2code + Exasim CPU
from source on a fresh Ubuntu 22.04 runner and runs `poisson2d_template`
as a smoke. Heavy deps are cached by source SHA.

**Status: disabled** as of HOT.7.9. The GitHub-hosted runner has been
flaky during the SymEngine + Kokkos build (intermittent OOM and
unexplained hangs). The `workflow_dispatch` trigger is preserved so
it can be invoked manually for diagnosis, and the job is gated by
`if: false`. Steps to re-enable once the runner is sorted out:

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

The fast iteration loop:

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

When a regression is suspected:

```bash
# Bisect candidate range:
git bisect start HEAD <good-sha>
git bisect run bash apps/library_example/validate_codegen.sh --facade poisson2d
```

Cross-arch validation (lanka):

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
ctest --test-dir build_cpu     -R poisson2d --output-on-failure
ctest --test-dir build_gpu     -R poisson2d --output-on-failure
ctest --test-dir build_mpi     -R poisson2d --output-on-failure
ctest --test-dir build_mpi_gpu -R poisson2d --output-on-failure
```
