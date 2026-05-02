# Catalog: how to run an Exasim app

Every app has at most three independent dimensions:

|                 | Author the model               | Author the runtime config        | Author the mesh                  |
| --------------- | ------------------------------ | -------------------------------- | -------------------------------- |
| Hand-written    | `my_model.hpp` (CRTP)          | C++ struct (PDE / InputParams)   | C++ flat arrays (`set_mesh*`)    |
| Convenience     | `pdemodel.txt` → `text2code`   | `pdeapp.txt`                     | `grid.bin`/`grid.txt` on disk    |

Each dimension is independent — you can mix authoring styles. Below is
the catalog of binaries Exasim ships and the command line for each.

For each "family" there are four backend variants:

| suffix             | what                                  |
| ------------------ | ------------------------------------- |
| (none)             | single-rank serial CPU                |
| `_gpu`             | single-rank CUDA or HIP               |
| `_mpi`             | multi-rank serial CPU                 |
| `_mpi_gpu`         | multi-rank CUDA or HIP                |

Which backend a binary supports is decided at *configure* time
(`-DEXASIM_CUDA=ON`, `-DEXASIM_MPI=ON`, etc.) — the four variants live
in four CMake build dirs (`build_cpu/`, `build_gpu/`, `build_mpi/`,
`build_mpi_gpu/`).

---

## Family A — codegen + `<exasim/run.hpp>` (legacy)

The legacy path. `text2code` emits `my_model.hpp` from
`pdemodel.txt`, and the binary's `main.cpp` is a 3-liner that calls
`exasim::run<GeneratedModel>(argc, argv)`. Runtime configuration
comes from a `pdeapp.txt` text file.

```
build_cpu/poisson2d_codegen           ./pdeapp.txt
build_gpu/poisson2d_codegen_gpu       ./pdeapp.txt
mpirun -np N build_mpi/poisson2d_codegen_mpi          ./pdeapp.txt
mpirun -np N build_mpi_gpu/poisson2d_codegen_mpi_gpu  ./pdeapp.txt
```

Outputs go to `dataout/outudg_np<r>.bin` etc., per the `pdeapp.txt`'s
`dataoutpath`. The runtime preprocessor writes intermediate
`datain/{app,master,mesh,sol}.bin` per rank (always — the file ABI is
load-bearing here).

12 codegen examples ship: `poisson2d`, `poisson3d`, `periodic`,
`naca0012steady`, `naca0012unsteady`, `lshape`, `isoq3d`, `isoq`,
`cone`, `orion`, `nsmach8`, `sharpb2`. Each gets four binaries.

## Family B — codegen + `ExasimSolver<M>` façade (HOT.7.9)

Same `my_model.hpp` and `pdeapp.txt` as Family A, but the binary's
main is `apps/library_example/main_facade.cpp` which routes through
`ExasimSolver<GeneratedModel>::load_pdeapp(filename) + .solve(...)`.
For tests, the externally-visible behavior is identical to Family A
(reads `pdeapp.txt`, writes `dataout/`).

```
build_cpu/poisson2d_facade            ./pdeapp.txt
build_gpu/poisson2d_facade_gpu        ./pdeapp.txt
mpirun -np N build_mpi/poisson2d_facade_mpi          ./pdeapp.txt
mpirun -np N build_mpi_gpu/poisson2d_facade_mpi_gpu  ./pdeapp.txt
```

Same 12 examples as Family A.

## Family C — hand-written model + façade (in-memory)

`apps/library_example/solve_square/solve_square_facade.cpp` is the
canonical example. The model, mesh (`set_mesh` or
`set_mesh_distributed`), boundary classifiers (`add_boundary` lambdas),
and physics params live entirely in C++. No `pdeapp.txt`, no
`pdemodel.txt`, no `grid.bin`, no `dataout/` writes.

```
build_cpu/solve_square_facade
build_gpu/solve_square_facade_gpu
mpirun -np N build_mpi/solve_square_facade_mpi
mpirun -np N build_mpi_gpu/solve_square_facade_mpi_gpu
```

Each takes no command-line arguments — the problem is hard-coded in
the source. Output is read from `solver.udg()` after `solver.solve()`.

A second variant `solve_square` (without `_facade` suffix) drives the
same problem manually via `CPreprocessing::take()` +
`CSolution<M>(Preprocessed&&, ...)` — useful for readers who want to
see what the façade hides.

## Family D — hand-written model + `<exasim/run.hpp>` + pdeapp.txt

Hand-written `my_model.hpp` but using the legacy file-driven runtime.
`apps/library_example/poisson2d/main.cpp` is the canonical example.

```
build_cpu/poisson2d_template          ./pdeapp.txt
build_gpu/poisson2d_template_gpu      ./pdeapp.txt
mpirun -np N build_mpi/poisson2d_template_mpi          ./pdeapp.txt
```

(No `_mpi_gpu` variant for poisson2d_template today; same shape as A
otherwise.)

---

## Validation

`apps/library_example/validate_codegen.sh` runs each binary and
compares output against `baseline/<name>_serial/`:

```
# Family A (codegen):
bash apps/library_example/validate_codegen.sh                              # cpu
bash apps/library_example/validate_codegen.sh --variant gpu      --build build_gpu
bash apps/library_example/validate_codegen.sh --variant mpi      --build build_mpi      --np 2
bash apps/library_example/validate_codegen.sh --variant mpi_gpu  --build build_mpi_gpu  --np 2

# Family B (facade) — same flags + --facade:
bash apps/library_example/validate_codegen.sh --facade
bash apps/library_example/validate_codegen.sh --facade --variant gpu      --build build_gpu
bash apps/library_example/validate_codegen.sh --facade --variant mpi --np 2 --build build_mpi
bash apps/library_example/validate_codegen.sh --facade --variant mpi_gpu --np 2 --build build_mpi_gpu

# Single example:
bash apps/library_example/validate_codegen.sh --facade poisson2d
```

The harness gates: `outqoi.txt` byte-identical (preferred), or
`md5(*.bin)` match, or rel-RMS of bin contents < 1e-6 (numerical
fallback).

---

## CTest integration

`CMakeLists.txt` registers the same gates as CTest tests:

```
cmake --build build_cpu
ctest --test-dir build_cpu              # runs every Family A + B test
ctest --test-dir build_cpu -R poisson2d # filter by name
ctest --test-dir build_cpu -j16         # parallel
ctest --test-dir build_mpi --output-on-failure
```

Each test wraps the binary invocation with the matching baseline diff.
The test names follow `<family>:<example>:<backend>`, e.g.,
`facade:poisson2d:mpi`.

For binaries that need preprocessed inputs (`datain/`, `my_model.hpp`),
each test depends on a `regenerate_<example>` fixture that runs
`regenerate.sh` once.
