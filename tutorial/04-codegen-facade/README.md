# 04 — Codegen + embedded `ExasimSolver<M>`

Authoring: SymEngine DSL in `pdemodel.txt` (same as path 03).
Driving: `ExasimSolver<GeneratedModel>::load_pdeapp(…)`. The facade
parses `pdeapp.txt` itself and runs in-memory — no `datain/*.bin`
round-trip.

The generic main lives at `apps/library_example/main_facade.cpp`:

```cpp
exasim::ExasimSolver<GeneratedModel> solver;
solver.load_pdeapp(argv[1], mpirank);
solver.solve(mpiprocs, mpirank);
const double* udg = solver.udg();
```

The same `main_facade.cpp` is reused by all 12 codegen examples;
each `<name>_facade` cmake target compiles it with that example's
`my_model.hpp` on the include path.

## Generate + run

```bash
cd apps/library_example/poisson2d_codegen
bash ../regenerate.sh poisson2d
cmake --build $EXASIM/build_cpu --target poisson2d_facade
$EXASIM/build_cpu/poisson2d_facade ./pdeapp.txt
```

For multi-rank, `pdeapp.txt`'s `gendatain=1` flag must be flipped to
`0` so the binary calls `ParallelPreprocessing` at runtime. The
test harness does this automatically (see `apps/library_example/validate_codegen.sh`);
for direct invocation:

```bash
sed 's/^gendatain *= *1/gendatain = 0/' pdeapp.txt > pdeapp_run.txt
mpirun -np 2 $EXASIM/build_mpi/poisson2d_facade_mpi ./pdeapp_run.txt
```
