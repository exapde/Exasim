# 03 — Codegen Model with the embedded `ExasimSolver<M>`

This section solves Poisson 2D on the unit square using the
codegen authoring path together with the embedded library API.
`text2code` generates the model header from `pdemodel.txt`;
`main.cpp` constructs an `ExasimSolver<GeneratedModel>`, calls
`load_pdeapp(argv[1])` to parse `pdeapp.txt`, and runs the
solver.

## Files

- `pdemodel.txt` is the SymEngine DSL describing the PDE
  (identical to the file in section 02).
- `main.cpp` constructs the embedded solver, calls
  `load_pdeapp`, and invokes `solve()`. Multi-rank builds bind
  one rank per device using `MPI_COMM_TYPE_SHARED`.
- `pdeapp.txt` is the runtime configuration that `load_pdeapp`
  parses.
- `grid.bin` is the binary mesh referenced by `pdeapp.txt`.
- `CMakeLists.txt` builds `tutorial_03_codegen_facade`. The build
  consumes the codegen output from this directory.

## Build and run

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-codegen-facade \
    $EXASIM/tutorial/03-codegen-facade/pdeapp.txt
cmake --build build --target tutorial_03_codegen_facade
cd $EXASIM/tutorial/03-codegen-facade
mkdir -p datain dataout
$EXASIM/build/tutorial_03_codegen_facade ./pdeapp.txt
```

For multi-rank builds, prepend `mpirun -np N` and use the matching
`build_mpi` or `build_mpi_gpu` directory.
