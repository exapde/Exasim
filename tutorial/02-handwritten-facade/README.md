# 02 — Hand-written model + embedded `ExasimSolver<M>`

Authoring: hand-written C++ struct `Poisson2D` in
`apps/library_example/solve_square/solve_square_facade.cpp`.
Driving: `ExasimSolver<Poisson2D>` — mesh, BCs, and parameters all
set in C++; the converged solution comes back through `solver.udg()`.

The full source is `apps/library_example/solve_square/solve_square_facade.cpp`.
The relevant body:

```cpp
exasim::ExasimSolver<Poisson2D> solver;
solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4);
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[1])     < 1e-8; });
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[0] - 1) < 1e-8; });
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[1] - 1) < 1e-8; });
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[0])     < 1e-8; });
solver.set_polynomial_order(3);
solver.set_quadrature_order(6);
solver.set_physics_params({1.0});
solver.solve();
const double* udg = solver.udg();
```

No `pdeapp.txt`. No `datain/*.bin`. No `dataout/*.bin`. The mesh is
built as flat C arrays in the program; the solution comes out via
`solver.udg()`.

## Run

```bash
cmake --build build_cpu --target solve_square_facade
$EXASIM/build_cpu/solve_square_facade
```

GPU and MPI variants:

```bash
$EXASIM/build_gpu/solve_square_facade_gpu
mpirun -np 2 $EXASIM/build_mpi/solve_square_facade_mpi
mpirun -np 2 $EXASIM/build_mpi_gpu/solve_square_facade_mpi_gpu
```

The MPI variants use `solver.set_mesh_distributed(…)` so each rank
provides its own slice of the global mesh. ParMETIS partitions for
load balance internally. See
[`../../docs/03-methods/embedded-facade.md`](../../docs/03-methods/embedded-facade.md)
§2 for the multi-rank pattern.
