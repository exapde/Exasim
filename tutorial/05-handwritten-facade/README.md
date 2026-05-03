# 05 — Hand-written Model with the embedded `ExasimSolver<M>` API

This section solves Poisson 2D on the unit square using a hand-
written `Poisson2D` C++ struct, and drives the solver through the
embedded library API. The mesh, boundary conditions, polynomial
order, and physics parameters are set in C++; the converged
solution is read back through `solver.udg()`.

## Files

- `my_model.hpp` defines the `Poisson2D` struct (the same struct
  used in section 04).
- `main.cpp` builds a 16x16 Cartesian quad mesh as flat C arrays,
  hands it to `ExasimSolver<Poisson2D>` via `set_mesh`, registers
  four boundary predicates that tag the four edges of the unit
  square, sets the polynomial order to 3 and the physics
  parameter `mu = 1`, calls `solver.solve()`, and reads the
  converged `udg` array.
- `CMakeLists.txt` builds the binary as
  `tutorial_05_handwritten_facade`.

## Build and run

```bash
cd $EXASIM
cmake --build build --target tutorial_05_handwritten_facade
EXASIM_DIR=$EXASIM $EXASIM/build/tutorial_05_handwritten_facade
```

The binary prints the maximum absolute value of the converged
`udg`. For the manufactured solution `u = sin(πx) sin(πy)`, the
expected value is approximately `3.14158`.

The `EXASIM_DIR` environment variable tells the runtime where to
find the master-element node tables under
`$EXASIM/backend/Preprocessing/`.

For a distributed-memory version of this section that uses
`set_mesh_distributed` instead of `set_mesh`, see
[`../06-handwritten-facade-mpi/`](../06-handwritten-facade-mpi/README.md).
