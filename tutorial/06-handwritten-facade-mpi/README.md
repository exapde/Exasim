# 06 — Hand-written Model with `ExasimSolver<M>::set_mesh_distributed`

This section solves Poisson 2D on the unit square through the
embedded library API in distributed-memory mode. Each rank builds
its own slice of the global Cartesian quad mesh, hands it to
`solver.set_mesh_distributed`, and calls
`solver.solve(mpiprocs, mpirank)`. ParMETIS migrates the slice and
repartitions for load balance internally.

## Files

- `my_model.hpp` defines the `Poisson2D` struct (the same struct
  used in section 04 and section 05).
- `main.cpp` initializes MPI, computes a contiguous slice of the
  global node and element ranges for this rank, fills `p_local`
  and `t_local` with global vertex IDs, hands them to
  `solver.set_mesh_distributed`, registers four boundary
  predicates that tag the four edges of the unit square, and
  runs the solve. After convergence each rank reads its locally
  owned part of `udg`, and an `MPI_Allreduce` collects the global
  maximum absolute value for printing on rank 0.
- `CMakeLists.txt` builds the binary as
  `tutorial_06_handwritten_facade_mpi`.

## Build and run

```bash
cd $EXASIM
cmake --build build_mpi --target tutorial_06_handwritten_facade_mpi
cd $EXASIM/tutorial/06-handwritten-facade-mpi
EXASIM_DIR=$EXASIM mpirun -np 2 \
    $EXASIM/build_mpi/tutorial_06_handwritten_facade_mpi
```

The expected output is a global `max|udg|` of approximately
`3.14158`. The MPI+GPU variant builds with `build_mpi_gpu` and the
binary suffix `_mpi_gpu`.
