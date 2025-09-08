# Exasim Backend

High-performance HDG/DG backend for PDE simulation with CPU/GPU portability via Kokkos.

- **Common/**: platform abstraction, BLAS/PBLAS stubs.
- **Discretization/**: geometry, kernels, assembly drivers.
- **Model/**: physics kernels (HDG/Kokkos), symbolic DSL → C++.
- **Solver/**: GMRES, PTC, system setup.
- **Preconditioning/**: preconditioners and apply routines.
- **Solution/**: state storage, time stepping, output.

See the docs for subsystem overviews.

