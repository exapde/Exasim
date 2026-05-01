# `<exasim/...>` — Exasim public C++ API

These headers are the supported surface for using Exasim as a header-only
template library from an external CMake project. After

    find_package(Exasim REQUIRED)
    target_link_libraries(my_solver PRIVATE Exasim::headers)

a downstream consumer pulls in everything below. There is no compiled
Exasim library to link — `Exasim::headers` is INTERFACE-only.

For a guided walkthrough, see `doc/header_only_api.md` (API reference)
and `doc/getting_started.md` (writing a new PDE from scratch).

## Layout

```
include/exasim/
  model.hpp          — Model contract: static constexpr ints +
                       `KOKKOS_INLINE_FUNCTION static` pointwise methods.
                       `ModelDefaults<Self>` CRTP base supplies zero-fill
                       defaults for optional methods.
  run.hpp            — `exasim::run<M>(argc, argv)` façade — MPI/Kokkos
                       init, argv parsing, DIRK time stepping, multi-domain
                       coupling, file I/O, finalize. User main.cpp's
                       collapse to a 3-line wrapper.
  drivers.hpp        — templated `*Driver<M>` wrappers that bridge the
                       FEM-internal call sites and the per-kernel templates.
  kernels/           — templated kernels that gather SoA inputs into a
                       pointwise buffer, call M::method(...) at the
                       quadrature point, and scatter the result back.
                       One header per kernel family (flux, source, boundary,
                       interface, init, tdfunc, eos, sourcew, qoi, output,
                       visualization).
  detail/abi_adapter.hpp        — sentinel type that makes the legacy
                                  libpdemodel.hpp ABI co-exist with the
                                  templated path.
  detail/driver_dispatch.hpp    — `EXASIM_DRIVER_CALL` macro. `if constexpr`
                                  on `M`: AbiAdapter routes to global
                                  non-templated `::Name(...)`; user M routes
                                  to `exasim::Name<M>(...)`.
  discretization.hpp — `CDiscretization<M>` (FEM internals)
  preconditioner.hpp — `CPreconditioner<M>`
  solver.hpp         — `CSolver<M>`
  solution.hpp       — `CSolution<M>` — top-level container; `exasim::run<M>`
                       instantiates this.
  visualization.hpp  — `CVisualization`
  pointlocator.hpp   — `CPointLocator`
  preprocessing.hpp  — `CPreprocessing` (text2code path)
  common.h           — types, macros, struct definitions
  cpuimpl.h          — CPU-side numerical helpers
  kokkosimpl.h       — Kokkos-side numerical helpers
  pblas.h            — BLAS / MPI wrappers
  libpdemodel.hpp    — legacy ~30-symbol ABI, kept for cput2cEXASIM /
                       cpumpit2cEXASIM and for the `EXASIM_DRIVER_CALL`
                       constexpr-else branch.
```

## What a user writes

The minimum is two files: `my_model.hpp` (the math) and `main.cpp`
(the entry point). See `apps/library_example/poisson2d/` for a hand-
written example, or any `apps/library_example/<name>_codegen/` for a
text2code-generated one. main.cpp in either case:

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<MyModel>(argc, argv);
}
```

`MyModel` is a struct inheriting `exasim::ModelDefaults<MyModel>` with
compile-time constants (`nd, ncu, ncw, nco, nparam, disc`) and the
required pointwise methods (`flux`, `source`, `fbou`, `fbou_hdg`,
`ubou`, `initu`, plus HDG Jacobians for HDG models). See
`doc/header_only_api.md` for the full Model contract.
