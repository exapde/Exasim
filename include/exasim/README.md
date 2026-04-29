# `<exasim/...>` — Exasim public C++ API

These headers are the supported surface for using Exasim as a library
from an external CMake project (`find_package(Exasim)` →
`Exasim::core`).

During the library port (see `LIBRARY_PLAN.md` and
`LIBRARY_PORT_INVENTORY.md` at the repo root), most files here are
**thin forwarding shims** that re-include the existing
`backend/<dir>/*.h` headers. As Phase 1.2 progresses, the bodies migrate
into proper `.cpp` translation units in `src/exasim_core/` and the
forwarding shims become the canonical headers.

## Layout

```
include/exasim/
  common.h           — types, macros, struct definitions (forwards backend/Common/common.h)
  cpuimpl.h          — CPU-side numerical helpers (declarations only after Phase 1.2)
  kokkosimpl.h       — Kokkos-side numerical helpers (declarations only after Phase 1.2)
  pblas.h            — BLAS/MPI wrappers
  discretization.hpp — CDiscretization
  preconditioner.hpp — CPreconditioner
  solver.hpp         — CSolver
  solution.hpp       — CSolution
  visualization.hpp  — CVisualization
  pointlocator.hpp   — CPointLocator
  preprocessing.hpp  — CPreprocessing (only if Exasim was built with text2code support)
  libpdemodel.hpp    — model ABI (~30 free functions a user fills in)
  run.hpp            — exasim::run(argc, argv) top-level driver
  model.hpp          — (Phase 3) model contract: static constexpr ints + KOKKOS_INLINE_FUNCTION static methods
  kernels/           — (Phase 3) templated kernels users invoke from their glue.cpp
```

## What a user writes

See `apps/library_example/poisson2d/` (added in Phase 4). Three files:

1. `poisson2d.hpp` — a struct with `static constexpr int nd, ncu, ncw, nparam;`
   and `KOKKOS_INLINE_FUNCTION static` methods (`flux`, `flux_jac_uq`,
   `source`, `ubou`, …). Plain pointwise math; Jacobians written by hand.
2. `glue.cpp` — one one-line `extern` definition per `libpdemodel.hpp`
   symbol, calling the corresponding `exasim::*_kernel<Poisson2D>(...)`.
3. `main.cpp` — three lines: `#include <exasim/run.hpp>` and call `exasim::run`.

No DSL, no codegen, no autodiff.
