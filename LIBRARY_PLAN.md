# Plan: Expose Exasim as a Reusable C++ Library

## Goal

Let a user write their own `main.cpp` and their own `MyModel` struct in their own
project, then link Exasim's FEM machinery (mesh, master, discretization, solver,
preconditioner) as a normal CMake library (`find_package(Exasim)`).

The user's model is **plain C++** — pointwise math + hand-written pointwise
Jacobians, exposed as `KOKKOS_INLINE_FUNCTION static` methods on a struct.
**No autodiff, no SymEngine, no DSL.** Templates appear only as the conduit
between the user's struct and Exasim's pre-shipped kernel headers; users do not
write any template code themselves.

`text2code` keeps working unchanged — it becomes one of two ways to populate
the model interface (the other being hand-written code).

## Non-goals

- Removing `text2code`, SymEngine, or the codegen flow.
- Templating the entire FEM library on `Model` (`CSolution<M>`, etc.). That is
  a possible future step; this plan stops at the `libpdemodel.hpp` ABI as the
  link boundary.
- Changing any numerics. Every step is a build / packaging refactor.

## Current state (anchor points)

- `backend/Main/main.cpp` `#include`s implementation `.cpp` files from
  `Discretization/`, `Solver/`, `Preconditioning/`, `Visualization/`,
  `PointLocator/`, `Solution/`, `Preprocessing/`. The runtime is one TU.
- The model is loaded as a shared library (`libpdemodelserial.{so,dylib}`,
  `libpdemodelcuda.*`, `libpdemodelhip.*`) produced by `text2code`. Its ABI
  is `backend/Model/libpdemodel.hpp` — ~30 free functions (`KokkosFlux`,
  `HdgFbou`, `KokkosSource`, …).
- `install/CMakeLists.txt` builds an executable matrix
  `{cpu,gpu} × {serial,mpi} × {builtin, t2c, plain}` that all link the same
  `main.cpp` with different compile definitions.
- Kokkos / METIS / SymEngine are vendored under `kokkos/`, `metis/`,
  `text2code/symengine/` and built separately.

---

## Phase 0 — Inventory & freeze (½ day)

### 0.1 Audit the runtime TU

List every file `backend/Main/main.cpp` pulls in via `#include "*.cpp"`. Confirm
each can be split into a `.hpp` + `.cpp` pair without touching numerics. Files
that are header-only (templates, inline-heavy) stay header-only.

### 0.2 Snapshot a regression baseline

Pick three examples — `apps/poisson/poisson2d`, `apps/poisson/poisson3d`,
`apps/navierstokes/naca0012steady` — and record their `dataout/` outputs after
running with the current `cput2cEXASIM` and `cpumpit2cEXASIM`. These are the
diff targets every later phase must reproduce bit-for-bit (or up to a tight
tolerance).

### 0.3 Fix the Kokkos backend layout

Decide that we ship **one static library per Kokkos backend**:
`libexasim_core_serial.a`, `libexasim_core_cuda.a`, `libexasim_core_hip.a`.
A CMake `Exasim::core` interface target chooses the right one based on which
Kokkos build was selected.

---

## Phase 1 — Split runtime into a real static library (3–5 days)

### 1.1 Move headers into a public include tree

Create `include/exasim/` and move *only the public-facing headers* there:

```
include/exasim/
  common.h          (was backend/Common/common.h)
  cpuimpl.h
  kokkosimpl.h
  pblas.h
  mesh.hpp          (was Preprocessing public surface)
  master.hpp
  discretization.hpp
  preconditioner.hpp
  solver.hpp
  solution.hpp
  visualization.hpp
  pointlocator.hpp
  libpdemodel.hpp   (was backend/Model/libpdemodel.hpp — the model ABI)
  run.hpp           (NEW — top-level driver, see 1.4)
```

Implementation-only headers (helpers, inline kernels currently in `Common/`)
stay under `backend/<dir>/internal/` and are **not** installed.

### 1.2 Convert `.cpp`-included-from-main into compiled TUs

For each subdirectory under `backend/` that's currently consumed by
`main.cpp`'s `#include` chain:

1. Create a `.cpp` that compiles standalone against its public `.hpp`.
2. Replace `#include "../X/y.cpp"` in `main.cpp` with `#include <exasim/y.hpp>`.
3. Add the new `.cpp` to the `exasim_core` library sources.

Subdirs to process, in dependency order: `Common` → `Preprocessing` →
`Discretization` → `Preconditioning` → `Solver` → `Postprocessing` →
`PointLocator` → `Visualization` → `Solution`.

Validation: after each subdir flip, the existing `cput2cEXASIM` /
`cpumpit2cEXASIM` builds must still produce binaries identical (modulo
timestamps) to the baseline outputs in 0.2.

### 1.3 Build `exasim_core` per backend

In `install/CMakeLists.txt`, add:

```cmake
add_library(exasim_core_serial STATIC ${CORE_SOURCES})
target_include_directories(exasim_core_serial
    PUBLIC  $<BUILD_INTERFACE:${EXASIM_DIR}/include>
            $<INSTALL_INTERFACE:include>)
target_link_libraries(exasim_core_serial PUBLIC Kokkos::kokkos ${lapackblas_libraries})
target_compile_features(exasim_core_serial PUBLIC cxx_std_17)

# CUDA / HIP variants gated on EXASIM_CUDA / EXASIM_HIP, with the appropriate
# Kokkos_DIR and target_compile_definitions(... PUBLIC _CUDA / _HIP).
```

The existing executable matrix at `install/CMakeLists.txt:115-365` collapses:
each executable becomes
`add_executable(NAME ${MAIN_SOURCES})` +
`target_link_libraries(NAME PRIVATE exasim_core_<backend> <model_lib>)`.
The `_MPI` / `_CUDA` / `_HIP` / `_TEXT2CODE` defines move to the library
target as `PUBLIC` interface definitions where they affect public headers,
`PRIVATE` otherwise.

### 1.4 Refactor `main.cpp` into `exasim::run`

The body of `backend/Main/main.cpp` (MPI init, Kokkos init, `parseInputFile`,
`CSolution` orchestration, finalize) moves into
`include/exasim/run.hpp` + `src/run.cpp` (compiled into `exasim_core`):

```cpp
// include/exasim/run.hpp
namespace exasim {
int run(int argc, char** argv);
}
```

`backend/Main/main.cpp` becomes:

```cpp
#include <exasim/run.hpp>
int main(int argc, char** argv) { return exasim::run(argc, argv); }
```

Phase-1 exit criterion: every existing executable variant builds, every
example reproduces baseline output. **No user-visible API yet** — this is
pure mechanical refactoring.

---

## Phase 2 — Install + `find_package(Exasim)` (1–2 days)

### 2.1 CMake install rules

```cmake
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS exasim_core_serial
        EXPORT  ExasimTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/exasim
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ExasimTargets
        NAMESPACE Exasim::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Exasim)

configure_package_config_file(
    cmake/ExasimConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Exasim)
```

`ExasimConfig.cmake.in` calls `find_dependency(Kokkos)`, `find_dependency(MPI)`
(if built with MPI), `find_dependency(BLAS)`, `find_dependency(LAPACK)`,
then includes `ExasimTargets.cmake`.

### 2.2 Smoke-test from outside the source tree

In a scratch directory:

```cmake
cmake_minimum_required(VERSION 3.16)
project(exasim_smoke CXX)
find_package(Exasim REQUIRED)
add_executable(smoke smoke.cpp)
target_link_libraries(smoke PRIVATE Exasim::core)
```

`smoke.cpp` just calls `exasim::run(argc, argv)`. Linking against a
text2code-built `libpdemodelserial` provides the model symbols at link time.
This proves the install layout works *before* we ship templated kernels.

---

## Phase 3 — Templated kernel headers (3–4 days)

### 3.1 The model contract

Document in `include/exasim/model.hpp`:

```cpp
// A user model must provide these as static members:
//
//   static constexpr int nd;          // spatial dimension
//   static constexpr int ncu;         // # primary unknowns
//   static constexpr int ncw;         // # auxiliary scalar fields (often 0)
//   static constexpr int nparam;      // # physics parameters
//
//   // Pointwise value functions: doubled as KOKKOS_INLINE_FUNCTION static
//   //   void flux   (double f [/*ncu*nd*/], const double x[/*nd*/],
//   //                const double uq[/*ncu*(1+nd)*/], const double w[/*ncw*/],
//   //                const double mu[/*nparam*/], double t);
//   //   void source (double s [/*ncu*/], ...);
//   //   void ubou   (double ub[/*ncu*/], int ib, ...);
//   //   void fbou   (double fb[/*ncu*/], int ib, ...);
//   //   void initu  (double ui[/*ncu*/], const double x[/*nd*/], ...);
//   //   ... (full list in this header)
//
//   // Hand-written Jacobians (only those your discretization uses):
//   //   void flux_jac_uq (double f_uq [/*ncu*nd*Nq*/], ...);   // Nq = ncu*(1+nd)
//   //   void flux_jac_w  (double f_w  [/*ncu*nd*ncw*/], ...);
//   //   void source_jac_uq (...);
//   //   void fbou_jac_uq (...);
//   //   void fbou_jac_uhg (...);
//   //   ...
```

Provide `exasim::ModelDefaults` as a CRTP base supplying zero-fills for every
optional method, so users only override what their PDE needs.

### 3.2 Templated kernels — one header per family

Create `include/exasim/kernels/`:

```
kernels/flux.hpp         // kokkos_flux_kernel<M>, hdg_flux_kernel<M>
kernels/source.hpp
kernels/sourcew.hpp
kernels/boundary.hpp     // fbou, ubou, fhat, uhat, stab kernels
kernels/init.hpp         // initu, initq, initwdg, initudg, initodg
kernels/tdfunc.hpp
kernels/eos.hpp          // eos, eos_du, eos_dw, avfield
kernels/output.hpp       // output, monitor
kernels/visualization.hpp// vis_scalars, vis_vectors, vis_tensors
kernels/qoi.hpp          // qoi_volume, qoi_boundary
```

Each header contains 1–2 templates: a Kokkos value kernel and (where the FEM
library needs Jacobians) an HDG variant that calls the user's
`*_jac_*` static methods. Every template does:

1. SoA gather: `xdg[k*ng+i] → x[k]`, `udg[k*ng+i] → uq[k]`, etc.
2. Call `M::foo(...)` with stack-allocated outputs.
3. SoA scatter: `out[k] → f[k*ng+i]`.

Reference shape (from `backend/Model/KokkosFlux1.cpp`):

```cpp
namespace exasim {
template <class M>
void kokkos_flux_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        const dstype* uinf, const dstype* param,
                        dstype t, int ng) {
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu*(1+nd);
    Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw];
        for (int k=0;k<nd;++k)  x[k]  = xdg[k*ng+i];
        for (int k=0;k<Nq;++k)  uq[k] = udg[k*ng+i];
        for (int k=0;k<ncw;++k) w[k]  = wdg[k*ng+i];
        double fi[ncu*nd];
        M::flux(fi, x, uq, w, param, t);
        for (int k=0;k<ncu*nd;++k) f[k*ng+i] = fi[k];
    });
}
} // namespace exasim
```

### 3.3 Validate kernels against the codegen path

For Poisson 2D: write `Poisson2D.hpp` by hand (struct with the static methods),
write a `glue.cpp` that defines each `libpdemodel.hpp` symbol as a one-line
call to the corresponding `exasim::*_kernel<Poisson2D>(...)`. Build a binary
that links `Exasim::core` + this `glue.cpp` instead of the text2code `.so`.

Run against the same `pdeapp.txt` as the baseline. Outputs must match.

---

## Phase 4 — User-facing example (1 day)

Create `apps/library_example/poisson2d/` to show the external-project
workflow. Note: this lives inside the Exasim repo for CI, but is structured
as if it were an external project (uses `find_package(Exasim)`, not relative
paths into `backend/`).

```
apps/library_example/poisson2d/
├── CMakeLists.txt
├── main.cpp
├── poisson2d.hpp
├── glue.cpp
└── pdeapp.txt        # same runtime config as apps/poisson/poisson2d
```

Contents are shown in the **End-state usage** section below.

CI: build this example, run it, diff against the baseline.

---

## Phase 5 — Documentation (1 day)

1. `doc/library_api.md`: the model contract, the kernel templates, the ABI
   between user `glue.cpp` and `exasim_core`, the SoA layout convention
   (`c*ng + i`), the Jacobian indexing convention (row-major `i*Nq + j`).
2. Update `README.md` with a third install option ("Use Exasim as a library")
   alongside the existing text2code and Julia/Python/Matlab paths.
3. List which `libpdemodel.hpp` symbols each scheme actually calls (DG vs HDG
   vs HDG with `w`-fields), so users can tell which `M::*` methods they must
   implement vs which they can leave as `ModelDefaults` no-ops.

---

## Phase 6 — Optional polish (later)

- `<exasim/bind_model.inl>`: a single include that, given
  `#define EXASIM_MODEL Poisson2D` first, stamps out all 30 `extern "C"`
  wrappers. Sugar for users who don't want to write `glue.cpp` by hand.
- Template `CSolution`/`CDiscretization` on `M`. Removes the `extern "C"`
  hop entirely — `exasim::run<Poisson2D>(argc, argv)`. Real refactor of
  ~10 classes; only worth doing if profiling shows the ABI hop matters
  (it almost certainly doesn't — these are coarse-grained calls).
- Ship a CMake function `exasim_add_model(target MODEL Poisson2D HEADER ...)`
  that wires up the glue file generation automatically.

---

## End-state usage — what a user does

After Phase 4, a user with Exasim installed (e.g. `cmake --install build`)
writes the following four files in their own project:

### `poisson2d.hpp` — the math

```cpp
#pragma once
#include <exasim/model.hpp>
#include <Kokkos_Core.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;
    static constexpr int ncu    = 1;
    static constexpr int ncw    = 0;
    static constexpr int nparam = 1;
    static constexpr int Nq     = ncu*(1+nd);    // 3

    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double[], const double uq[],
              const double[], const double mu[], double) {
        f[0] = mu[0] * uq[1];
        f[1] = mu[0] * uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double[], const double[],
                     const double[], const double mu[], double) {
        for (int k = 0; k < ncu*nd*Nq; ++k) f_uq[k] = 0.0;
        f_uq[0*Nq + 1] = mu[0];
        f_uq[1*Nq + 2] = mu[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double[], const double[],
                const double[], double) {
        s[0] = 2.0 * M_PI * M_PI
             * Kokkos::sin(M_PI * x[0])
             * Kokkos::sin(M_PI * x[1]);
    }

    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/, const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], double) {
        ub[0] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double[], const double[], const double[]) {
        ui[0] = 0.0;
    }
    // fbou, fhat, vis_*, qoi_* inherited as no-ops from ModelDefaults
};
```

### `glue.cpp` — bind the model into the ABI (one line per kernel)

```cpp
#include <exasim/libpdemodel.hpp>
#include <exasim/kernels/flux.hpp>
#include <exasim/kernels/source.hpp>
#include <exasim/kernels/boundary.hpp>
#include <exasim/kernels/init.hpp>
#include "poisson2d.hpp"

using M = Poisson2D;

void KokkosFlux(dstype* f, const dstype* x, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uinf, const dstype* param,
                dstype t, int, int ng, int, int, int, int, int, int) {
    exasim::kokkos_flux_kernel<M>(f, x, udg, odg, wdg, uinf, param, t, ng);
}

void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg,
             const dstype* x, const dstype* udg, const dstype* odg,
             const dstype* wdg, const dstype* uinf, const dstype* param,
             dstype t, int, int ng, int, int, int, int, int, int) {
    exasim::hdg_flux_kernel<M>(f, f_udg, f_wdg, x, udg, odg, wdg, uinf, param, t, ng);
}

void KokkosSource(dstype* f, const dstype* x, const dstype* udg, const dstype* odg,
                  const dstype* wdg, const dstype* uinf, const dstype* param,
                  dstype t, int, int ng, int, int, int, int, int, int) {
    exasim::kokkos_source_kernel<M>(f, x, udg, odg, wdg, uinf, param, t, ng);
}

void KokkosUbou(dstype* f, const dstype* x, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* uhg, const dstype* nlg,
                const dstype* tau, const dstype* uinf, const dstype* param,
                dstype t, int, int ib, int ng, int, int, int, int, int, int) {
    exasim::kokkos_ubou_kernel<M>(f, x, udg, odg, wdg, uhg, nlg, tau, uinf, param, t, ib, ng);
}

void KokkosInitu(dstype* f, const dstype* x, const dstype* uinf,
                 const dstype* param, int, int ng, int, int, int, int) {
    exasim::kokkos_initu_kernel<M>(f, x, uinf, param, ng);
}

// ... one one-liner per remaining libpdemodel.hpp symbol.
// Optional kernels (vis_*, qoi_*, fhat, etc.) get no-op wrappers when the
// model uses ModelDefaults — or omit the wrappers and provide weak symbols
// in exasim_core (decision in 3.1).
```

### `main.cpp` — five lines

```cpp
#include <exasim/run.hpp>
int main(int argc, char** argv) {
    return exasim::run(argc, argv);
}
```

The user can also write a richer main that constructs an `exasim::Solution`
directly and drives time stepping themselves; `exasim::run` is just the
default convenience.

### `CMakeLists.txt` — the user's build

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_poisson CXX)

find_package(Exasim REQUIRED)

add_executable(poisson_solver main.cpp glue.cpp)
target_link_libraries(poisson_solver PRIVATE Exasim::core)
```

`Exasim::core` transitively pulls in Kokkos, BLAS/LAPACK, MPI (if Exasim was
built with MPI), and exposes `<exasim/...>` headers.

### Running

```
cd my_poisson/build
cmake .. && cmake --build .
mpirun -np 4 ./poisson_solver pdeapp.txt
```

`pdeapp.txt` is the same runtime config file used by the codegen path —
porder, GMRES tolerances, BC tags, mesh path. It has no model math in it;
that lives in `poisson2d.hpp`.

---

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| `#include "*.cpp"` pattern hides ODR / duplicate-symbol issues that surface only when files become independent TUs | Phase 1 splits one subdir at a time and rebuilds; catch issues incrementally rather than at the end |
| Existing `_MPI`/`_CUDA`/`_HIP`/`_TEXT2CODE` defines need to reach the right TUs | Make them `PUBLIC` interface definitions on the library target so anything linking `Exasim::core_*` inherits them. Audit every `#ifdef` site once |
| Per-Kokkos-backend libraries triple build time | One library per backend is unavoidable (Kokkos device code differs). Use `ccache`; treat the three backends as separate CMake build trees if needed |
| Multiple PDE models in one run (`nummodels > 1` in `main.cpp`) | Out of scope for v1. Document as "single model per binary"; multi-model users keep using the `text2code` path until Phase 6 templating |
| Optional kernels (`vis_*`, `qoi_*`, `fhat`) — does user need to define every symbol? | Provide weak no-op definitions in `exasim_core` (or a separate `exasim_default_kernels` object library users can opt out of). User's `glue.cpp` only defines the symbols they actually need |
| Numerical regression from the refactor | Phase 0.2 baselines + CI diff on `apps/library_example/` keep us honest |

---

## Effort estimate

| Phase | Effort | Touches |
|---|---|---|
| 0. Inventory & baseline | 0.5 d | none |
| 1. Static library | 3–5 d | `backend/`, `install/CMakeLists.txt` |
| 2. Install + find_package | 1–2 d | `install/CMakeLists.txt`, new `cmake/ExasimConfig.cmake.in` |
| 3. Templated kernel headers | 3–4 d | new `include/exasim/kernels/` |
| 4. User-facing example | 1 d | new `apps/library_example/` |
| 5. Documentation | 1 d | `doc/`, `README.md` |
| **Total to "Exasim as a library"** | **~10 working days** | |
| 6. Optional polish | 3–5 d | as needed |

Phases 1, 2, 4 are independently shippable. Phase 3 unlocks the
"no-codegen" model-authoring workflow and is the headline value — but
phases 1+2 alone already let users link `Exasim::core` against a
text2code-generated `.so`, which is itself a useful intermediate state.
