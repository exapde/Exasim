# 05 — Handwritten model with the embedded solver (single rank)

The user writes the model header (`my_model.hpp`) by hand, then
constructs an `ExasimSolver<Poisson2D>` from C++ and configures
everything in code: the mesh is built as flat C arrays inside
`main.cpp` and handed to `solver.set_mesh(...)`, and the boundary
conditions are registered with predicate lambdas. There is no
`pdeapp.txt`, no `grid.bin`, and no `text2code` step. The
converged solution is read back through `solver.udg()` for
post-processing inside the same program.

This is the embedded path for users who own the mesh data
structures: optimization loops, reduced-order modeling,
PDE-constrained inverse problems, and multi-physics couplings.

This section is **single rank only**: `set_mesh` expects the full
mesh on every rank, no partitioning happens, and `main.cpp` does
not initialize MPI. The CMake build skips the target when
`EXASIM_MPI=ON`. For a true distributed-memory variant where each
rank owns a slice of the mesh and ParMETIS repartitions inside
`solve()`, see [section 06](../06-handwritten-distributed/README.md).

## Files

- `my_model.hpp` — C++ struct `Poisson2D` implementing the
  `<exasim/model.hpp>` contract.
- `main.cpp` — embedded driver. Builds a 16×16 Cartesian quad
  mesh on `[0,1]²`, configures the solver, calls `solve()`,
  prints `max|udg|`. No MPI calls.
- `CMakeLists.txt` — out-of-tree build for
  `tutorial_05_handwritten_embedded`. Only added to the build when
  `EXASIM_MPI=OFF`.

The runtime expects to find the master-element node tables under
`$EXASIM/backend/Preprocessing/`, so the binary needs the
environment variable `EXASIM_DIR` to be set when it runs.

## Build and run

### CPU

```bash
cd $EXASIM
cmake --build build --target tutorial_05_handwritten_embedded
EXASIM_DIR=$EXASIM $EXASIM/build/tutorial_05_handwritten_embedded
```

### GPU

```bash
cd $EXASIM
cmake --build build_gpu --target tutorial_05_handwritten_embedded
EXASIM_DIR=$EXASIM $EXASIM/build_gpu/tutorial_05_handwritten_embedded
```

The expected output is `max|udg| = 3.14158` (the maximum absolute
value of the converged mixed state for the manufactured solution
`u = sin(πx) sin(πy)`).

## File contents

### `my_model.hpp`

```cpp
#pragma once

#include <exasim/model.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;
    static constexpr int ncu    = 1;
    static constexpr int ncw    = 0;
    static constexpr int nco    = 0;
    static constexpr int nparam = 1;

    static constexpr auto disc = exasim::Discretization::HDG;
    static constexpr int Nq = ncu * (1 + nd);

    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double /*x*/[], const double uq[],
              const double /*v*/[], const double /*w*/[], const double mu[],
              const double /*uinf*/[], double /*t*/) {
        const double mu0  = mu[0];
        const double udg2 = uq[1];
        const double udg3 = uq[2];
        f[0] = mu0 * udg2;
        f[1] = mu0 * udg3;
    }

    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double /*uq*/[],
                const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        const double xdg1 = x[0];
        const double xdg2 = x[1];
        s[0] = Kokkos::sin(xdg1 * 3.141592653589793)
             * Kokkos::sin(xdg2 * 3.141592653589793)
             * 1.973920880217872E+1;
    }

    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double /*x*/[],
               const double /*uinf*/[], const double /*mu*/[]) {
        ui[0] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double /*x*/[],
                     const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[],
                     const double mu[], const double /*uinf*/[],
                     double /*t*/) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[],  const double /*uq*/[],
                  const double /*v*/[],  const double /*w*/[],  const double uh[],
                  const double /*n*/[],  const double tau[],
                  const double /*mu*/[], const double /*uinf*/[],
                  double /*t*/) {
        fb[0] = -tau[0] * uh[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double tau[],
                         const double /*mu*/[], const double /*uinf*/[],
                         double /*t*/) {
        fb_uh[0] = -tau[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/,
              const double x[],  const double uq[],
              const double v[],  const double w[],  const double uh[],
              const double n[],  const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }

    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[],
              double /*t*/) {
        ub[0] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(double s[], const double /*x*/[], const double uq[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[0];
        s[1] = uq[1] + uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(double s[], const double /*x*/[], const double uq[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[1];
        s[1] = uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(double s[], const double x[], const double uq[],
                    const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                    const double /*uinf*/[], double /*t*/) {
        const double t1 = 3.141592653589793;
        const double t2 = Kokkos::sin(t1 * x[0]);
        const double t3 = Kokkos::sin(t1 * x[1]);
        const double uexact = t2 * t3;
        s[0] = (uq[0] - uexact) * (uq[0] - uexact);
        s[1] = uq[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(double fb[], int /*ib*/,
                      const double x[],  const double uq[],
                      const double v[],  const double w[],  const double uh[],
                      const double n[],  const double tau[],
                      const double mu[], const double uinf[],
                      double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }
};
```

`Poisson2D` inherits `exasim::ModelDefaults<Poisson2D>` (CRTP),
which supplies zero-fill defaults for every method on the model
contract. The compile-time constants drive how the templated FEM
internals dispatch. The HDG branch consumes `flux`, `source`,
`flux_jac_uq`, `fbou_hdg`, and `fbou_hdg_jac_uh`. The `qoi_volume`
QoI returns `(u − u_exact)²` so the runtime can integrate the
squared L² error.

The numeric constant `1.973920880217872E+1` in `source` is `2π²`
to enough digits to match `text2code`'s SymEngine output exactly.

### `main.cpp`

```cpp
#include <exasim/run.hpp>            // pulls common preamble + kokkos + namespace std
#include <exasim/solver_facade.hpp>
#include <exasim/model.hpp>

#include "my_model.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, char** argv) {
    (void)argc; (void)argv;

    Kokkos::initialize();
    {
        // Build a Cartesian quad mesh on [0,1]^2 with 16x16 elements.
        const int n  = 16;
        const int nv = (n + 1) * (n + 1);
        const int ne = n * n;

        std::vector<double> p(2 * nv);
        std::vector<int>    t(4 * ne);
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                int idx = j * (n + 1) + i;
                p[2 * idx + 0] = double(i) / n;
                p[2 * idx + 1] = double(j) / n;
            }
        }
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int e   = j * n + i;
                int v00 = j * (n + 1) + i;
                t[4 * e + 0] = v00;
                t[4 * e + 1] = v00 + 1;
                t[4 * e + 2] = v00 + 1 + (n + 1);
                t[4 * e + 3] = v00 + (n + 1);
            }
        }

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

        const double* udg   = solver.udg();
        const Int     udg_n = solver.udg_size();

        double maxabs = 0.0;
        for (Int i = 0; i < udg_n; ++i) {
            double v = std::abs(udg[i]);
            if (v > maxabs) maxabs = v;
        }
        std::printf("[tutorial_05] udg: %lld doubles, max|udg| = %.5f\n",
                    static_cast<long long>(udg_n), maxabs);
    }
    Kokkos::finalize();
    return 0;
}
```

`Kokkos::initialize()` brings up the device runtime; everything
between it and `Kokkos::finalize()` runs inside an explicit scope
so all device-allocating objects are destroyed before finalize.

The two nested `for` loops on `(j, i)` build the mesh data:
`p[2*idx + d]` is the `d`th coordinate of node `idx`, and
`t[4*e + k]` is the `k`th vertex of element `e`. Vertices are
laid out in row-major order on a `(n+1) × (n+1)` lattice.

`exasim::ExasimSolver<Poisson2D>` is the templated embedded
solver. `set_mesh(p, t, np, ne, nve)` hands the mesh to the
solver. `add_boundary(tag, predicate)` tags every face whose
midpoint satisfies the predicate with the given boundary id; the
four `add_boundary` calls cover the bottom (`y = 0`), right
(`x = 1`), top (`y = 1`), and left (`x = 0`) edges of the unit
square. `set_polynomial_order(3)` and `set_quadrature_order(6)`
match the `porder` and `pgauss` fields of the file-driven sections.
`set_physics_params({1.0})` sets `mu[0] = 1`. `solve()` (no
arguments — single rank) runs the Newton/GMRES loop in place.

After `solve()`, `solver.udg()` returns a pointer to the
converged mixed state, and `solver.udg_size()` returns its
element count. The driver scans it to find the maximum absolute
value and prints the result.

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)

set(_target tutorial_05_handwritten_embedded)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")
```

`tutorial_configure_target` is a helper defined in
`tutorial/CMakeLists.txt` that adds the right backend defines and
libraries for the active build variant: `_CUDA` and the CUDA
runtime/cuBLAS libraries on `build_gpu` and `build_mpi_gpu`,
`_HIP` on AMD GPUs, `_MPI` on MPI-enabled builds, and
`_TEXT2CODE` everywhere.

The target compiles `main.cpp` against the Exasim headers and
links it against Kokkos, BLAS/LAPACK, ParMETIS (via
`link_metis`), and the `libpdemodelserial.{so,dylib}` placeholder
under `backend/Model/`. The hand-written `Poisson2D` struct is
what's actually compiled into the binary.
