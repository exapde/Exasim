# Getting started — writing a new PDE in Exasim

This walks through building a Poisson 2D solver from scratch:

- the math you write in `my_model.hpp`
- the runtime configuration in `pdeapp.txt`
- the 3-line `main.cpp`
- how to build it in-tree (use one of the
  `apps/library_example/<name>/` slots) and out-of-tree (with
  `find_package(Exasim)`).

For the full `Model` contract reference (every method, every signature,
when it's called) see `doc/header_only_api.md`.

For the codegen authoring path (write `pdemodel.txt` in the SymEngine
DSL, run `text2code`, get `my_model.hpp` for free) see
`apps/library_example/README.md`. This document covers the
hand-written path.

---

## The PDE

Poisson with a manufactured solution:

    -∇·(μ ∇u) = 2π² sin(πx) sin(πy)   on (0,1)²
              u = 0                    on the boundary

Exact solution: `u(x,y) = sin(πx) sin(πy)`.

In the HDG / mixed framework Exasim uses, primary unknown `u` is paired
with its gradient `q = ∇u`. The discrete state `uq` carries both:
`uq = [u, q_1, q_2, …]` of total length `Nq = ncu * (1 + nd)`. For this
problem `ncu = 1, nd = 2`, so `Nq = 3` and `uq[0] = u, uq[1] = du/dx,
uq[2] = du/dy`.

The HDG flux is `f = μ ∇u`, i.e. `f[0] = μ uq[1]`, `f[1] = μ uq[2]`.

## Step 1: `my_model.hpp`

Inherit `exasim::ModelDefaults<Self>` (CRTP) and override the methods
your PDE needs. Methods you don't define get zero-fill defaults from
the base class.

```cpp
#pragma once
#include <exasim/model.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    // ---- Compile-time configuration ---------------------------
    static constexpr int nd     = 2;   // spatial dimension
    static constexpr int ncu    = 1;   // # primary unknowns (u)
    static constexpr int ncw    = 0;   // # auxiliary scalar fields (none)
    static constexpr int nco    = 0;   // # auxiliary "other DG" fields (none)
    static constexpr int nparam = 1;   // # physics parameters (μ)

    static constexpr auto disc  = exasim::Discretization::HDG;
    static constexpr int Nq     = ncu * (1 + nd);   // 3

    // ---- Volume terms -----------------------------------------

    // f = μ ∇u
    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double /*x*/[], const double uq[],
              const double /*v*/[], const double /*w*/[],
              const double mu[],
              const double /*uinf*/[], double /*t*/) {
        f[0] = mu[0] * uq[1];
        f[1] = mu[0] * uq[2];
    }

    // s = 2π² sin(πx) sin(πy)
    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double /*uq*/[],
                const double /*v*/[], const double /*w*/[],
                const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        constexpr double pi = 3.141592653589793;
        s[0] = 2.0 * pi * pi
             * Kokkos::sin(pi * x[0]) * Kokkos::sin(pi * x[1]);
    }

    // ---- HDG flux Jacobian ∂f/∂uq ------------------------------
    //
    // Layout (the kernels in <exasim/kernels/flux.hpp> assume this):
    //   f_uq[j * (ncu*nd) + i] = ∂f[i] / ∂uq[j]
    //
    // For Poisson 2D (ncu*nd = 2, Nq = 3, total = 6 entries):
    //   ∂f[0]/∂uq[1] = μ  →  index 1*2+0 = 2
    //   ∂f[1]/∂uq[2] = μ  →  index 2*2+1 = 5
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[],
                     const double /*x*/[], const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[],
                     const double mu[],
                     const double /*uinf*/[], double /*t*/) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }

    // ---- HDG boundary residual + its Jacobian wrt uhat --------
    //
    // The HDG trace ûh is pinned by the trace-side residual.
    // For a Dirichlet u = 0 boundary:  fb = -τ ûh.
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[], const double /*uq*/[],
                  const double /*v*/[], const double /*w*/[],
                  const double uh[],
                  const double /*n*/[], const double tau[],
                  const double /*mu*/[], const double /*uinf*/[],
                  double /*t*/) {
        fb[0] = -tau[0] * uh[0];
    }

    // ∂fbou_hdg/∂uh = -τ (1×1 block).
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[], const double /*uq*/[],
                         const double /*v*/[], const double /*w*/[],
                         const double /*uh*/[],
                         const double /*n*/[], const double tau[],
                         const double /*mu*/[], const double /*uinf*/[],
                         double /*t*/) {
        fb_uh[0] = -tau[0];
    }

    // ---- Initial condition (zero) ------------------------------
    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double /*x*/[],
               const double /*uinf*/[], const double /*mu*/[]) {
        ui[0] = 0.0;
    }

    // ---- LDG-style boundary flux for residual evaluation -------
    // (Exact form from the pdemodel.txt: f·n + τ(uq[0] - uh[0]).)
    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/,
              const double x[], const double uq[],
              const double v[], const double w[], const double uh[],
              const double n[], const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0]*n[0] + f_local[1]*n[1]
              + tau[0]*(uq[0] - uh[0]);
    }

    // ---- ubou: Dirichlet value (here: zero) --------------------
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[], const double /*uq*/[],
              const double /*v*/[], const double /*w*/[],
              const double /*uh*/[],
              const double /*n*/[], const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[],
              double /*t*/) {
        ub[0] = 0.0;
    }
};
```

What you don't write: `tdfunc`, `sourcew`, `eos`, `vis_*`, `qoi_*`,
`fhat`, `uhat`, `stab`, the LDG-path Jacobians, the `*_jac_w` blocks,
the `init{q,udg,wdg,odg}` initializers — `ModelDefaults<Poisson2D>` zero-
fills them. The runtime never calls into a method whose default makes
the math wrong; defaults are deliberately "absent" semantics.

## Step 2: `pdeapp.txt`

Runtime configuration. Sets discretization, mesh file, boundary
conditions, solver tolerances, output paths. See
`apps/poisson/poisson2d/pdeapp.txt` for the canonical example. The
key fields:

```
model           = "ModelD";
modelfile       = "pdemodel.txt";        # not used by hand-written
                                         # path; only by text2code
meshfile        = "grid.bin";
discretization  = "hdg";
platform        = "cpu";

mpiprocs        = 1;
gendatain       = 1;                     # 0 = use existing datain/, 1 = regenerate
porder          = 3;                     # polynomial order
pgauss          = 6;                     # quadrature order

ncu             = 1;                     # # primary unknowns (must match
                                         # MyModel::ncu)
ncw             = 0;
torder          = 1;                     # 1 = steady, 2/3/… = DIRK
dt              = [0];                   # ignored when torder=1
NewtonTol       = 1e-06;
GMREStol        = 1e-08;
physicsparam    = [1];                   # → mu[] in your Model methods
boundaryconditions   = [1, 1, 1, 1];     # boundary tags (your model's
                                         # `int ib` discriminator)
boundaryexpressions  = ["abs(y)<1e-8", "abs(x-1)<1e-8",
                        "abs(y-1)<1e-8", "abs(x)<1e-8"];
```

## Step 3: `main.cpp` (3 lines of code)

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<Poisson2D>(argc, argv);
}
```

`exasim::run<Poisson2D>(argc, argv)` does the MPI / Kokkos init,
parses `pdeapp.txt`, builds the `CSolution<Poisson2D>` object,
runs the Newton/GMRES solve (or DIRK time loop), saves outputs,
and finalizes.

## Step 4: Building (in-tree)

Drop the three files into `apps/library_example/poisson2d_v2/`,
copy `grid.bin` from `apps/poisson/poisson2d/grid.bin`, and add to
`install/CMakeLists.txt`:

```cmake
add_codegen_example(poisson2d_v2)
```

Even though this is hand-written, the `add_codegen_example()` helper
is the right CMake hook — it sets up the same target shape. Then
`cmake --build build --target poisson2d_v2_codegen`.

## Step 5: Building (out-of-tree)

Once Exasim is installed (`cmake --install build --prefix /opt/exasim`),
your project's CMakeLists.txt:

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_solver CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Kokkos REQUIRED PATHS /opt/exasim/external/kokkos NO_DEFAULT_PATH)
find_package(MPI    REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Exasim REQUIRED PATHS /opt/exasim NO_DEFAULT_PATH)

add_executable(my_solver main.cpp)
target_link_libraries(my_solver PRIVATE
    Exasim::headers Kokkos::kokkos MPI::MPI_CXX
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
```

`Exasim::headers` is INTERFACE-only — no compiled library is shipped.
The runtime exists in your TU after the templated kernels inline your
pointwise math at build time.

## Running

```bash
cd /path/to/your_app_dir
my_solver ./pdeapp.txt
```

The binary writes `dataout/outudg_np<rank>.bin`, `dataout/outuhat_np<rank>.bin`,
optionally `dataout/outqoi.txt` (if you defined `qoi_volume` /
`qoi_boundary` in `MyModel`), and `outvis_*.vtu` files.

## Validating

For a manufactured solution, the `qoi_volume` method can return
`(u - uexact)²` — Exasim integrates that over the domain and you read
the squared L² error from `outqoi.txt`. See `apps/library_example/poisson2d/my_model.hpp`'s
`qoi_volume` method for a working example.

For comparing your hand-written model to a text2code-generated one,
or comparing two builds bit-for-bit, see
`apps/library_example/validate_codegen.sh`.

## Common gotchas

**Newton diverges with garbage residual at iter 1:** check that your
Jacobians are layout-correct (column-major, j outer / i inner; see the
inline comment on `flux_jac_uq` above). The runtime can converge the
first iteration on a bad Jacobian and then diverge.

**`ncu_runtime != ncu` assertion at iter 0:** the compile-time `ncu`
in your `Model` doesn't match the runtime's view from `pdeapp.txt`.
Make sure both agree.

**Initial state is all zeros, then NaN:** if your problem requires a
non-trivial initial condition (NS freestream, e.g.), make sure you
override `initu(ui, x, uinf, mu)` to set it. The `ModelDefaults` default
zero-fills.

**Build failure mentioning `cout`:** `<exasim/run.hpp>` opens with
`using namespace std;` at file scope before pulling in the backend
headers (which use unqualified `cout`, `vector`, etc.). If you need a
clean global namespace, write your own MPI/Kokkos init and instantiate
`CSolution<MyModel>` directly instead of going through the façade.
