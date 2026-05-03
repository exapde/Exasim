# 04 — Handwritten model with the CLI wrapper

The user writes the model header (`my_model.hpp`) by hand as a C++
struct that satisfies the `<exasim/model.hpp>` contract: flux,
source, boundary residual, initial condition, Jacobians, and
optional QoI/visualization hooks. There is no DSL, no autodiff,
and no codegen step — the hand-written struct is consumed directly
by the templated FEM internals.

The driver is the same three-line `main.cpp` as section 02:
`exasim::run<Poisson2D>(argc, argv)` reads `pdeapp.txt`, runs
preprocessing, solves, and writes outputs. The model template
parameter is the user's `Poisson2D` struct rather than a
`text2code`-generated one.

## Files

- `my_model.hpp` — C++ struct `Poisson2D` implementing the
  `<exasim/model.hpp>` contract.
- `pdeapp.txt` — runtime configuration consumed by `exasim::run`.
- `main.cpp` — three-line entry point.
- `CMakeLists.txt` — out-of-tree build for `tutorial_04_handwritten_cli`.
- `grid.bin` — binary mesh for the unit square.
- `pdemodel.txt` — present in this directory only because the
  `text2code` toolchain expects to find one when it generates the
  shared library `libpdemodelserial.{so,dylib}` that the legacy
  ABI links against. The hand-written `Poisson2D` struct is the
  one actually compiled into the binary.

## Build and run

The build invokes `text2code` first to populate
`libpdemodelserial.{so,dylib}` (used as a link-time placeholder),
then compiles `main.cpp` against the hand-written header.

### CPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/04-handwritten-cli \
    $EXASIM/tutorial/04-handwritten-cli/pdeapp.txt
cmake --build build --target tutorial_04_handwritten_cli
cd $EXASIM/tutorial/04-handwritten-cli
mkdir -p datain dataout
$EXASIM/build/tutorial_04_handwritten_cli ./pdeapp.txt
```

### GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/04-handwritten-cli \
    $EXASIM/tutorial/04-handwritten-cli/pdeapp.txt
cmake --build build_gpu --target tutorial_04_handwritten_cli
cd $EXASIM/tutorial/04-handwritten-cli
mkdir -p datain dataout
$EXASIM/build_gpu/tutorial_04_handwritten_cli ./pdeapp.txt
```

### MPI

Edit `pdeapp.txt` so `mpiprocs = N` matches the rank count.

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/04-handwritten-cli \
    $EXASIM/tutorial/04-handwritten-cli/pdeapp.txt
cmake --build build_mpi --target tutorial_04_handwritten_cli
cd $EXASIM/tutorial/04-handwritten-cli
mkdir -p datain dataout
mpirun -np 2 $EXASIM/build_mpi/tutorial_04_handwritten_cli ./pdeapp.txt
```

### MPI+GPU

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/04-handwritten-cli \
    $EXASIM/tutorial/04-handwritten-cli/pdeapp.txt
cmake --build build_mpi_gpu --target tutorial_04_handwritten_cli
cd $EXASIM/tutorial/04-handwritten-cli
mkdir -p datain dataout
mpirun -np 2 $EXASIM/build_mpi_gpu/tutorial_04_handwritten_cli ./pdeapp.txt
```

## File contents

### `my_model.hpp`

```cpp
#pragma once

#include <exasim/model.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    // ---- Compile-time configuration -------------------------------
    static constexpr int nd     = 2;   // spatial dimension
    static constexpr int ncu    = 1;   // # primary unknowns
    static constexpr int ncw    = 0;   // # auxiliary scalar fields
    static constexpr int nco    = 0;   // # other DG fields (`v` / `odg`)
    static constexpr int nparam = 1;   // # physics parameters

    static constexpr auto disc = exasim::Discretization::HDG;

    // Derived: Nq = ncu*(1 + nd) = 3
    static constexpr int Nq = ncu * (1 + nd);

    // ---- Volume terms ---------------------------------------------

    // Flux f = μ ∇u
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

    // Source s = 2π² sin(πx) sin(πy)
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

    // ---- Required initial condition -------------------------------

    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double /*x*/[],
               const double /*uinf*/[], const double /*mu*/[]) {
        ui[0] = 0.0;
    }

    // ---- HDG path -------------------------------------------------

    // Volume-term Jacobian ∂f/∂uq (column-major, uq outer):
    //   ∂f[0]/∂uq[1] = mu  ->  index 1*2+0 = 2
    //   ∂f[1]/∂uq[2] = mu  ->  index 2*2+1 = 5
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

    // FbouHdg = -τ * uhat   (Dirichlet via stabilization)
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[],  const double /*uq*/[],
                  const double /*v*/[],  const double /*w*/[],  const double uh[],
                  const double /*n*/[],  const double tau[],
                  const double /*mu*/[], const double /*uinf*/[],
                  double /*t*/) {
        fb[0] = -tau[0] * uh[0];
    }

    // ∂fbou_hdg/∂uh = -τ
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double tau[],
                         const double /*mu*/[], const double /*uinf*/[],
                         double /*t*/) {
        fb_uh[0] = -tau[0];
    }

    // ---- LDG boundary flux Fbou (used by some assembly paths) -----
    // Fbou = f·n + τ*(uq[0] - uhat[0])
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

    // ---- Ubou: zero Dirichlet -------------------------------------
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[],
              double /*t*/) {
        ub[0] = 0.0;
    }

    // ---- Visualization & QoI --------------------------------------

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
contract. Only the methods relevant to Poisson 2D are overridden
here; the rest (`tdfunc`, `sourcew`, `eos`, `vis_*` extensions,
LDG-path Jacobians, `init{q,udg,wdg,odg}`, etc.) come from the
base class.

The compile-time constants (`nd`, `ncu`, `ncw`, `nco`, `nparam`,
`disc`) drive how the templated FEM internals dispatch; the
runtime relies on them to skip calls whose output would be unused.
The HDG branch (`disc = HDG`) drives the assembly through `flux`,
`source`, `flux_jac_uq`, `fbou_hdg`, and `fbou_hdg_jac_uh`. The
LDG branch (`fbou`) is also implemented because some assembly
paths still touch it.

The numeric constant `1.973920880217872E+1` in `source` is `2π²`
to enough digits to match `text2code`'s SymEngine output exactly,
so this hand-written model can reproduce the existing baseline
numerics bit-for-bit when the runtime context is the same.

### `pdeapp.txt`

```
model = "ModelD";
modelfile = "pdemodel.txt";
meshfile = "grid.bin";
discretization = "hdg";
platform = "cpu";

builtinmodelID = 1;
gendatain = 0;
mpiprocs = 1;
debugmode = 0;
runmode = 0;
modelnumber = 0;

porder = 3;
pgauss = 6;
torder = 1;
nstage = 1;

ncu = 1;
ncw = 0;
neb = 4096;
nfb = 8192;

saveSolBouFreq = 1;
ibs = 1;

NewtonIter = 20;
NewtonTol = 1e-08;
GMRESiter = 50;
GMRESrestart = 25;
GMREStol = 1e-08;
GMRESortho = 1;
ppdegree = 1;
RBdim = 0;
matvecorder = 1;
matvectol = 0.001;
precMatrixType = 0;
preconditioner = 1;

time = 0;
tau = [1];
dt = [1];
physicsparam = [1.0];
externalparam = [0, 0];

boundaryconditions = [1];
boundaryexpressions = ["abs(y)<1e-8 || abs(y-1)<1e-8 || abs(x-1)<1e-8 || abs(x)<1e-8"];
curvedboundaries = [0];
curvedboundaryexprs = [""];
periodicboundaries1 = [];
periodicexprs1 = [];
periodicboundaries2 = [];
periodicexprs2 = [];
interfaceconditions = [];
```

`exasim::run<Poisson2D>(argc, argv)` reads this file from
`argv[1]`. `model = "ModelD"` selects the steady-state HDG class.
`porder = 3` and `pgauss = 6` set discretization order and
quadrature order. `physicsparam = [1.0]` sets `mu[0] = 1`, which
the hand-written `flux` and `flux_jac_uq` read from the `mu`
parameter pack. `boundaryconditions = [1]` and
`boundaryexpressions` tag every face on the four edges of the
unit square with id 1, which the model's `fbou_hdg` and `fbou`
treat as zero-Dirichlet. For multi-rank runs `mpiprocs = N` must
match the rank count.

### `main.cpp`

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<Poisson2D>(argc, argv);
}
```

`exasim::run<M>(argc, argv)` does the MPI / Kokkos init, parses
`argv[1]` as a `pdeapp.txt`, builds the `CSolution<M>` object,
runs the Newton/GMRES solve (or DIRK time loop), saves outputs,
and finalizes. The template parameter `Poisson2D` selects the
hand-written struct from `my_model.hpp`.

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)

set(_target tutorial_04_handwritten_cli)

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
`_TEXT2CODE` everywhere. It also picks `libpdemodelserial`,
`libpdemodelcuda`, or `libpdemodelhip` depending on the variant.

The hand-written `Poisson2D` struct is what's actually compiled
into the binary; the link-time `libpdemodel*` library is present
only because the legacy ABI plumbing expects it.

### `grid.bin`

A binary blob in the format `[size(p), size(t), p(:), t(:)]`. See
[section gridbin](../00-overview/README.md#gridbin) of the
overview for the format and origin.
