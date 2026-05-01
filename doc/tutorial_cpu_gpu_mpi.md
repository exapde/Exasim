# Tutorial — Solving a PDE with Exasim

End-to-end walkthrough: bootstrap dependencies, write a `Model` in C++,
solve on CPU / GPU / MPI / MPI+GPU.

The primary path here is the **programmatic** one — you write a
`Model` struct, build the mesh as flat arrays in C++, and call
`exasim::ExasimSolver<MyModel>`. The runtime never touches `pdeapp.txt`
or `datain/`/`dataout/`. The codegen path (write `pdemodel.txt`, run
`text2code`, get `my_model.hpp` for free) is documented in §6 as a
convenience for users who'd rather author math symbolically — both
paths produce identical solver behavior.

The four supported backends are configured by which Kokkos build
prefix the binary links against and which `EXASIM_*` flags CMake
sees. They are picked at compile time, not runtime; one source file
plus four CMake targets is the standard setup:

| variant            | what it is                          | binary suffix       |
| ------------------ | ----------------------------------- | ------------------- |
| CPU                | single-rank serial                  | _facade_            |
| MPI                | multi-rank serial CPU               | _facade_mpi_        |
| GPU                | single-rank CUDA or HIP             | _facade_gpu_        |
| MPI + GPU          | multi-rank GPU (one rank per device)| _facade_mpi_gpu_    |

Tested on Ubuntu 22.04 with CUDA 12.x on V100s. HIP follows the same
script with `kokkos/buildhip` and `EXASIM_HIP=ON` substituted.

---

## 1. Install dependencies

Three deps are vendored in-tree (`kokkos/`, `text2code/symengine/`,
`metis/`); two are external (BLAS/LAPACK, optional MPI, optional CUDA
or ROCm). The vendored builds are how Exasim ships its known-good
versions — you don't need a system Kokkos.

### 1a. System packages

```bash
# Ubuntu 22.04
sudo apt install -y build-essential cmake git \
                    libblas-dev liblapack-dev \
                    libmkl-dev    # optional but faster
```

For GPU builds, install the CUDA toolkit (`nvcc` on PATH) or ROCm
(`hipcc` on PATH) ahead of time.

### 1b. OpenMPI with CUDA awareness

If you'll do MPI + GPU, OpenMPI must be built `--with-cuda` so ranks
can pass device pointers directly. The system OpenMPI almost
certainly is not.

```bash
EXASIM=$PWD                   # repo root
MPI_PREFIX=$EXASIM/openmpi    # or wherever you want it

cd $EXASIM
mkdir build_openmpi && cd build_openmpi
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz
tar xf openmpi-4.1.6.tar.gz && cd openmpi-4.1.6

PATH=/usr/local/cuda/bin:$PATH ./configure \
    --prefix=$MPI_PREFIX \
    --with-cuda=/usr/local/cuda \
    --disable-mpi-fortran --disable-oshmem
make -j16 install

export PATH=$MPI_PREFIX/bin:$PATH
ompi_info | grep -i cuda           # must list `cuda` under MPI extensions
```

For CPU-only MPI, drop `--with-cuda=…` and `PATH=/usr/local/cuda/bin:`.

### 1c. Kokkos (serial + CUDA prefixes)

```bash
cd $EXASIM/kokkos

# Serial — CPU + MPI builds use this
mkdir buildserial && cd buildserial
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildserial \
         -DKokkos_ENABLE_SERIAL=ON
make -j16 install
cd ..

# CUDA — GPU + MPI+GPU builds use this. ARCH_VOLTA70 = V100; pick
# yours from kokkos/cmake/kokkos_arch.cmake.
mkdir buildcuda && cd buildcuda
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildcuda \
         -DKokkos_ENABLE_SERIAL=ON \
         -DKokkos_ENABLE_CUDA=ON \
         -DKokkos_ARCH_VOLTA70=ON \
         -DCMAKE_CXX_COMPILER=$EXASIM/kokkos/buildcuda/bin/nvcc_wrapper
make -j16 install
```

For HIP, the second build uses `Kokkos_ENABLE_HIP=ON
Kokkos_ARCH_VEGA90A=ON` (or your arch) and `CMAKE_CXX_COMPILER=hipcc`.

### 1d. METIS / GKlib / ParMETIS

```bash
cd $EXASIM/metis

cd GKlib && make config prefix=$EXASIM/metis/GKlib && make -j8 install && cd ..
cd METIS && make config gklib_path=$EXASIM/metis/GKlib prefix=$EXASIM/metis/METIS \
                  shared=1 r64=1 i64=1 && make -j8 install && cd ..

# Only needed for MPI variants.
PATH=$MPI_PREFIX/bin:$PATH cd ParMETIS && \
    make config gklib_path=$EXASIM/metis/GKlib \
                metis_path=$EXASIM/metis/METIS \
                prefix=$EXASIM/metis/ParMETIS \
                cc=mpicc cxx=mpicxx && \
    make -j8 install
```

### 1e. SymEngine + text2code (optional — only for the codegen path)

You only need these if you're going to use the §6 convenience path
(write `pdemodel.txt`, generate `my_model.hpp` with `text2code`).
The programmatic path doesn't link SymEngine or text2code.

```bash
cd $EXASIM/text2code/symengine
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/text2code/symengine/build \
         -DBUILD_SHARED_LIBS=OFF -DWITH_GMP=OFF -DINTEGER_CLASS=boostmp
make -j16 install

cd $EXASIM/text2code/text2code
cmake -S . -B build -DEXASIM_ROOT=$EXASIM
cmake --build build -j8
```

---

## 2. Configure and build Exasim

Four build directories, one per backend. Each links the appropriate
Kokkos prefix and sets the right `EXASIM_*` flags. Same source code
in every case — the differences live in `CMakeLists.txt` flags.

```bash
cd $EXASIM

# CPU only
cmake -S install -B build_cpu \
    -DEXASIM_CUDA=OFF -DEXASIM_HIP=OFF \
    -DEXASIM_MPI=OFF  -DEXASIM_NOMPI=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildserial/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
cmake --build build_cpu -j16

# GPU (single rank)
cmake -S install -B build_gpu \
    -DEXASIM_CUDA=ON -DEXASIM_NOMPI=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildcuda/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
PATH=/usr/local/cuda/bin:$PATH cmake --build build_gpu -j8

# MPI (multi-rank CPU)
PATH=$MPI_PREFIX/bin:$PATH cmake -S install -B build_mpi \
    -DEXASIM_CUDA=OFF -DEXASIM_HIP=OFF \
    -DEXASIM_MPI=ON   -DEXASIM_NOMPI=OFF \
    -DWITH_PARMETIS=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildserial/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON \
    -DMPI_CXX_COMPILER=$MPI_PREFIX/bin/mpicxx \
    -DMPI_C_COMPILER=$MPI_PREFIX/bin/mpicc
PATH=$MPI_PREFIX/bin:$PATH cmake --build build_mpi -j16

# MPI + GPU (multi-rank GPU)
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake -S install -B build_mpi_gpu \
    -DEXASIM_CUDA=ON -DEXASIM_MPI=ON -DEXASIM_NOMPI=OFF \
    -DWITH_PARMETIS=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildcuda/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build build_mpi_gpu -j8
```

Each build dir produces (among other targets) the matching
`solve_square_facade*` example used in §3.

To install Exasim system-wide:

```bash
cmake --install build_cpu --prefix /opt/exasim
cmake --install build_gpu --prefix /opt/exasim   # adds GPU libs
```

The prefix layout:

```
/opt/exasim/
├── bin/text2code            # convenience codegen tool (§6)
├── include/exasim/          # public template headers
├── include/backend/         # backend headers exasim/* depends on
└── lib/cmake/Exasim/        # ExasimConfig.cmake for find_package
```

---

## 3. Solve a PDE — the programmatic path

The whole flow lives in one `main()`. No `pdeapp.txt`, no
`pdemodel.txt`, no `datain/` writes, no `dataout/` writes. Output
comes back as a host pointer.

We solve Poisson 2D with manufactured solution
`u(x,y) = sin(πx) sin(πy)` on the unit square, source
`s = 2π² sin(πx) sin(πy)`, Dirichlet `u = 0` on the boundary.

### 3.1 Write the Model (HDG)

Inherit `exasim::ModelDefaults<Self>` (CRTP) and override only what
your PDE needs. Methods you don't supply default to zero-fill.

```cpp
#include <exasim/run.hpp>
#include <exasim/model.hpp>
#include <exasim/solver_facade.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;       // spatial dim
    static constexpr int ncu    = 1;       // # primary unknowns (u)
    static constexpr int ncw    = 0;       // no auxiliary scalars
    static constexpr int nco    = 0;
    static constexpr int nparam = 1;       // physics μ
    static constexpr auto disc  = exasim::Discretization::HDG;
    static constexpr int Nq     = ncu * (1 + nd);

    // f = μ ∇u
    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double[], const double uq[],
              const double[], const double[], const double mu[],
              const double[], double) {
        f[0] = mu[0] * uq[1];
        f[1] = mu[0] * uq[2];
    }
    // s = 2π² sin(πx) sin(πy)
    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double[],
                const double[], const double[], const double[],
                const double[], double) {
        constexpr double pi = 3.141592653589793;
        s[0] = 2.0 * pi * pi
             * Kokkos::sin(pi * x[0]) * Kokkos::sin(pi * x[1]);
    }
    // ∂f/∂uq — ncu*nd × Nq Jacobian, column-major
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double[], const double[],
                     const double[], const double[], const double mu[],
                     const double[], double) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }
    // HDG trace Dirichlet: fb = -τ ûh
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int, const double[], const double[],
                  const double[], const double[], const double uh[],
                  const double[], const double tau[],
                  const double[], const double[], double) {
        fb[0] = -tau[0] * uh[0];
    }
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int, const double[],
                         const double[], const double[], const double[],
                         const double[], const double[], const double tau[],
                         const double[], const double[], double) {
        fb_uh[0] = -tau[0];
    }
    // residual side fbou + Dirichlet ubou + zero IC initu — see
    // doc/getting_started.md for full annotations.
    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double[], const double[], const double[]) { ui[0] = 0.0; }
    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int, const double x[], const double uq[],
              const double v[], const double w[], const double uh[],
              const double n[], const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0]*n[0] + f_local[1]*n[1] + tau[0]*(uq[0] - uh[0]);
    }
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int, const double[], const double[],
              const double[], const double[], const double[],
              const double[], const double[], const double[],
              const double[], double) { ub[0] = 0.0; }
};
```

### 3.2 Build the mesh as flat arrays

`p` is `nd × np` doubles (column-major); `t` is `nve × ne` ints
(column-major, **0-based**). For a quad mesh on `[0,1]²` at
resolution `n × n`:

```cpp
const int n  = 16;
const int nv = (n + 1) * (n + 1);
const int ne = n * n;
std::vector<double> p(2 * nv);
std::vector<int>    t(4 * ne);

for (int j = 0; j <= n; ++j)
    for (int i = 0; i <= n; ++i) {
        int idx = j * (n + 1) + i;
        p[2*idx + 0] = double(i) / n;
        p[2*idx + 1] = double(j) / n;
    }
for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
        int e = j * n + i;
        int v00 = j * (n + 1) + i;       // 0-based
        t[4*e + 0] = v00;
        t[4*e + 1] = v00 + 1;
        t[4*e + 2] = v00 + 1 + (n + 1);
        t[4*e + 3] = v00 + (n + 1);
    }
```

### 3.3 Drive the solver

`ExasimSolver<MyModel>` collects the inputs and runs the pipeline
(preprocessing → discretization → Newton/GMRES). Picks CPU vs CUDA
vs HIP at compile time from the binary's defines; same source builds
for all four backends.

```cpp
exasim::ExasimSolver<Poisson2D> solver;
solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4);
solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[1])     < 1e-8; });  // y=0
solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[0] - 1) < 1e-8; });  // x=1
solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[1] - 1) < 1e-8; });  // y=1
solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[0])     < 1e-8; });  // x=0
solver.set_polynomial_order(3);
solver.set_quadrature_order(6);
solver.set_physics_params({1.0});
solver.solve();
```

Predicates are typed `std::function<bool(const double*)>` —
no tinyexpr strings. Curved boundaries take a level-set instead:

```cpp
solver.add_curved_boundary(/*tag=*/1,
    /*classify=*/  [](const double* x){ return std::sqrt(x[0]*x[0]+x[1]*x[1]) < 1.001; },
    /*level_set=*/ [](const double* x){ return x[0]*x[0]+x[1]*x[1] - 1.0; });
```

The projector does one Newton step `x ← x - f·∇f/|∇f|²` toward
`f=0` (∇f from finite differences), so any function whose zero set
is the curved surface works — SDF or otherwise.

### 3.4 Read the solution from memory

```cpp
const double* udg = solver.udg();      // GPU runs: lazy device→host copy
Int           n   = solver.udg_size(); // npe * nc * ne
// udg layout per element: [u₀..u_{npe-1}, du/dx₀..du/dx_{npe-1}, du/dy₀..]
```

`udg()` / `uhat()` / `wdg()` are valid until the `ExasimSolver`
destructs.

### 3.5 The `main()` boilerplate

The binary needs to initialize MPI (when built with `_MPI`) and
Kokkos in the right order. With CUDA the per-rank GPU binding has
to come **before** `Kokkos::initialize`:

```cpp
int main(int argc, char** argv) {
    int mpiprocs = 1, mpirank = 0, shmrank = 0;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
#endif
#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
#ifdef HAVE_HIP
    int nd; hipGetDeviceCount(&nd); hipSetDevice(shmrank % nd);
#endif
    Kokkos::initialize();
    {
        // ...build mesh, instantiate solver, call solver.solve()...
    }
    Kokkos::finalize();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}
```

The full `solve_square_facade.cpp` source (~170 lines) lives in
`apps/library_example/solve_square/`.

---

## 4. Cross-arch — same source, four binaries

The same `solve_square_facade.cpp` is built four times with
different defines/links. CMake handles it; you just pick which
build directory to invoke.

```bash
build_cpu/solve_square_facade                                          # CPU
build_gpu/solve_square_facade_gpu                                      # GPU
mpirun -np 2 build_mpi/solve_square_facade_mpi                         # 2-rank CPU
mpirun -np 2 build_mpi_gpu/solve_square_facade_mpi_gpu                 # 2-rank GPU
```

`EXASIM_DIR` must point at the source/install tree (it tells
`setcommonstruct` where to find `backend/Preprocessing/{master,gauss}nodes.bin`):

```bash
EXASIM_DIR=$EXASIM build_cpu/solve_square_facade
```

Working directory ends with zero files written. The runs are
numerically equivalent across backends — Newton converges
`0.211548 → ≈5e-12` in 1 iteration, max\|udg\| = π (the analytical
gradient of `sin(πx)sin(πy)`).

---

## 5. External CMake project — full standalone application

This section walks a complete `solve_square_facade` consumer
**outside** the Exasim tree. The directory layout:

```
my_solver/
├── CMakeLists.txt
└── solve_square_facade.cpp     # ~170 lines, see apps/library_example/solve_square/
```

That's it. No `pdeapp.txt`, no `pdemodel.txt`, no `grid.bin`, no
`datain/`/`dataout/`. The model, mesh, boundary classifiers, and
physics params all live in the C++ source.

### 5.1 Install Exasim

Run §2's commands once with `cmake --install` and an explicit prefix:

```bash
cmake --install build_cpu --prefix /opt/exasim
cmake --install build_gpu --prefix /opt/exasim
cmake --install build_mpi --prefix /opt/exasim
cmake --install build_mpi_gpu --prefix /opt/exasim
```

This populates:

```
/opt/exasim/include/exasim/      # public template headers + run.hpp + solver_facade.hpp
/opt/exasim/include/backend/     # backend headers exasim/* depends on
/opt/exasim/lib/cmake/Exasim/    # ExasimConfig.cmake (find_package target)
/opt/exasim/bin/text2code        # codegen tool (only needed for §6 path)
```

The four Kokkos prefixes (`/opt/exasim/external/kokkos/build*`) need
to be reachable via `find_package(Kokkos PATHS …)`. If you used a
different layout, adjust the paths below.

### 5.2 `my_solver/CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_solver CXX)
set(CMAKE_CXX_STANDARD 17)

option(USE_GPU "build for CUDA/HIP" OFF)
option(USE_MPI "build for MPI"      OFF)

# Kokkos: pick the prefix matching the backend.
if(USE_GPU)
    set(KOKKOS_PREFIX "/opt/exasim/external/kokkos/buildcuda")
else()
    set(KOKKOS_PREFIX "/opt/exasim/external/kokkos/buildserial")
endif()
find_package(Kokkos REQUIRED PATHS ${KOKKOS_PREFIX}/lib/cmake/Kokkos NO_DEFAULT_PATH)

if(USE_MPI)
    find_package(MPI REQUIRED)
endif()

find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Exasim REQUIRED PATHS /opt/exasim/lib/cmake/Exasim NO_DEFAULT_PATH)

add_executable(my_solver solve_square_facade.cpp)

if(USE_GPU)
    set_source_files_properties(solve_square_facade.cpp PROPERTIES
        COMPILE_OPTIONS "-ffast-math;-O3;-DNDEBUG;-fno-unroll-loops")
endif()

target_link_libraries(my_solver PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# `_TEXT2CODE` enables the templated-kernel path inside
# `<exasim/run.hpp>` even though no codegen actually runs.
target_compile_definitions(my_solver PRIVATE _TEXT2CODE)

if(USE_MPI)
    target_compile_definitions(my_solver PRIVATE _MPI)
    target_link_libraries(my_solver PRIVATE MPI::MPI_CXX)
endif()
if(USE_GPU)
    target_compile_definitions(my_solver PRIVATE _CUDA)   # or _HIP
endif()

# `libpdemodel{serial,cuda,hip}.so` must be on the rpath if your
# CSolution<M> ever falls through to the `AbiAdapter` branch (it
# doesn't, in the programmatic path — but the link line is the
# same). Point at the same Model_LIB_DIR Exasim's own builds use.
target_link_directories(my_solver PRIVATE /opt/exasim/lib)
set_target_properties(my_solver PROPERTIES BUILD_RPATH /opt/exasim/lib)
```

`Exasim::headers` is INTERFACE-only — the user's pointwise math
from `MyModel` is inlined into the templated kernels at *the
consumer's* compile time, so each backend variant gets fully
specialized code.

### 5.3 Build the four variants

```bash
cd my_solver

cmake -B build_cpu       -DUSE_GPU=OFF -DUSE_MPI=OFF
cmake -B build_gpu       -DUSE_GPU=ON  -DUSE_MPI=OFF
cmake -B build_mpi       -DUSE_GPU=OFF -DUSE_MPI=ON
cmake -B build_mpi_gpu   -DUSE_GPU=ON  -DUSE_MPI=ON

PATH=/usr/local/cuda/bin:/opt/exasim/openmpi/bin:$PATH cmake --build build_cpu      -j16
PATH=/usr/local/cuda/bin:/opt/exasim/openmpi/bin:$PATH cmake --build build_gpu      -j8
PATH=/usr/local/cuda/bin:/opt/exasim/openmpi/bin:$PATH cmake --build build_mpi      -j16
PATH=/usr/local/cuda/bin:/opt/exasim/openmpi/bin:$PATH cmake --build build_mpi_gpu  -j8
```

### 5.4 Run

`EXASIM_DIR` points at the install prefix so the runtime can locate
`backend/Preprocessing/{master,gauss}nodes.bin`:

```bash
EXASIM_DIR=/opt/exasim ./build_cpu/my_solver
EXASIM_DIR=/opt/exasim ./build_gpu/my_solver
EXASIM_DIR=/opt/exasim mpirun -np 2 ./build_mpi/my_solver
EXASIM_DIR=/opt/exasim mpirun -np 2 ./build_mpi_gpu/my_solver
```

Each writes nothing to the working directory (`saveOutputs = 0` by
default in the façade); read the converged state from the `udg()`
accessor inside `main()`.

> **MPI status (HOT.7.7):** the in-memory ABI is wired for CPU + GPU.
> MPI variants currently error out with a "use file path" message
> because `meshFromArrays` produces a serial mesh and ParMETIS
> expects distributed metadata (`elemGlobalID`, `elmdist`,
> `ne_global`). The proper fix — a `solver.set_mesh_distributed(…)`
> API where each rank provides its slice and ParMETIS repartitions —
> is queued as HOT.7.8. For multi-rank runs today, drop down to
> the legacy file-driven `CSolution<M>(filein, ...)` constructor
> (§6) and use `mpirun` with the `_codegen_mpi` /
> `_codegen_mpi_gpu` binaries.

---

## 6. Convenience — `text2code` codegen path

If you'd rather author the math symbolically, `text2code` consumes
a `pdemodel.txt` (SymEngine DSL) and emits the same `MyModel` struct
the §3 path expects. This is what most existing apps in
`apps/poisson/` and `apps/navierstokes/` use. Both paths produce
identical solver behavior — choose by preference.

```
my_app/
├── pdemodel.txt        # SymEngine DSL
├── pdeapp.txt          # runtime config (boundary conditions etc.)
├── grid.bin            # mesh in legacy binary format
└── main.cpp            # 3-line entrypoint: exasim::run<GeneratedModel>(argc, argv)
```

```cpp
// main.cpp — entire body:
#include <exasim/run.hpp>
#include "my_model.hpp"      // generated by text2code
int main(int argc, char** argv) {
    return exasim::run<GeneratedModel>(argc, argv);
}
```

```bash
text2code my_app/pdeapp.txt        # writes my_model.hpp, app.bin, ...
build_cpu/my_app pdeapp.txt        # runs and writes outudg_np0.bin etc.
```

The runtime here goes through `pdeapp.txt` parsing, ParMETIS
partitioning of `grid.bin`, the `datain/{app,master,mesh,sol}.bin`
ABI between preprocessing and solver, and writes outputs to
`dataout/`. It does the same numerical work as the §3 path, just
with file-mediated handoffs in between.

When the codegen path is what you want:

- pointwise math is easier to write in `pdemodel.txt` (no need to
  hand-derive Jacobians; SymEngine differentiates symbolically)
- you want the legacy `dataout/*.bin` files for downstream tooling
  (Paraview converters etc.)
- you're porting an existing pdeapp.txt-based workflow without
  changes

When the §3 programmatic path is what you want:

- you're embedding Exasim in another program (UQ outer loop,
  optimization, parameter sweep) and want to call `solver.solve()`
  in a hot loop
- you need the converged state in memory rather than on disk
- you want a typed boundary-classification API
- you don't want the SymEngine + text2code build dependency

---

## 7. Common gotchas

- **MPI+GPU segfault in `mca_btl_vader`**: OpenMPI was not built
  `--with-cuda`. See §1b. Symptom: address `0x7f…` in the failing
  signal handler is a CUDA/HIP virtual address.
- **`nvcc: command not found` from `nvcc_wrapper`**: PATH doesn't
  have `/usr/local/cuda/bin`. Set it before *configure* and *build*;
  it has to be visible to subprocess invocations of `nvcc_wrapper`.
- **All ranks bind to GPU 0 on multi-GPU nodes**: the binary's
  per-rank `cudaSetDevice(shmrank % deviceCount)` requires
  `MPI_COMM_TYPE_SHARED` to actually return a per-node `shmrank`.
  If your MPI launcher gives every rank the same node-shared rank,
  set `CUDA_VISIBLE_DEVICES=$LOCAL_RANK` per rank instead.
- **Solver runs but `udg()` returns nullptr**: you forgot to call
  `solve()`, or the solve raised an exception. Check the return
  value of `solver.solve()` and the runtime's stderr.
- **`ExasimSolver::solve(): mesh not set`**: `set_mesh()` was never
  called, or `np` was 0. Confirm your `nv`/`ne` arithmetic.
- **GPU view destructor segfault at program exit**: Kokkos finalize
  ordering. Make sure `solver` (the `ExasimSolver` object) is
  destroyed *before* `Kokkos::finalize()`. The pattern above
  scopes `solver` inside an inner block to enforce this.
- **Curved-boundary projection moves nodes off the boundary by a
  large amount**: your level-set function probably doesn't have
  zero gradient norm at the surface. Newton step is `x - f·∇f/|∇f|²`
  — when |∇f| is tiny, the step explodes. Pick a function with
  bounded gradient near the surface, or run multiple projector
  iterations (Exasim runs once; for sharper curvature, project to
  the boundary in pre-processing).
