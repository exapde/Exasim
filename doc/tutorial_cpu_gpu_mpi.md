# Tutorial — CPU / GPU / MPI from scratch on a fresh machine

End-to-end walkthrough that builds Exasim, generates a model with
`text2code`, and runs it on every supported backend:

| variant            | what it is                          | binary suffix       |
| ------------------ | ----------------------------------- | ------------------- |
| CPU                | single-rank serial                  | `_codegen`          |
| GPU                | single-rank CUDA or HIP (Kokkos)    | `_codegen_gpu`      |
| MPI                | multi-rank serial CPU               | `_codegen_mpi`      |
| MPI + GPU          | multi-rank GPU (one rank per device)| `_codegen_mpi_gpu`  |

Tested on Ubuntu 22.04 with CUDA 12.x on V100s. HIP follows the same
script with `kokkos/buildhip` and `EXASIM_HIP=ON` substituted.

---

## 1. Install deps

Three deps are vendored in-tree (`kokkos/`, `text2code/symengine/`,
`metis/`); two are external (`BLAS/LAPACK`, optional `MPI`, optional
`CUDA`/`ROCm`). The vendored builds are how Exasim ships its known-good
versions — you don't need a system Kokkos.

### 1a. System packages

```bash
# Ubuntu 22.04
sudo apt install -y build-essential cmake git \
                    libblas-dev liblapack-dev \
                    libmkl-dev    # optional but faster
```

For GPU builds, install the CUDA toolkit (`nvcc` on PATH) or the ROCm
toolchain (`hipcc` on PATH) ahead of time. CUDA 12.x and Kokkos 4.7 are
the tested combination.

### 1b. OpenMPI with CUDA awareness

If you're going to do MPI + GPU, OpenMPI must be built `--with-cuda` so
ranks can pass device pointers directly. The system OpenMPI almost
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
ompi_info | grep -i cuda           # should print: MPI extensions: ... cuda ...
```

For CPU-only MPI (no GPU), drop `--with-cuda=...` and `PATH=/usr/local/cuda/bin:`.

### 1c. Kokkos (vendored)

Build serial and CUDA backends as separate prefixes. The CMake
configuration picks one based on `EXASIM_CUDA` / `EXASIM_HIP`.

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

For HIP, replace the second build with `Kokkos_ENABLE_HIP=ON
Kokkos_ARCH_VEGA90A=ON` (or your arch) and `CMAKE_CXX_COMPILER=hipcc`.

### 1d. METIS / GKlib (vendored)

```bash
cd $EXASIM/metis

cd GKlib && make config prefix=$EXASIM/metis/GKlib && make -j8 install && cd ..
cd METIS && make config gklib_path=$EXASIM/metis/GKlib prefix=$EXASIM/metis/METIS \
                  shared=1 r64=1 i64=1 && make -j8 install && cd ..

# Only needed for MPI variants — does not need r64=1/i64=1.
PATH=$MPI_PREFIX/bin:$PATH cd ParMETIS && \
    make config gklib_path=$EXASIM/metis/GKlib \
                metis_path=$EXASIM/metis/METIS \
                prefix=$EXASIM/metis/ParMETIS \
                cc=mpicc cxx=mpicxx && \
    make -j8 install
```

### 1e. SymEngine (vendored, builds with text2code)

`text2code` links SymEngine statically. The bundled SymEngine builds
without GMP and uses Boost.Multiprecision so there are no external
deps to chase.

```bash
cd $EXASIM/text2code/symengine
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/text2code/symengine/build \
         -DBUILD_SHARED_LIBS=OFF \
         -DWITH_GMP=OFF \
         -DINTEGER_CLASS=boostmp
make -j16 install
```

### 1f. text2code

```bash
cd $EXASIM/text2code/text2code
cmake -S . -B build -DEXASIM_ROOT=$EXASIM
cmake --build build -j8
```

This produces `$EXASIM/build/text2code` (binary). It auto-runs once on
the bundled `pdeapp.txt` as a smoke test — that's expected.

---

## 2. Configure + install Exasim

Exasim ships as four flavors of the same template library. You pick
which by configuring different build dirs against different Kokkos
prefixes.

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

Each build produces:

| build dir          | targets built                                |
| ------------------ | -------------------------------------------- |
| `build_cpu/`       | `<name>_codegen` for every example           |
| `build_gpu/`       | `<name>_codegen_gpu`                         |
| `build_mpi/`       | `<name>_codegen_mpi`                         |
| `build_mpi_gpu/`   | `<name>_codegen_mpi_gpu`                     |

To install Exasim system-wide:

```bash
cmake --install build_cpu --prefix /opt/exasim
cmake --install build_gpu --prefix /opt/exasim   # adds GPU libs to same prefix
```

After `--install`, the prefix has:

```
/opt/exasim/
├── bin/text2code           # generates my_model.hpp from pdemodel.txt
├── include/exasim/         # public template headers
├── include/backend/        # internal headers exasim/* depends on
└── lib/cmake/Exasim/       # ExasimConfig.cmake for find_package
```

---

## 3. Run a model

Pick `poisson2d` (smallest example, 2D Poisson with manufactured
solution). Each example has a runnable `pdeapp.txt` driving the
solver.

```bash
cd $EXASIM/apps/library_example/poisson2d_codegen

# Generate my_model.hpp from pdemodel.txt
$EXASIM/build/text2code pdeapp.txt          # or `text2code pdeapp.txt` if installed

# Or use the convenience driver
bash $EXASIM/apps/library_example/regenerate.sh poisson2d
```

Now run any of the four variants. They all consume the same `pdeapp.txt`:

```bash
# CPU
$EXASIM/build_cpu/poisson2d_codegen pdeapp.txt

# GPU
$EXASIM/build_gpu/poisson2d_codegen_gpu pdeapp.txt

# MPI (e.g. 2 ranks)
PATH=$MPI_PREFIX/bin:$PATH \
    mpirun -np 2 $EXASIM/build_mpi/poisson2d_codegen_mpi pdeapp.txt

# MPI + GPU (one rank per GPU; for 2 GPUs use -np 2)
PATH=$MPI_PREFIX/bin:$PATH \
    mpirun -np 2 $EXASIM/build_mpi_gpu/poisson2d_codegen_mpi_gpu pdeapp.txt
```

Each writes `dataout/outudg_np<rank>.bin` (the discrete state),
`dataout/outuhat_np<rank>.bin` (HDG trace), and `dataout/outqoi.txt`
(scalar QoI integrand if the model defines `qoi_volume`/`qoi_boundary`).

The bundled validation harness exercises every variant against the
recorded baseline:

```bash
cd $EXASIM
bash apps/library_example/validate_codegen.sh                     # CPU
bash apps/library_example/validate_codegen.sh --variant gpu      --build build_gpu
bash apps/library_example/validate_codegen.sh --variant mpi      --build build_mpi      --np 2
bash apps/library_example/validate_codegen.sh --variant mpi_gpu  --build build_mpi_gpu  --np 2
```

---

## 4. External app via CMake

The interesting case is consuming Exasim from a separate project,
building against an installed Exasim.

The runtime needs *some* type that satisfies the `Model` contract from
`<exasim/model.hpp>` — that's it. You have two ways to produce one:

- **Codegen path** — write `pdemodel.txt` (SymEngine DSL), run
  `text2code` to emit `my_model.hpp`. This is the path most apps in
  `apps/poisson/` and `apps/navierstokes/` take.
- **Hand-written path** — write `my_model.hpp` directly. Inherit
  `exasim::ModelDefaults<Self>` and override only the methods your
  PDE needs. See `doc/getting_started.md` for a worked example. No
  `pdemodel.txt`, no `text2code`, no SymEngine at runtime.

Both paths produce the same `Model` template type that
`exasim::run<M>(argc, argv)` instantiates. The choice is purely about
how *you* prefer to author the pointwise math.

Project layout (codegen path):

```
my_solver/
├── CMakeLists.txt
├── pdemodel.txt        # SymEngine DSL — your PDE (codegen path only)
├── pdeapp.txt          # runtime config (always needed)
├── grid.bin            # mesh
└── main.cpp            # 3 lines (see below)
```

For the hand-written path, drop `pdemodel.txt` and add `my_model.hpp`
to the source tree directly.

### `main.cpp`

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"     // hand-written, OR generated by text2code

int main(int argc, char** argv) {
    return exasim::run<GeneratedModel>(argc, argv);
    // ↑ name matches whatever `struct ... : ModelDefaults<...>` you
    // wrote in my_model.hpp. text2code emits `GeneratedModel`;
    // hand-written code can name it whatever (Poisson2D, NSEuler, etc.).
}
```

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)
project(my_solver CXX)
set(CMAKE_CXX_STANDARD 17)

option(USE_GPU "build for CUDA/HIP" OFF)
option(USE_MPI "build for MPI"      OFF)

# Pick the matching Kokkos prefix for your backend.
if(USE_GPU)
    find_package(Kokkos REQUIRED PATHS /opt/exasim/external/kokkos/buildcuda
                 NO_DEFAULT_PATH)
else()
    find_package(Kokkos REQUIRED PATHS /opt/exasim/external/kokkos/buildserial
                 NO_DEFAULT_PATH)
endif()

if(USE_MPI)
    find_package(MPI REQUIRED)
endif()

find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Exasim REQUIRED PATHS /opt/exasim NO_DEFAULT_PATH)

# Run text2code at configure time to regenerate my_model.hpp from pdemodel.txt.
add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/my_model.hpp
    COMMAND /opt/exasim/bin/text2code ${CMAKE_CURRENT_SOURCE_DIR}/pdeapp.txt
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pdemodel.txt
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(gen_model DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/my_model.hpp)

add_executable(my_solver main.cpp)
add_dependencies(my_solver gen_model)
target_include_directories(my_solver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(my_solver PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(USE_MPI)
    target_compile_definitions(my_solver PRIVATE _MPI)
    target_link_libraries(my_solver PRIVATE MPI::MPI_CXX)
endif()

if(USE_GPU)
    target_compile_definitions(my_solver PRIVATE _CUDA)   # or _HIP
endif()
```

### Build matrix

```bash
cmake -B build_cpu       -DUSE_GPU=OFF -DUSE_MPI=OFF
cmake -B build_gpu       -DUSE_GPU=ON  -DUSE_MPI=OFF
cmake -B build_mpi       -DUSE_GPU=OFF -DUSE_MPI=ON
cmake -B build_mpi_gpu   -DUSE_GPU=ON  -DUSE_MPI=ON
cmake --build build_cpu -j16
# ... same for the others
```

### Run

```bash
./build_cpu/my_solver                                    pdeapp.txt
./build_gpu/my_solver                                    pdeapp.txt
mpirun -np 4 ./build_mpi/my_solver                       pdeapp.txt
mpirun -np 4 ./build_mpi_gpu/my_solver                   pdeapp.txt   # 4 GPUs
```

`Exasim::headers` is INTERFACE-only — no Exasim runtime library exists.
The user's pointwise math from `my_model.hpp` is inlined into the
templated kernels at the consumer's compile time.

---

## 5. No text2code: drive Exasim internals from your `main()`

The codegen path is convenient but optional. Nothing in the runtime
requires `text2code` at build or run time — it's a code generator,
not a runtime dependency.

This section walks the *fully embedded* path: one `main()` that

1. generates a square quad mesh and writes it as `grid.txt`,
2. populates the `PDE` / `InputParams` / `ParsedSpec` structs in C++,
3. drives `CPreprocessing` and `CSolution<MyModel>` directly (the
   same internals `exasim::run<>` calls — without going through
   `<exasim/run.hpp>`),
4. reads the QoI back.

The hand-written model (`MyModel`) lives in the same translation unit;
no external `pdeapp.txt`, no `pdemodel.txt`, no `text2code` invocation,
no separate mesh-builder tool. The *only* file the program writes
on its way to the solve is the mesh.

### `solve_square.cpp`

```cpp
// solve_square.cpp — single-file Exasim consumer.
//
// Builds a Cartesian quad mesh, lays it out as the runtime expects,
// then drives CPreprocessing + CSolution<Poisson2D> directly. This is
// what <exasim/run.hpp> does internally; we replace it with code we
// can edit, instrument, or wrap in a tighter outer loop.

#include <Kokkos_Core.hpp>

#include <exasim/preprocessing.hpp>      // CPreprocessing
#include <exasim/solution.hpp>           // CSolution<M>
#include <exasim/model.hpp>              // exasim::ModelDefaults<>

#include <filesystem>
#include <fstream>
#include <string>

// ------------------------------------------------------------------
// 1. The PDE: Poisson 2D, hand-written. No text2code involved.
//    See doc/getting_started.md for a full annotation of every method.
// ------------------------------------------------------------------
struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;
    static constexpr int ncu    = 1;
    static constexpr int ncw    = 0;
    static constexpr int nco    = 0;
    static constexpr int nparam = 1;
    static constexpr auto disc  = exasim::Discretization::HDG;
    static constexpr int Nq     = ncu * (1 + nd);

    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double[], const double uq[],
              const double[], const double[], const double mu[],
              const double[], double) {
        f[0] = mu[0] * uq[1];
        f[1] = mu[0] * uq[2];
    }
    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double[],
                const double[], const double[], const double[],
                const double[], double) {
        constexpr double pi = 3.141592653589793;
        s[0] = 2.0 * pi * pi
             * Kokkos::sin(pi * x[0]) * Kokkos::sin(pi * x[1]);
    }
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double[], const double[],
                     const double[], const double[], const double mu[],
                     const double[], double) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }
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
    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double[], const double[], const double[]) {
        ui[0] = 0.0;
    }
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
              const double[], double) {
        ub[0] = 0.0;
    }
};

// ------------------------------------------------------------------
// 2. Generate an n×n quad mesh on the unit square and write it in the
//    text format CPreprocessing accepts. Format:
//        nd np nve ne
//        <np × nd doubles>     vertex coords, column-major
//        <ne × nve ints>       element conn, 1-based, column-major
//    Quad corners CCW; nve=4 → Exasim infers elemtype=1.
// ------------------------------------------------------------------
static void write_square_mesh(const std::string& path, int n) {
    int nv = (n + 1) * (n + 1);
    int ne = n * n;
    std::ofstream o(path);
    o << "2 " << nv << " 4 " << ne << "\n";
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i)
            o << double(i) / n << " " << double(j) / n << "\n";
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            int v00 = j * (n + 1) + i + 1;
            o << v00 << " " << v00 + 1 << " "
              << v00 + 1 + (n + 1) << " " << v00 + (n + 1) << "\n";
        }
}

// ------------------------------------------------------------------
// 3. Build the runtime config in C++. No pdeapp.txt — we populate the
//    PDE / InputParams / ParsedSpec structs directly and hand them to
//    the programmatic CPreprocessing constructor.
// ------------------------------------------------------------------
static void build_runtime_config(PDE& pde, InputParams& params,
                                 ParsedSpec& spec)
{
    pde.discretization = "hdg";
    pde.platform       = "cpu";
    pde.datapath       = ".";
    pde.datainpath     = "./datain";
    pde.dataoutpath    = "./dataout";
    pde.exasimpath     = ".";
    pde.meshfile       = "grid.txt";
    pde.modelfile      = "";       // unused on the hand-written path
    pde.gendatain      = 1;
    pde.builtinmodelID = 1;
    pde.porder = 3;
    pde.pgauss = 6;
    pde.torder = 1;
    pde.nstage = 1;
    pde.tdep   = 0;
    pde.nd     = 2;       pde.ncu = 1;       pde.ncw = 0;
    pde.nc     = 3;       // ncu * (1 + nd) for Poisson2D
    pde.neb    = 4096;    pde.nfb = 8192;
    pde.ibs    = 1;
    pde.NewtonIter   = 20;     pde.NewtonTol   = 1e-6;
    pde.GMRESiter    = 200;    pde.GMRESrestart = 50;
    pde.GMREStol     = 1e-8;
    pde.tau          = {1.0};
    pde.dt           = {0.0};
    pde.physicsparam = {1.0};

    params.boundaryConditions = {1, 1, 1, 1};
    params.boundaryExprs      = {"abs(y)<1e-8",   "abs(x-1)<1e-8",
                                 "abs(y-1)<1e-8", "abs(x)<1e-8"};
    params.tau                = {1.0};
    params.dt                 = {0.0};
    params.physicsParam       = {1.0};

    // ParsedSpec is left empty: it's only consulted to *override* the
    // pde dimensions from a parsed pdemodel.txt. With nothing to
    // override, the values we set on `pde` above stay in effect.
    (void)spec;
}

int main(int argc, char** argv) {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    int mpiprocs = 1, mpirank = 0;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
#else
    int mpiprocs = 1, mpirank = 0;
    (void)argc; (void)argv;
#endif

    Kokkos::initialize();
    {
        std::filesystem::create_directories("datain");
        std::filesystem::create_directories("dataout");
        write_square_mesh("grid.txt", /*n=*/32);

        // Build runtime config in memory.
        PDE         pde;
        InputParams params;
        ParsedSpec  spec;
        build_runtime_config(pde, params, spec);

        // a) preprocessing — programmatic constructor, no pdeapp.txt.
        //    Reads grid.txt, writes datain/{app,master,mesh,sol}.bin.
        CPreprocessing preproc(pde, params, spec, mpirank, mpiprocs);
        if (mpiprocs == 1) preproc.SerialPreprocessing();
#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
        else preproc.ParallelPreprocessing(EXASIM_COMM_LOCAL);
#endif

        // b) instantiate the solver. CSolution<M> reads the binaries
        //    CPreprocessing just wrote.
        std::string filein  = pde.datainpath + "/";
        std::string fileout = pde.dataoutpath + "/out";
        int fileoffset = 0, gpuid = 0, backend = 0;
        CSolution<Poisson2D> sol(filein, fileout, pde.exasimpath,
                                 mpiprocs, mpirank, fileoffset,
                                 gpuid, backend, pde.builtinmodelID);

        // CSolution expects nomodels/ncarray/udgarray to be populated
        // even for a single-model run — see run.hpp around line 340.
        sol.disc.common.nomodels = 1;
        sol.disc.common.ncarray  = new Int[1]{ sol.disc.common.nc };
        sol.disc.sol.udgarray    = new dstype*[1]{ &sol.disc.sol.udg[0] };

        // c) solve.
        std::ofstream resnorms;   // unused; run.hpp only opens this
                                  // file when pde.saveResNorm == 1.
        sol.SolveProblem(resnorms, backend);

        delete[] sol.disc.common.ncarray;
        delete[] sol.disc.sol.udgarray;
    }
    Kokkos::finalize();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
```

### Build

```cmake
cmake_minimum_required(VERSION 3.16)
project(solve_square CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Kokkos REQUIRED PATHS /opt/exasim/external/kokkos/buildserial NO_DEFAULT_PATH)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Exasim REQUIRED PATHS /opt/exasim NO_DEFAULT_PATH)

add_executable(solve_square solve_square.cpp)
target_link_libraries(solve_square PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
```

```bash
cmake -B build && cmake --build build -j
./build/solve_square
```

The binary writes `dataout/outqoi.txt` (the squared L² error vs. the
manufactured solution if `qoi_volume` is defined; trivially zero
otherwise) plus `dataout/outudg_np0.bin` (the discrete state). No
`text2code`, no `pdeapp.txt`, no `pdemodel.txt`.

### What still touches the disk and why

The mesh goes to `grid.txt` because `initializeMesh()` reads its
input from a file path. The four runtime ABI files
(`datain/{app,master,mesh,sol}.bin`) are written by `CPreprocessing`
and read back by `CSolution<M>` — that's the boundary between the
preprocessor and the solver. Both could be exposed as in-memory APIs
in the future; for now the disk hop is small and only happens once
per run.

For multi-rank, the same `solve_square` binary launched under `mpirun
-np N` partitions the global mesh with ParMETIS via the
`ParallelPreprocessing(EXASIM_COMM_LOCAL)` branch above.

## Common gotchas

- **MPI+GPU segfault in `mca_btl_vader`**: OpenMPI was not built
  `--with-cuda`. See §1b. Symptom: address `0x7f...` in the failing
  signal handler is a CUDA virtual address.
- **`nvcc: command not found` from `nvcc_wrapper`**: PATH doesn't have
  `/usr/local/cuda/bin`. Set it before *configure* and *build*; it has
  to be visible to subprocess invocations of `nvcc_wrapper`.
- **`gendatain=1` overwrites your mesh**: `gendatain=1` re-runs the mesh
  generator/partitioner before each solve. Set to 0 once `datain/` is
  populated to avoid re-partitioning on every run.
- **All ranks bind to GPU 0**: Exasim's device-selection logic does
  `shmrank % deviceCount`, so `mpirun -np N --bind-to none` (or letting
  OpenMPI's default placement decide) gets one rank per GPU. If your
  scheduler gives every rank the same GPU, set
  `CUDA_VISIBLE_DEVICES=$LOCAL_RANK` per rank.
- **Segfault at exit (after solve completes)**: under multi-rank GPU,
  the Kokkos View destructors run on rank 0 after another rank has
  already torn down its CUDA context, which can dereference an
  unmapped IPC handle. The QoI / output bins are written before
  teardown, so this does not affect correctness. Suppress with
  `mpirun --output-filename` or run with one rank to confirm.
