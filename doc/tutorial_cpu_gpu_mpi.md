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

The interesting case is consuming Exasim from a separate project — your
PDE, your `pdemodel.txt`, your `main.cpp`, building against an
installed Exasim.

Project layout:

```
my_solver/
├── CMakeLists.txt
├── pdemodel.txt        # SymEngine DSL — your PDE
├── pdeapp.txt          # runtime config
├── grid.bin            # mesh
└── main.cpp            # 3 lines (see below)
```

### `main.cpp`

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"     // generated by text2code from pdemodel.txt

int main(int argc, char** argv) {
    return exasim::run<GeneratedModel>(argc, argv);
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
