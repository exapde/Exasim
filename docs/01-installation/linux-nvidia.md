# Installation — Linux x86_64 with NVIDIA GPU

For machines with NVIDIA GPUs. Adds GPU and MPI+GPU variants on
top of the CPU baseline.

## System packages

```bash
sudo apt install -y build-essential cmake git \
                    libblas-dev liblapack-dev
```

CUDA toolkit (12.x or later): install via the NVIDIA CUDA repository.
After install, `nvcc --version` should report your toolkit version.

## CUDA-aware OpenMPI

Required for MPI+GPU. The system OpenMPI is almost never built with
CUDA awareness; build from source:

```bash
EXASIM=$PWD
MPI_PREFIX=$EXASIM/openmpi

mkdir -p $EXASIM/build_openmpi && cd $EXASIM/build_openmpi
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz
tar xf openmpi-4.1.6.tar.gz && cd openmpi-4.1.6

PATH=/usr/local/cuda/bin:$PATH ./configure \
    --prefix=$MPI_PREFIX \
    --with-cuda=/usr/local/cuda \
    --disable-mpi-fortran --disable-oshmem
make -j16 install

export PATH=$MPI_PREFIX/bin:$PATH
ompi_info | grep -i cuda
```

`ompi_info` should list `cuda` under MPI extensions. For CPU-only
MPI, omit the CUDA bits.

## Vendored dependencies

Follow [`common.md`](common.md) for GKlib + METIS + ParMETIS,
Kokkos serial, SymEngine + text2code.

## Build Kokkos CUDA

```bash
cd $EXASIM/kokkos
mkdir buildcuda && cd buildcuda
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildcuda \
         -DKokkos_ENABLE_SERIAL=ON \
         -DKokkos_ENABLE_CUDA=ON \
         -DKokkos_ARCH_VOLTA70=ON \
         -DCMAKE_CXX_COMPILER=$EXASIM/kokkos/buildcuda/bin/nvcc_wrapper
make -j16 install
```

`Kokkos_ARCH_*` selects the GPU architecture. Common values:

| Flag | Hardware |
|---|---|
| `Kokkos_ARCH_VOLTA70` | V100 |
| `Kokkos_ARCH_AMPERE80` | A100 |
| `Kokkos_ARCH_HOPPER90` | H100 |
| `Kokkos_ARCH_TURING75` | T4, RTX 20-series |

Full list in `kokkos/cmake/kokkos_arch.cmake`.

## Build Exasim

GPU (single-rank):

```bash
cd $EXASIM
cmake -S install -B build_gpu \
    -DEXASIM_CUDA=ON -DEXASIM_NOMPI=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildcuda/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
PATH=/usr/local/cuda/bin:$PATH cmake --build build_gpu -j8
```

CPU + MPI:

```bash
PATH=$MPI_PREFIX/bin:$PATH \
    cmake -S install -B build_mpi \
    -DEXASIM_CUDA=OFF -DEXASIM_HIP=OFF \
    -DEXASIM_MPI=ON   -DEXASIM_NOMPI=OFF \
    -DWITH_PARMETIS=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildserial/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON \
    -DMPI_CXX_COMPILER=$MPI_PREFIX/bin/mpicxx \
    -DMPI_C_COMPILER=$MPI_PREFIX/bin/mpicc
PATH=$MPI_PREFIX/bin:$PATH cmake --build build_mpi -j16
```

GPU + MPI:

```bash
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake -S install -B build_mpi_gpu \
    -DEXASIM_CUDA=ON -DEXASIM_MPI=ON -DEXASIM_NOMPI=OFF \
    -DWITH_PARMETIS=ON \
    -DKokkos_DIR=$EXASIM/kokkos/buildcuda/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build build_mpi_gpu -j8
```

## Verification

```bash
ctest --test-dir build_gpu     -R "poisson2d:.*:gpu"     -j 4
ctest --test-dir build_mpi     -R "poisson2d:.*:mpi"     -j 4
ctest --test-dir build_mpi_gpu -R "poisson2d:.*:mpi_gpu" -j 4
```
