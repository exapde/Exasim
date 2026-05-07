# Installation — Linux x86_64 with AMD GPU

For machines with AMD GPUs running ROCm. Adds HIP and MPI+HIP
variants on top of the CPU baseline.

## System packages

```bash
sudo apt install -y build-essential cmake git \
                    libblas-dev liblapack-dev
```

ROCm: install via the AMD ROCm repository for your distribution.
After install, `hipcc --version` should report your ROCm version.

## ROCm-aware OpenMPI

Required for MPI+HIP. Build OpenMPI with UCX + ROCm support so ranks
can pass device pointers directly:

```bash
EXASIM=$PWD
MPI_PREFIX=$EXASIM/openmpi

mkdir -p $EXASIM/build_openmpi && cd $EXASIM/build_openmpi
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.gz
tar xf openmpi-4.1.6.tar.gz && cd openmpi-4.1.6

./configure \
    --prefix=$MPI_PREFIX \
    --with-rocm=/opt/rocm \
    --with-ucx \
    --disable-mpi-fortran --disable-oshmem
make -j16 install

export PATH=$MPI_PREFIX/bin:$PATH
ompi_info --param btl smcuda --level 9   # ROCm/UCX flags appear here
```

UCX (built with ROCm support) handles the GPU-aware transport
layer. If your system already has UCX with ROCm, OpenMPI picks it
up automatically. For CPU-only MPI, omit `--with-rocm` and `--with-ucx`.

## Vendored dependencies

Follow [`common.md`](common.md) for GKlib + METIS + ParMETIS,
Kokkos serial, SymEngine + text2code.

## Build Kokkos HIP

```bash
cd $EXASIM/kokkos
mkdir buildhip && cd buildhip
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildhip \
         -DKokkos_ENABLE_SERIAL=ON \
         -DKokkos_ENABLE_HIP=ON \
         -DKokkos_ARCH_VEGA90A=ON \
         -DCMAKE_CXX_COMPILER=hipcc
make -j16 install
```

`Kokkos_ARCH_*` for AMD GPUs:

| Flag | Hardware |
|---|---|
| `Kokkos_ARCH_VEGA90A` | MI100, MI250 |
| `Kokkos_ARCH_AMD_GFX942` | MI300A, MI300X |
| `Kokkos_ARCH_NAVI1030` | RX 6900-series |

Full list in `kokkos/cmake/kokkos_arch.cmake`.

## Build Exasim

GPU (single-rank):

```bash
cd $EXASIM
cmake -S install -B build_gpu \
    -DEXASIM_HIP=ON -DEXASIM_NOMPI=ON \
    -DCMAKE_CXX_COMPILER=hipcc \
    -DKokkos_DIR=$EXASIM/kokkos/buildhip/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
cmake --build build_gpu -j8
```

CPU + MPI: same as in [`linux-cpu.md`](linux-cpu.md), but using the
ROCm-aware `MPI_PREFIX` if you want to share one MPI install.

GPU + MPI:

```bash
PATH=$MPI_PREFIX/bin:$PATH \
    cmake -S install -B build_mpi_gpu \
    -DEXASIM_HIP=ON -DEXASIM_MPI=ON -DEXASIM_NOMPI=OFF \
    -DWITH_PARMETIS=ON \
    -DCMAKE_CXX_COMPILER=hipcc \
    -DKokkos_DIR=$EXASIM/kokkos/buildhip/lib/cmake/Kokkos \
    -DEXASIM_BUILD_LIBRARY_EXAMPLES=ON
PATH=$MPI_PREFIX/bin:$PATH cmake --build build_mpi_gpu -j8
```

## Verification

```bash
ctest --test-dir build_gpu     -R "poisson2d:.*:gpu"     -j 4
ctest --test-dir build_mpi_gpu -R "poisson2d:.*:mpi_gpu" -j 4
```
