# 1. Installation

Build prerequisites, vendored dependencies, and the four Exasim
backend variants. Tested on Ubuntu 22.04 with CUDA 12.x; HIP follows
the same pattern with `kokkos/buildhip` substituted for
`kokkos/buildcuda`.

After cloning, the repository folder MUST be named `Exasim`. If you
renamed it, rename it back — Kokkos will not build properly when the
path contains spaces or differs from the expected layout.

## System packages

```bash
# Ubuntu 22.04
sudo apt install -y build-essential cmake git \
                    libblas-dev liblapack-dev \
                    libmkl-dev    # optional, faster
```

For GPU builds install the CUDA toolkit (`nvcc` on PATH) or ROCm
(`hipcc` on PATH) ahead of time.

## Vendored dependencies

Exasim ships its own known-good builds of Kokkos, METIS / ParMETIS,
SymEngine, and text2code. You don't need a system Kokkos.

### OpenMPI with CUDA awareness (only for MPI+GPU)

For MPI + GPU runs OpenMPI must be built `--with-cuda` so ranks can
pass device pointers directly. The system OpenMPI almost certainly
is not built that way.

```bash
EXASIM=$PWD                    # repo root
MPI_PREFIX=$EXASIM/openmpi     # or wherever you want it

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

For CPU-only MPI, drop `--with-cuda=…` and the CUDA `PATH` prefix.

### Kokkos (serial + CUDA prefixes)

```bash
cd $EXASIM/kokkos

# Serial — CPU + MPI builds use this
mkdir buildserial && cd buildserial
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildserial \
         -DKokkos_ENABLE_SERIAL=ON
make -j16 install
cd ..

# CUDA — GPU + MPI+GPU builds use this. ARCH_VOLTA70 = V100; pick
# your arch from kokkos/cmake/kokkos_arch.cmake.
mkdir buildcuda && cd buildcuda
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildcuda \
         -DKokkos_ENABLE_SERIAL=ON \
         -DKokkos_ENABLE_CUDA=ON \
         -DKokkos_ARCH_VOLTA70=ON \
         -DCMAKE_CXX_COMPILER=$EXASIM/kokkos/buildcuda/bin/nvcc_wrapper
make -j16 install
```

For HIP, the second build uses `-DKokkos_ENABLE_HIP=ON
-DKokkos_ARCH_VEGA90A=ON` (or your arch) and
`-DCMAKE_CXX_COMPILER=hipcc`.

### METIS / GKlib / ParMETIS

```bash
cd $EXASIM/metis

cd GKlib && make config prefix=$EXASIM/metis/GKlib && \
            make -j8 install && cd ..
cd METIS && make config gklib_path=$EXASIM/metis/GKlib \
                        prefix=$EXASIM/metis/METIS \
                        shared=1 r64=1 i64=1 && \
            make -j8 install && cd ..

# Only needed for MPI variants.
PATH=$MPI_PREFIX/bin:$PATH cd ParMETIS && \
    make config gklib_path=$EXASIM/metis/GKlib \
                metis_path=$EXASIM/metis/METIS \
                prefix=$EXASIM/metis/ParMETIS \
                cc=mpicc cxx=mpicxx && \
    make -j8 install
```

### SymEngine + text2code (only for the codegen path)

You only need these if you're going to use the codegen authoring
path (write `pdemodel.txt`, generate `my_model.hpp` with text2code).
The hand-written-`Model` path doesn't link SymEngine or text2code.

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

The text2code binary lands at `$EXASIM/build/text2code` (the
`CMAKE_RUNTIME_OUTPUT_DIRECTORY` setting in the text2code project
points there explicitly so the rest of the harness can find it).

## Build Exasim

Four CMake configure lines, one per backend. Same source code; only
the flags differ.

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

Each build dir produces dozens of binaries — the in-tree example
apps (`apps/library_example/<name>_codegen/main.cpp` and
`<name>_facade` variants) plus the legacy `cputc2EXASIM`,
`gpuEXASIM`, and friends.

## Install system-wide

To consume Exasim as a header-only library from another project:

```bash
cmake --install build_cpu --prefix /opt/exasim
cmake --install build_gpu --prefix /opt/exasim   # adds GPU libs
```

Layout:

```
/opt/exasim/
├── include/exasim/          # public headers
├── lib/cmake/exasim/        # find_package(Exasim) machinery
└── lib/                     # text2code-generated dylibs (optional)
```

Consumers do `find_package(Exasim REQUIRED)` and link
`Exasim::headers Kokkos::kokkos MPI::MPI_CXX`. See
[`03-methods/hand-written-model.md`](03-methods/hand-written-model.md)
for the consumer CMake template.

## Verification

After building any variant, run the in-tree validation:

```bash
ctest --test-dir build_cpu -j 8 -R "poisson2d:.*:cpu"
```

Should report all 4 `poisson2d` tests pass in under a minute. If the
build dir doesn't have any tests, your configure didn't pass
`-DEXASIM_BUILD_LIBRARY_EXAMPLES=ON`. See
[`04-internals/testing.md`](04-internals/testing.md) for the
ctest matrix details.

## Troubleshooting

- **`Kokkos_DIR not found`** — the install prefix you pointed at
  doesn't have `lib/cmake/Kokkos/`. Verify the Kokkos build's
  `make install` ran cleanly and you used the right `--prefix`.
- **`Could NOT find ParMETIS`** — pass `-DPARMETIS_ROOT=$EXASIM/metis/ParMETIS`
  to your `cmake -S install -B build_mpi` line, or check that the
  ParMETIS make ran with the matching MPI compiler.
- **Build hangs at `Building CXX object` for >5 min on a small file** —
  almost always Kokkos's nvcc_wrapper hitting an OOM. Drop
  `-j16` to `-j4` and retry.
- **Linker can't find `libpdemodelserial.so`** — text2code wasn't
  run for the example you're building. Run
  `bash apps/library_example/regenerate.sh <name>` first.
