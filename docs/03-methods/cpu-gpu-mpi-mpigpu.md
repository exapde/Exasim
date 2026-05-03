# Backend variants — CPU, GPU, MPI, MPI+GPU

Backend is a compile-time decision. The same source file produces
four binaries; each lives in its own build directory.

| Variant | Build dir | Binary suffix |
|---|---|---|
| CPU | `build_cpu/` | (none) |
| GPU | `build_gpu/` | `_gpu` |
| MPI | `build_mpi/` | `_mpi` |
| MPI+GPU | `build_mpi_gpu/` | `_mpi_gpu` |

Configure lines per platform: [`../01-installation.md`](../01-installation.md).

## Running

Single-rank:

```bash
$EXASIM/build_cpu/poisson2d_codegen ./pdeapp.txt
$EXASIM/build_gpu/poisson2d_codegen_gpu ./pdeapp.txt
```

Multi-rank:

```bash
mpirun -np N $EXASIM/build_mpi/poisson2d_codegen_mpi ./pdeapp.txt
mpirun -np N $EXASIM/build_mpi_gpu/poisson2d_codegen_mpi_gpu ./pdeapp.txt
```

The MPI binaries are linked against the MPI compiler picked at
configure time. `mpirun` from the same MPI install runs them.

## GPU device binding

Multi-rank GPU code pins one rank per device using
`MPI_COMM_TYPE_SHARED`:

```cpp
MPI_Comm shmcomm;
MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
int shmrank;
MPI_Comm_rank(shmcomm, &shmrank);

#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
#ifdef HAVE_HIP
    int nd; hipGetDeviceCount(&nd); hipSetDevice(shmrank % nd);
#endif
```

Multi-host runs with N GPUs per host map cleanly: each host's
shmcomm has N ranks, each pinning a different device.

## Building all four at once

```bash
EXASIM=$PWD
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build $EXASIM/build_cpu     -j16 &
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build $EXASIM/build_gpu     -j8  &
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build $EXASIM/build_mpi     -j16 &
PATH=/usr/local/cuda/bin:$MPI_PREFIX/bin:$PATH \
    cmake --build $EXASIM/build_mpi_gpu -j8  &
wait
```

## Validation

Run the test matrix per build dir:

```bash
ctest --test-dir build_cpu     -j 8
ctest --test-dir build_gpu     -j 8
ctest --test-dir build_mpi     -j 8
ctest --test-dir build_mpi_gpu -j 8
```

The four ctests share `apps/library_example/<name>_codegen/` for
text2code regen, so they should run sequentially.

## Tutorial

[`../../tutorial/`](../../tutorial/README.md) walks Poisson 2D
through every authoring × driving combination. The
`tutorial:run-all` ctest runs the four supported paths.
