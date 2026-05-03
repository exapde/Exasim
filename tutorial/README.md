# Tutorial

Solve the same problem (Poisson 2D on the unit square, manufactured
solution `u = sin(πx)·sin(πy)`) through every authoring × driving
combination, ordered from **most automated** (drop two text files
and run a prebuilt binary) to **least automated** (manage mesh,
boundaries, and partitioning by hand from C++).

| Section | Authoring | Driving | Walkthrough |
|---|---|---|---|
| 01 | `pdemodel.txt` + text2code | `cput2cEXASIM` (prebuilt AbiAdapter binary) | [`01-codegen-abi-adapter/`](01-codegen-abi-adapter/README.md) |
| 02 | `pdemodel.txt` + text2code | `exasim::run<GeneratedModel>(argc, argv)` (legacy CLI) | [`02-codegen-cli/`](02-codegen-cli/README.md) |
| 03 | `pdemodel.txt` + text2code | `ExasimSolver<GeneratedModel>::load_pdeapp` | [`03-codegen-facade/`](03-codegen-facade/README.md) |
| 04 | hand-written `Poisson2D` struct | `exasim::run<M>(argc, argv)` (legacy CLI) | [`04-handwritten-cli/`](04-handwritten-cli/README.md) |
| 05 | hand-written `Poisson2D` struct | `ExasimSolver<M>::set_mesh` (embedded) | [`05-handwritten-facade/`](05-handwritten-facade/README.md) |
| 06 | hand-written `Poisson2D` struct | `ExasimSolver<M>::set_mesh_distributed` (MPI embedded) | [`06-handwritten-facade-mpi/`](06-handwritten-facade-mpi/README.md) |

Every section produces the same converged solution `udg`. The
`tutorial/run-all.sh` script runs each section against the chosen
backend variant and verifies the result.

## Prerequisites

Build Exasim per [`../docs/01-installation.md`](../docs/01-installation.md).
Run the tutorial test with:

```bash
ctest --test-dir build -R tutorial
```

## Backend variants

Each section's `README.md` shows the CPU run. Substitute the matching
binary suffix for GPU / MPI / MPI+GPU:

| Backend | Binary suffix |
|---|---|
| CPU | (none) |
| GPU | `_gpu` |
| MPI | `_mpi` (run with `mpirun -np N`) |
| MPI+GPU | `_mpi_gpu` (run with `mpirun -np N`) |

Section 06 only runs under the MPI variants; the other five run
under every variant.

See [`../docs/03-methods/cpu-gpu-mpi-mpigpu.md`](../docs/03-methods/cpu-gpu-mpi-mpigpu.md)
for the configure lines per backend.
