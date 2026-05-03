# Tutorial

Solve the same problem (Poisson 2D on the unit square, manufactured
solution `u = sin(πx)·sin(πy)`) through every authoring × driving
combination.

| Path | Authoring | Driving | Walkthrough |
|---|---|---|---|
| 01 | hand-written `Poisson2D` struct | `exasim::run<M>(argc, argv)` (legacy CLI) | [`01-handwritten-cli/`](01-handwritten-cli/README.md) |
| 02 | hand-written `Poisson2D` struct | `ExasimSolver<M>` (embedded) | [`02-handwritten-facade/`](02-handwritten-facade/README.md) |
| 03 | `pdemodel.txt` + text2code | `exasim::run<GeneratedModel>(argc, argv)` | [`03-codegen-cli/`](03-codegen-cli/README.md) |
| 04 | `pdemodel.txt` + text2code | `ExasimSolver<GeneratedModel>::load_pdeapp` | [`04-codegen-facade/`](04-codegen-facade/README.md) |
| 05 | `pdemodel.txt` + text2code | `cput2cEXASIM` (legacy AbiAdapter, deprecated) | [`05-codegen-abi-adapter/`](05-codegen-abi-adapter/README.md) |

Every path produces the same converged solution `udg`. The
`tutorial/run-all.sh` script runs all four supported paths
(01–04) and verifies they agree to within the element-L2
threshold.

## Prerequisites

Build Exasim per [`../docs/01-installation.md`](../docs/01-installation.md).
Each tutorial path uses the corresponding example apps shipped in
`apps/library_example/`. Run the tutorial test with:

```bash
ctest --test-dir build_cpu -R tutorial
```

## Backend variants

Each path's `README.md` shows the CPU run. Substitute the matching
binary suffix for GPU / MPI / MPI+GPU:

| Backend | Binary suffix |
|---|---|
| CPU | (none) |
| GPU | `_gpu` |
| MPI | `_mpi` (run with `mpirun -np N`) |
| MPI+GPU | `_mpi_gpu` (run with `mpirun -np N`) |

See [`../docs/03-methods/cpu-gpu-mpi-mpigpu.md`](../docs/03-methods/cpu-gpu-mpi-mpigpu.md)
for the configure lines per backend.
