# 0. Introduction

Exasim is a C++ library for solving partial differential equations
(PDEs) with high-order discontinuous Galerkin (DG) methods. It
combines a templated solver core (HDG / LDG, GMRES, Newton, DIRK
time stepping) with optional symbolic code generation and Kokkos-
based backend portability, so the same model definition runs on
CPU, GPU, MPI, and MPI+GPU without changing the math.

## What you get

- **Spatial discretizations** — local DG and hybridized DG; hex /
  tet / quad / tri elements; arbitrary polynomial order.
- **Time stepping** — diagonally implicit Runge–Kutta (DIRK) with
  configurable stage count.
- **Solvers** — Newton outer + GMRES inner, multiple
  preconditioners (block ILU, additive Schwarz, RB acceleration).
- **Multi-physics** — monolithic HDG with multi-domain coupling.
- **Backend portability** — Kokkos picks CPU / NVIDIA GPU (CUDA) /
  AMD GPU (HIP) at compile time.
- **Distributed-memory** — MPI partitioning via ParMETIS.

## How you use it

Two independent dimensions: how you **author** the math, and how you
**drive** the solver.

### Authoring

| Path | What you write | When to pick |
|---|---|---|
| **Hand-written `Model`** | A C++ struct with `KOKKOS_INLINE_FUNCTION static` methods (Flux, Source, BCs, ICs, …) | Full control, custom math, hand-derived Jacobians, no DSL learning curve |
| **Symbolic codegen** | `pdemodel.txt` in the SymEngine DSL; text2code emits `my_model.hpp` | Concise math, automatic CSE, derivative generation (Jacobians for free) |

Both produce the same `my_model.hpp` shape — the templated core
doesn't know which path generated it.

### Driving

| Path | Entry point | When to pick |
|---|---|---|
| **Legacy CLI** | `exasim::run<Model>(argc, argv)` reading `pdeapp.txt` from disk | Standalone simulation binaries; HPC job scripts |
| **Embedded library** | `ExasimSolver<Model>` — set mesh / boundaries / params programmatically; read solution from accessors | Optimization loops, reduced-order modeling, coupled multi-physics, anywhere Exasim is one component of a larger C++ program |

### Backends

Compile-time choice driven by CMake flags:

| Variant | Kokkos build | EXASIM flags |
|---|---|---|
| CPU | `kokkos/buildserial` | `EXASIM_NOMPI=ON` |
| GPU | `kokkos/buildcuda` (or `buildhip`) | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON` |
| MPI | `kokkos/buildserial` | `EXASIM_MPI=ON` |
| MPI+GPU | `kokkos/buildcuda` | `EXASIM_MPI=ON, EXASIM_CUDA=ON` |

The same source code produces all four binaries — only the build
configuration differs. See
[`02-methods-overview.md`](02-methods-overview.md) for the full
matrix of authoring × driving × backend combinations.

## Repository layout

```
Exasim/
├── apps/                    # in-tree example PDEs
│   ├── poisson/             # poisson2d, poisson3d, periodic, ...
│   ├── navierstokes/        # naca0012steady, cone, sharpb2, ...
│   └── library_example/     # codegen + facade test apps + harness
├── backend/                 # implementation: Common, Discretization,
│   ├── Common/              # Preprocessing, Solution, Postprocessing,
│   ├── Discretization/      # Model
│   ├── Preprocessing/
│   ├── Solution/
│   └── Postprocessing/
├── include/exasim/          # public API (headers consumed by users)
├── kokkos/                  # vendored Kokkos
├── metis/                   # vendored METIS + ParMETIS + GKlib
├── text2code/               # vendored SymEngine + text2code generator
├── baseline/                # recorded reference outputs for tests
├── docs/                    # this documentation tree
└── doc/                     # legacy docs (some still authoritative)
```

## Where to go from here

- **Set up your environment** → [01-installation.md](01-installation.md)
- **Pick how you'll use Exasim** → [02-methods-overview.md](02-methods-overview.md)
- **Write your first solver** → [03-methods/hand-written-model.md](03-methods/hand-written-model.md)
- **Embed Exasim in a larger program** → [03-methods/embedded-facade.md](03-methods/embedded-facade.md)
- **Generate models symbolically** → [03-methods/codegen-text2code.md](03-methods/codegen-text2code.md)
- **Run on CPU / GPU / MPI** → [03-methods/cpu-gpu-mpi-mpigpu.md](03-methods/cpu-gpu-mpi-mpigpu.md)
