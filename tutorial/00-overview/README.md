# 00 — Overview and naming convention

Exasim has two independent decisions that the user has to make
when wiring up a solver:

1. **How is the math authored?** Either the user describes the PDE
   symbolically in `pdemodel.txt` and lets `text2code` emit the C++
   model header, or the user writes the model header
   (`my_model.hpp`) by hand against the `<exasim/model.hpp>`
   contract.
2. **How is the runtime driven?** Either the runtime is a prebuilt
   binary that consumes `pdeapp.txt` from disk, or the user writes
   their own `main.cpp` that links against the templated solver
   (with several variants for how the mesh is supplied).

The tutorial covers six combinations of these two decisions,
ordered from most automated to least automated. Each section is a
self-contained subdirectory with its own `main.cpp` (when
applicable), `CMakeLists.txt`, `pdeapp.txt`, `pdemodel.txt`, and
mesh — there are no cross-section dependencies.

## Authoring vocabulary

| Word | Meaning |
|---|---|
| **Generated** | The model header is emitted by `text2code` from a `pdemodel.txt` description. The user writes math, the symbolic toolchain writes C++. |
| **Handwritten** | The user writes the model header (`my_model.hpp`) directly as a C++ struct that satisfies the `<exasim/model.hpp>` contract. No DSL, no autodiff, no codegen. |

## Driving vocabulary

| Word | Meaning |
|---|---|
| **Prebuilt** | The runtime is one of the binaries Exasim ships under `WITH_TEXT2CODE=ON` (`cput2cEXASIM`, `gput2cEXASIM`, `cpumpit2cEXASIM`, `gpumpit2cEXASIM`). The user writes only text files (`pdemodel.txt`, `pdeapp.txt`) and points the binary at them. |
| **CLI** | The user writes a tiny `main.cpp` that calls `exasim::run<MyModel>(argc, argv)`. Everything else (mesh, boundaries, solver options) comes from `pdeapp.txt`. |
| **Embedded** | The user constructs an `ExasimSolver<MyModel>` and configures it from C++. The simplest variant (`load_pdeapp`) still parses `pdeapp.txt`, but the solver lives inside the user's program rather than as a standalone binary. |
| **Distributed** | An embedded solver where each MPI rank builds its own slice of the mesh and hands it to `set_mesh_distributed`. ParMETIS repartitions inside `solve()`. |

## The six sections

| Section | Authoring | Driving | What the user writes |
|---|---|---|---|
| [01](../01-generated-prebuilt/README.md) | Generated | Prebuilt | `pdemodel.txt`, `pdeapp.txt` |
| [02](../02-generated-cli/README.md)      | Generated | CLI | `pdemodel.txt`, `pdeapp.txt`, 3-line `main.cpp`, `CMakeLists.txt` |
| [03](../03-generated-embedded/README.md) | Generated | Embedded | `pdemodel.txt`, `pdeapp.txt`, ~30-line `main.cpp`, `CMakeLists.txt` |
| [04](../04-handwritten-cli/README.md)    | Handwritten | CLI | `my_model.hpp`, `pdeapp.txt`, 3-line `main.cpp`, `CMakeLists.txt` |
| [05](../05-handwritten-embedded/README.md) | Handwritten | Embedded (in-memory) | `my_model.hpp`, ~80-line `main.cpp`, `CMakeLists.txt` |
| [06](../06-handwritten-distributed/README.md) | Handwritten | Distributed | `my_model.hpp`, ~140-line `main.cpp`, `CMakeLists.txt` |

There are eight authoring × driving combinations in principle.
The tutorial does not include "Handwritten + Prebuilt" because the
prebuilt binaries dispatch through the `libpdemodel.{so,dylib}`
ABI that `text2code` emits — there is no way to use them with a
hand-coded `my_model.hpp`. It also does not include separate
"Generated + In-memory" or "Generated + Distributed" sections,
because the embedded API is identical between generated and
handwritten models (you can swap `Poisson2D` for `GeneratedModel`
in section 05 or section 06 with no other change).

## Choosing a section

- **You want to run a PDE without writing any C++** → section 01.
- **You want to embed the solver in a larger C++ program** →
  section 03 (generated) or section 05 / 06 (handwritten).
- **You want full control over the math, including hand-derived
  Jacobians** → sections 04, 05, or 06.
- **You want to run on more than one rank with mesh data assembled
  in memory** → section 06.
- **You want to learn how the pieces fit together** → read all six
  in order.

## Backend variants

Each section's instructions show the build and run commands for
all four supported backend variants:

| Backend | Build directory | Notes |
|---|---|---|
| CPU | `build` | Single-rank, CPU only. |
| GPU | `build_gpu` | Single-rank, NVIDIA (CUDA) or AMD (HIP). |
| MPI | `build_mpi` | Multi-rank, CPU only. |
| MPI+GPU | `build_mpi_gpu` | Multi-rank, GPU per rank. |

The tutorial CMake targets carry no variant suffix; the build
directory disambiguates which backend the binary uses. The section
01 binaries are an exception: their names encode the variant
(`cput2cEXASIM` vs `gput2cEXASIM` vs `cpumpit2cEXASIM` vs
`gpumpit2cEXASIM`) because they live in the same build dir.

## grid.bin

`grid.bin` is the binary mesh file referenced by the `meshfile =
"grid.bin"` line in each section's `pdeapp.txt`. The format is a
flat blob of the form `[size(p), size(t), p(:), t(:)]`, where `p`
is the node-position array (`np × nd`) and `t` is the
element-vertex incidence array (`ne × nve`).

The tutorial directories ship a prebuilt `grid.bin` for the unit
square, copied from `apps/poisson/poisson2d/grid.bin`. That file
was originally produced by the legacy MATLAB / Julia / Python
frontends — for example
[`examples/Poisson/poisson2d/squaregrid.m`](../../examples/Poisson/poisson2d/squaregrid.m)
calls `squaremesh(n, n, 1, elemtype)` to generate `(p, t)` and
then `writebin("grid.bin", [size(p), size(t), p(:)', t(:)'])` to
serialize it.

Sections 05 and 06 sidestep `grid.bin` entirely by building the
same mesh as flat C arrays in `main.cpp` and handing them to
`solver.set_mesh(...)` or `solver.set_mesh_distributed(...)`.
