# 00 — Overview

To solve a PDE with Exasim the user must specify three things:

| Axis | What it is |
|---|---|
| **PDE math** | The equations: flux, source, boundary residual, initial condition, and (for HDG) the trace residual `FbouHdg`. |
| **Solver setup** | Discretization order, quadrature order, Newton/GMRES tolerances, physics parameters, boundary tagging. |
| **Mesh** | Node coordinates `p` and element-vertex incidence `t` for the spatial discretization. |

For each of the three, the user picks one of two forms:

| Axis | "built" form (external artifact) | "C++" form (in your program) |
|---|---|---|
| **PDE math** | `pdemodel.txt` (a SymEngine DSL; `text2code` emits the C++ header) | `my_model.hpp` (a struct against `<exasim/model.hpp>`) |
| **Solver setup** | `pdeapp.txt` (a key-value text file) | `ExasimSolver<M>` API calls |
| **Mesh** | `grid.bin` (a binary blob; produced by a small generator that the tutorial ships) | flat C arrays passed to `set_mesh` / `set_mesh_distributed` |

Six sections walk through the lattice, ordered by how much of the
specification has been pulled from external artifacts into C++:

| Section | PDE | Solver | Mesh | Notes |
|---|---|---|---|---|
| [01](../01-generated-prebuilt/README.md)    | built | built | built | The user writes no C++ at all; runs a prebuilt Exasim binary (`cput2cEXASIM` etc.) on the two text files. |
| [02](../02-generated-cli/README.md)         | built | built | built | The user writes a 3-line `main.cpp` that calls `exasim::run<GeneratedModel>(argc, argv)` and owns the resulting binary. |
| [03](../03-generated-embedded/README.md)    | built | built | built | The solver lives as a C++ object you own (`ExasimSolver<GeneratedModel>`), but its parameters still come from `pdeapp.txt` via `load_pdeapp`. |
| [04](../04-handwritten-cli/README.md)       | C++   | built | built | Same shape as section 02, but the model is a hand-written `Poisson2D` struct instead of a generated one. |
| [05](../05-handwritten-embedded/README.md)  | C++   | C++   | C++   | Single rank. The mesh, boundary predicates, polynomial order, and physics parameters are all set in C++; the converged solution is read back through `solver.udg()`. |
| [06](../06-handwritten-distributed/README.md) | C++ | C++   | C++   | Distributed across MPI ranks. Each rank builds its own slice of the mesh and hands it to `set_mesh_distributed`; ParMETIS repartitions inside `solve()`. |

There are eight points in the 2×2×2 lattice in principle. Sections
01 and 02 share the same axes and differ only in who owns the
binary. Sections 05 and 06 share the same axes and differ only in
whether the mesh is global or per-rank. The tutorial does not
include "handwritten + prebuilt" because the prebuilt binaries
dispatch through the `libpdemodel.{so,dylib}` ABI that `text2code`
emits, and that ABI is incompatible with a hand-coded
`my_model.hpp`.

## Choosing a section

- **You want to run a PDE without writing any C++** → section 01.
- **You want to embed the solver in a larger C++ program** →
  section 03 (generated PDE) or section 05 / 06 (handwritten PDE).
- **You want full control over the math, including hand-derived
  Jacobians** → sections 04, 05, or 06.
- **You want a true distributed-memory solve where each rank only
  ever materializes its share of the mesh** → section 06.
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
directory disambiguates which backend the binary uses. The
section 01 binaries are an exception: their names encode the
variant (`cput2cEXASIM` vs `gput2cEXASIM` vs `cpumpit2cEXASIM` vs
`gpumpit2cEXASIM`) because they live in the same build dir.

## grid.bin

`grid.bin` is the binary mesh file referenced by `meshfile =
"grid.bin"` in the `pdeapp.txt` of sections 01–04. The format is
`[size(p), size(t), p(:), t(:)]` written column-major as 64-bit
floats — `nd, np, nve, ne` followed by `np × nd` node coordinates
and `ne × nve` 1-based vertex IDs.

The tutorial does not check in any `grid.bin`. The script
[`tutorial/tools/squaregrid.py`](../tools/squaregrid.py)
regenerates it on demand: it imports `SquareMesh` from
`frontends/Python/Mesh/squaremesh.py` (pure NumPy; no `gmsh`
dependency for the unit square) and writes the legacy
`writebin`-compatible binary blob. `tutorial/run-all.sh` invokes
the script before running each section that needs it.

To regenerate by hand:

```bash
python3 $EXASIM/tutorial/tools/squaregrid.py 16 \
    $EXASIM/tutorial/01-generated-prebuilt/grid.bin
```

The first argument is the per-side element count `n`; the script
emits an `n × n` Cartesian quad mesh on the unit square.

Sections 05 and 06 sidestep `grid.bin` entirely by building the
same mesh as flat C arrays in `main.cpp` and handing them to
`solver.set_mesh(...)` or `solver.set_mesh_distributed(...)`.
