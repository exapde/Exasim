# Codegen — `pdemodel.txt` + text2code

Author the math symbolically in the SymEngine DSL. text2code reads
`pdemodel.txt` and emits:

- `my_model.hpp` — a templated `Model` struct (`GeneratedModel`)
- `libpdemodel{serial,cuda,hip}.{so,dylib}` — the same kernels
  pre-compiled for the legacy AbiAdapter path

The generated `my_model.hpp` satisfies the same contract as a
hand-written one ([`model-contract.md`](model-contract.md)) and is
consumed by `<exasim/run.hpp>` and `ExasimSolver<M>` identically.

## Authoring

`pdemodel.txt` declares the equation in a small DSL:

```
outputs Flux, Source, Tdfunc, Ubou, Fbou, FbouHdg, Initu, VisScalars

function Flux(x, uq, v, w, eta, mu, t)
    output_size(f) = 2;
    f[0] = mu[0] * uq[1];
    f[1] = mu[0] * uq[2];
end

function Source(x, uq, v, w, eta, mu, t)
    output_size(s) = 1;
    s[0] = 2.0 * pi * pi * sin(pi * x[0]) * sin(pi * x[1]);
end

…
```

text2code generates the Jacobians via SymEngine's symbolic
differentiation; you don't write them by hand.

`pdeapp.txt` declares the runtime config (porder, tau, NewtonTol,
boundary conditions, mesh file, etc.). Layout:

```
my_app/
├── pdemodel.txt       # SymEngine DSL → math
├── pdeapp.txt         # runtime config
├── grid.bin           # mesh
└── main.cpp           # 3 lines
```

`main.cpp`:

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"
int main(int argc, char** argv) {
    return exasim::run<GeneratedModel>(argc, argv);
}
```

## Generation

`apps/library_example/regenerate.sh <name>` runs text2code:

```bash
cd apps/library_example/poisson2d_codegen
bash ../regenerate.sh poisson2d
ls
# my_model.hpp                    ← include it from main.cpp
# libpdemodelserial.dylib         ← linked at runtime
# Code2Cpp.cpp + SymbolicFunctions.{cpp,hpp}  ← intermediate
# datain/                          ← serialized mesh
```

text2code writes everything into the example's own directory; nothing
goes to a shared `backend/Model/`.

## Run

```bash
cmake --build build_cpu --target poisson2d_codegen
cd apps/library_example/poisson2d_codegen
$EXASIM/build_cpu/poisson2d_codegen ./pdeapp.txt
```

For GPU / MPI / MPI+GPU substitute the matching binary suffix and
build directory. See [`cpu-gpu-mpi-mpigpu.md`](cpu-gpu-mpi-mpigpu.md).

## With the embedded API

The same generated `my_model.hpp` works with `ExasimSolver<M>`. The
generic main lives at `apps/library_example/main_facade.cpp`:

```cpp
exasim::ExasimSolver<GeneratedModel> solver;
solver.load_pdeapp(argv[1], mpirank);
solver.solve(mpiprocs, mpirank);
```

Build target: `<name>_facade`. Each codegen example reuses the same
generic `main_facade.cpp` with the example's `my_model.hpp` on the
include path.

## In-tree examples

12 codegen apps under `apps/library_example/<name>_codegen/`:

```
poisson2d  poisson3d  periodic  naca0012steady  naca0012unsteady
lshape  isoq3d  isoq  cone  orion  nsmach8  sharpb2
```

Each gets eight binaries (4 backends × `_codegen` and `_facade`).

## See also

- Hand-written authoring: [`hand-written-model.md`](hand-written-model.md)
- Embedded driving: [`embedded-facade.md`](embedded-facade.md)
- Legacy AbiAdapter driving: [`abi-adapter.md`](abi-adapter.md)
- Tutorial: [`../../tutorial/`](../../tutorial/README.md)
