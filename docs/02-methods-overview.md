# Methods overview

Three independent choices when running Exasim.

| Dimension | Options |
|---|---|
| Authoring | hand-written `Model` struct, or symbolic `pdemodel.txt` + text2code |
| Driving | legacy CLI `exasim::run<M>(argc, argv)`, or embedded library `ExasimSolver<M>`, or AbiAdapter (text2code-generated `libpdemodel*.so` loaded via `cput2cEXASIM` and friends) |
| Backend | CPU, GPU (CUDA / HIP), MPI, MPI+GPU |

Every authoring × driving × backend combination works. The model
contract, runtime kernels, and pdeapp/grid contract are identical
across them.

## Authoring

### Hand-written `Model`

Write a C++ struct with `KOKKOS_INLINE_FUNCTION static` methods
(Flux, Source, Ubou, Initu, Tau, …). Full contract in
[`03-methods/model-contract.md`](03-methods/model-contract.md).

### Codegen via `pdemodel.txt`

Write the math symbolically in the SymEngine DSL. `text2code` reads
`pdemodel.txt`, generates `my_model.hpp`, and compiles
`libpdemodel{serial,cuda,hip}.{so,dylib}`. The generated `my_model.hpp`
satisfies the same contract a hand-written struct does.

Walkthrough: [`03-methods/codegen-text2code.md`](03-methods/codegen-text2code.md).

## Driving

### Legacy CLI — `exasim::run<M>(argc, argv)`

`main.cpp` is three lines:
```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"
int main(int argc, char** argv) {
    return exasim::run<MyModel>(argc, argv);
}
```

The binary reads `pdeapp.txt` from disk, runs preprocessing into
`datain/*.bin`, solves, writes `dataout/*.bin`.

Walkthrough: [`03-methods/hand-written-model.md`](03-methods/hand-written-model.md).

### Embedded library — `ExasimSolver<M>`

```cpp
exasim::ExasimSolver<MyModel> solver;
solver.set_mesh(p.data(), t.data(), nv, ne, nve);
solver.add_boundary(tag, [](const double* x){ … });
solver.set_polynomial_order(3);
solver.solve();
const double* udg = solver.udg();
```

The mesh, boundaries, polynomial order, and physics parameters
are set in C++. The converged solution is returned through
accessors on the solver object. The driver keeps its data in
memory across the entire solve.

Walkthrough: [`03-methods/embedded-facade.md`](03-methods/embedded-facade.md).

### AbiAdapter — `cput2cEXASIM` / `gpumpit2cEXASIM` and friends

The original Exasim driving path. `backend/Main/main.cpp` calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`; the AbiAdapter
dispatches every kernel call through `dlsym` lookup against
`libpdemodel{serial,cuda,hip}.{so,dylib}`. text2code generates the
library; the binary loads it at runtime.

Walkthrough: [`03-methods/abi-adapter.md`](03-methods/abi-adapter.md).

## Backend

Compile-time choice driven by CMake flags. Same source code, four
binaries:

| Variant | Kokkos build | EXASIM flags | Binary suffix |
|---|---|---|---|
| CPU | `kokkos/buildserial` | `EXASIM_NOMPI=ON` | (none) |
| GPU | `kokkos/buildcuda` or `buildhip` | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON` (or `EXASIM_HIP=ON`) | `_gpu` |
| MPI | `kokkos/buildserial` | `EXASIM_MPI=ON` | `_mpi` |
| MPI+GPU | `kokkos/buildcuda` or `buildhip` | `EXASIM_MPI=ON, EXASIM_CUDA=ON` | `_mpi_gpu` |

Build instructions per platform in
[`01-installation.md`](01-installation.md). Backend-specific runtime
notes in [`03-methods/cpu-gpu-mpi-mpigpu.md`](03-methods/cpu-gpu-mpi-mpigpu.md).

## See also

For a guided tour through every authoring × driving combination on
the same Poisson problem, see
[`tutorial/`](../tutorial/README.md).
