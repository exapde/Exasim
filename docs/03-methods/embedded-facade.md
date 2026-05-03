# Embedded — `ExasimSolver<M>`

`ExasimSolver<M>` is the C++ API for driving Exasim from a containing
program. The mesh, boundaries, polynomial order, and physics
parameters are set on the solver object; the converged solution is
returned through accessors. The mesh and parameters live in C++
values you pass to the solver, and the output stays in memory.

```cpp
#include <exasim/solver_facade.hpp>
#include "my_model.hpp"

int main() {
    Kokkos::initialize();
    {
        exasim::ExasimSolver<MyModel> solver;
        solver.set_mesh(p.data(), t.data(), nv, ne, nve);
        solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[1]) < 1e-8; });
        solver.set_polynomial_order(3);
        solver.set_physics_params({1.0});
        solver.solve();
        const double* udg = solver.udg();
    }
    Kokkos::finalize();
}
```

Model contract: [`model-contract.md`](model-contract.md).
Codegen authoring: [`codegen-text2code.md`](codegen-text2code.md).

## Sections

1. [Single-rank (CPU or GPU)](#1-single-rank-cpu-or-gpu)
2. [Multi-rank (MPI)](#2-multi-rank-mpi)
3. [Reading the converged state](#3-reading-the-converged-state)
4. [File-mesh fallback (`load_pdeapp`)](#4-file-mesh-fallback-load_pdeapp)
5. [The four backend variants](#5-the-four-backend-variants)
6. [Testing your model](#6-testing-your-model)
7. [Diagnosing test failures](#7-diagnosing-test-failures)

---

## 1. Single-rank (CPU or GPU)

```cpp
exasim::ExasimSolver<MyModel> solver;

// Mesh as flat C arrays.
//   p : nd   × np   doubles   (column-major: p[d + nd*j] = coord d of vertex j)
//   t : nve  × ne   ints      (column-major, 0-based vertex IDs per element)
solver.set_mesh(p.data(), t.data(),
                /*np=*/nv, /*ne=*/ne, /*nve=*/4);

// Boundary tagging. The predicate runs at every boundary face's
// centroid; `tag` indexes into the boundaryConditions array (added
// below). Multiple add_boundary calls with the same tag union.
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[1])     < 1e-8; });   // y=0
solver.add_boundary(/*tag=*/1,
    [](const double* x){ return std::abs(x[1] - 1) < 1e-8; });   // y=1

// Discretization.
solver.set_polynomial_order(3);    // p
solver.set_quadrature_order(6);    // 2p in 1D, picked per-element-type
solver.set_physics_params({1.0});  // arbitrary doubles → MyModel sees them as `mu`

// Run. Returns 0 on success (Newton converged within tolerance);
// throws on hard failure.
solver.solve();   // single-rank → mpiprocs=1 implied

// Solution accessors return host-side pointers. For GPU backends
// the first call after solve() copies device buffers to a member-
// cached host vector.
const double* udg = solver.udg();
int64_t       n   = solver.udg_size();
```

Curved boundaries:

```cpp
solver.add_curved_boundary(/*tag=*/2, [](const double* x){
    return x[0]*x[0] + x[1]*x[1] - 1.0;   // unit circle level-set
});
```

The runtime projects mesh nodes onto the curve via Newton step
`x ← x - f·∇f / |∇f|²`.

## 2. Multi-rank (MPI)

```cpp
MPI_Init(&argc, &argv);
EXASIM_COMM_WORLD = MPI_COMM_WORLD;
EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
int mpiprocs, mpirank;
MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
```

The `EXASIM_COMM_*` globals are the communicator Exasim uses
internally. Assign `MPI_COMM_WORLD` or any sub-communicator.

Each rank gives the solver its own slice of the global mesh, with
`t_local` referring to vertices by global node IDs. ParMETIS
repartitions for load balance internally; the initial split can be
any reasonable contiguous-range layout:

```cpp
int np_local, ne_local, node_off, elem_off;
compute_local_range(nv, mpiprocs, mpirank, np_local, node_off);
compute_local_range(ne, mpiprocs, mpirank, ne_local, elem_off);

std::vector<double> p_local(2 * np_local);
std::vector<int>    t_local(4 * ne_local);
// ... fill p_local with this rank's contiguous range of nodes
// ... fill t_local with this rank's contiguous range of elements
//     (the t_local entries are global vertex IDs; ParMETIS will
//     fetch the corresponding nodes from whatever rank holds them
//     during partition migration)

solver.set_mesh_distributed(
    p_local.data(), t_local.data(),
    /*np_local=*/np_local, /*ne_local=*/ne_local, /*nve=*/4,
    /*np_global=*/nv, /*ne_global=*/ne);

solver.add_boundary(...);
// add_boundary, set_polynomial_order, etc. are identical to single-rank

solver.solve(mpiprocs, mpirank);
```

The facade calls `CPreprocessing::takeParallel(comm)` internally:
ParMETIS partitions, DMD is built, ghost rows are sendrecv'd, and
the four runtime structs (app / master / mesh / sol) are built in
memory.

`EXASIM_FACADE_INMEMORY_MPI=0` env var routes through the file ABI
(per-rank `datain/{app,master,mesh,sol}.bin`) for debugging.

`add_boundary` predicates are evaluated on every face including
faces partitioned to other ranks; a predicate on rank 0 tags faces
that ParMETIS later assigns to any rank.

## 3. Reading the converged state

After `solver.solve()`:

```cpp
const double* udg     = solver.udg();        // (npe × nc × ne_local)
int64_t       udg_n   = solver.udg_size();

const double* uhat    = solver.uhat();       // (ncu × npf × nf)
int64_t       uhat_n  = solver.uhat_size();

const double* wdg     = solver.wdg();        // (npe × ncw × ne_local) if ncw>0
int64_t       wdg_n   = solver.wdg_size();
```

Layout: `udg` is column-major in `[npe, nc, ne_local]` order, so
`udg[i + npe*(c + nc*e)]` is field `c` at node `i` of local element
`e`. For HDG, `nc = ncu * (1 + nd)` — first `ncu` slots are `u`,
remaining are gradients `q`.

Multi-rank: each rank's `udg` is its locally-owned elements. For
post-processing, iterate per-rank and aggregate by global element ID.
The runtime exposes that mapping via `mesh.elempart[0..ne_local)`,
which the test harness uses (see `outelemid_np<r>.bin` sidecar in
§6).

## 4. File-mesh fallback (`load_pdeapp`)

If you have an existing `pdeapp.txt` workflow and want to drive it
through the facade without rewriting:

```cpp
exasim::ExasimSolver<MyModel> solver;
solver.load_pdeapp(argv[1], mpirank);  // parses pdeapp.txt
                                       //   + initializePDE
                                       //   + TextParser::parseFile
solver.solve(mpiprocs, mpirank);
```

`load_pdeapp` parses `pdeapp.txt` and runs the same in-memory
pipeline as the programmatic path. `apps/library_example/main_facade.cpp`
uses this; the same generic `main` is reused by all 12 codegen
examples.

## 5. The four backend variants

Backend is a compile-time decision. Same source file, four binaries:

| Variant | Kokkos lib | EXASIM flags | Binary suffix |
|---|---|---|---|
| CPU | `kokkos/buildserial` | `EXASIM_NOMPI=ON` | (none) |
| GPU | `kokkos/buildcuda` | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON` | `_gpu` |
| MPI | `kokkos/buildserial` | `EXASIM_MPI=ON` | `_mpi` |
| MPI+GPU | `kokkos/buildcuda` | `EXASIM_MPI=ON, EXASIM_CUDA=ON` | `_mpi_gpu` |

Configure lines per platform: [`cpu-gpu-mpi-mpigpu.md`](cpu-gpu-mpi-mpigpu.md)
and [`../01-installation.md`](../01-installation.md).

Multi-rank GPU device binding:

```cpp
#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
```

`shmrank` is the rank within the host-shared-memory communicator
(`MPI_COMM_TYPE_SHARED`).

## 6. Testing your model

Exasim's test harness compares each run's output to a recorded
baseline using a **DG element-L2 relative norm** that is invariant
to MPI partitioning and tolerant to platform numerics drift (Apple
Accelerate vs Intel MKL):

```
sqrt( sum_e ||u_run - u_baseline||² / sum_e ||u_baseline||² )
```

passes if relative L2 < 1e-3.

The runtime writes two files per rank into `dataout/`:

- `outudg_np<r>.bin` — solution values, `[npe, nc, ne_local]` doubles
- `outelemid_np<r>.bin` — `ne_local` int64s, the **global element
  ID** of each local element in `outudg`'s order

The validate harness (`apps/library_example/validate_codegen.sh`)
aggregates across ranks, sorts by global ID, and computes the
element-L2 norm against the baseline (a serial `np=1` recording).
Different `np`, different ParMETIS partition, different platform —
all wash out as long as the math is the same.

```bash
# Run a single example
bash apps/library_example/validate_codegen.sh poisson2d
bash apps/library_example/validate_codegen.sh --variant mpi --np 4 poisson2d

# Quick subset (Poisson2D + periodic + naca0012steady, ~30s on Mac)
bash apps/library_example/validate_codegen.sh --quick

# Full ctest matrix in one variant
ctest --test-dir build_cpu -j 8
```

To record a new baseline for your own model:

```bash
# Run your binary serially, point dataout at a fresh dir
./my_solver

# Copy the relevant files
mkdir -p baseline/my_app_serial
cp dataout/outudg_np0.bin    baseline/my_app_serial/
cp dataout/outelemid_np0.bin baseline/my_app_serial/

# Optional: outqoi.txt if your model defines QoIvolume / QoIboundary
cp dataout/outqoi.txt        baseline/my_app_serial/
```

Now any subsequent run (any backend, any partition) gates against
that baseline via `element_l2_diff.py`. See
`apps/library_example/element_l2_diff.py` for the comparison logic.

The full test infrastructure is documented in
[`../04-internals/testing.md`](../04-internals/testing.md), and the
baseline format in [`../04-internals/baselines.md`](../04-internals/baselines.md).

## 7. Diagnosing test failures

Failure messages contain the relative element-L2 magnitude:

```
[FAIL] foo — element-L2: outudg: relative element-L2 = 3.258e-01 (>1e-03) over 64 elements
```

| Relative L2 | Reading |
|---|---|
| 1e-12 to 1e-7 | FP-order drift; passes |
| 1e-3 to 1e-2 | borderline; check `NewtonTol`, `GMREStol` |
| > 1e-1 | real divergence; check partition / boundary handling |

Identical L2 magnitude on both `codegen:X:V` and `facade:X:V`
means the divergence is shared between paths.

Open known-failures: [`../04-internals/known-divergences.md`](../04-internals/known-divergences.md).

## See also

- Model contract: [`model-contract.md`](model-contract.md)
- Codegen authoring: [`codegen-text2code.md`](codegen-text2code.md)
- Backend builds: [`cpu-gpu-mpi-mpigpu.md`](cpu-gpu-mpi-mpigpu.md)
- Test harness: [`../04-internals/testing.md`](../04-internals/testing.md)
- Architecture: [`../04-internals/architecture.md`](../04-internals/architecture.md)
