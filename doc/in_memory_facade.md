# In-memory Exasim — `ExasimSolver<M>` walkthrough

An end-to-end guide to running Exasim *without* `pdeapp.txt` /
`datain/` / `dataout/` round-trips on disk. The math, the mesh, the
parameters, and the solution all live in C++ values you pass to and
read from a façade object.

This is the path most users should take for new code:

```cpp
#include <exasim/solver_facade.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    Kokkos::initialize();
    {
        exasim::ExasimSolver<MyModel> solver;
        solver.set_mesh(p.data(), t.data(), nv, ne, nve);
        solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[1]) < 1e-8; });
        solver.set_polynomial_order(3);
        solver.set_physics_params({1.0});
        solver.solve();

        const double* udg = solver.udg();
        // ... do something with the solution
    }
    Kokkos::finalize();
}
```

The legacy `<exasim/run.hpp>` + `pdeapp.txt` + `cmd_line` driver stays
supported and is a fine choice when you want a CLI binary that ingests
text files. The facade is the right choice when Exasim is *embedded*
in another C++ program — design optimization loops, reduced-order
modeling, coupled multi-physics, anywhere you want to drive Exasim
from a containing program rather than from a shell.

For the underlying `Model` contract (the per-element pointwise math
the FEM internals call) see `doc/header_only_api.md`. For the codegen
authoring path (write `pdemodel.txt`, get `my_model.hpp` for free)
see `apps/library_example/README.md`. This doc is about the *driver*
side — how to feed mesh + state to the solver and pull results back
out.

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

Typed boundary predicates beat the legacy tinyexpr strings in two
ways. (a) compile-time checked — typos and wrong field accesses
fail to build, not at runtime parse time. (b) full C++ — capture
state, call any function, no string DSL.

Polynomial-defined boundary geometry uses
`add_curved_boundary(tag, level_set_fn)` — same shape, takes a
level-set function `f(x)` that's zero on the curve. The runtime does
a Newton step `x ← x - f·∇f / |∇f|²` to project mesh nodes onto the
curve.

## 2. Multi-rank (MPI)

```cpp
MPI_Init(&argc, &argv);
EXASIM_COMM_WORLD = MPI_COMM_WORLD;
EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
int mpiprocs, mpirank;
MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
```

You can call MPI either before or after `Kokkos::initialize()` — both
are tested. The two `EXASIM_COMM_*` globals are the communicator
Exasim uses internally; assigning `MPI_COMM_WORLD` (or any sub-
communicator you've made) is enough.

Each rank then gives the solver its own slice of the global mesh,
with `t_local` referring to vertices by **global** node IDs. The
partitioner (ParMETIS) repartitions for load balance, so you can
hand it any reasonable initial split — contiguous stripes work
fine:

```cpp
int np_local, ne_local, node_off, elem_off;
compute_local_range(nv, mpiprocs, mpirank, np_local, node_off);
compute_local_range(ne, mpiprocs, mpirank, ne_local, elem_off);

std::vector<double> p_local(2 * np_local);
std::vector<int>    t_local(4 * ne_local);
// ... fill p_local with this rank's contiguous range of nodes
// ... fill t_local with this rank's contiguous range of elements
//     (vertex IDs are GLOBAL — ParMETIS doesn't need the rank to
//     own all referenced vertices)

solver.set_mesh_distributed(
    p_local.data(), t_local.data(),
    /*np_local=*/np_local, /*ne_local=*/ne_local, /*nve=*/4,
    /*np_global=*/nv, /*ne_global=*/ne);

solver.add_boundary(...);
// add_boundary, set_polynomial_order, etc. are identical to single-rank

solver.solve(mpiprocs, mpirank);
```

Under the hood the facade calls `CPreprocessing::takeParallel(comm)`,
which runs ParMETIS, builds the DMD (distributed mesh decomposition),
sendrecvs ghost rows, and builds the four runtime structs (app /
master / mesh / sol) **all in memory**. No per-rank `datain/*.bin`
files are written or read. (HOT.7.13 made this the default. The
`EXASIM_FACADE_INMEMORY_MPI=0` env var falls back to the legacy
file ABI as an escape hatch.)

`add_boundary` predicates run on every face the solver inspects —
including faces owned by other ranks before partitioning. That's
how a single boundary predicate on rank 0 can tag faces that ParMETIS
later assigns to rank 5.

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

`load_pdeapp` reads the same `pdeapp.txt` the legacy `<exasim/run.hpp>`
binary would — `meshfile`, `boundaryconditions`, `boundaryexpressions`,
all of it. The facade then reads the mesh from disk and runs the
same in-memory pipeline as the programmatic path.

This is what `apps/library_example/main_facade.cpp` uses — one
generic main reusable across all 12 codegen examples.

## 5. The four backend variants

Backend is a **compile-time** decision, not a runtime one. The same
source file produces four binaries:

| variant         | Kokkos lib          | EXASIM flags                     | binary suffix    |
| --------------- | ------------------- | -------------------------------- | ---------------- |
| CPU             | `kokkos/buildserial` | `EXASIM_NOMPI=ON`               | (none)           |
| GPU             | `kokkos/buildcuda`  | `EXASIM_NOMPI=ON, EXASIM_CUDA=ON` | `_gpu`          |
| MPI             | `kokkos/buildserial` | `EXASIM_MPI=ON`                  | `_mpi`           |
| MPI+GPU         | `kokkos/buildcuda`  | `EXASIM_MPI=ON, EXASIM_CUDA=ON`  | `_mpi_gpu`       |

All four CMake builds are independent dirs that all consume the same
sources; switching is just a different `cmake -B build_<variant>`
configure plus the matching Kokkos `-DKokkos_DIR`. See
`doc/tutorial_cpu_gpu_mpi.md` for the exact configure lines and
mpirun usage.

For multi-rank GPU, each rank pins one device:

```cpp
#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
```

`shmrank` is the rank within the host-shared-memory communicator
(`MPI_COMM_TYPE_SHARED`), so multi-host runs with N GPUs per host
land cleanly.

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

The full test infrastructure is documented in `doc/testing.md`.

## 7. Diagnosing test failures

When a test fails, the harness prints a relative element-L2 magnitude:

```
[FAIL] foo — element-L2: outudg: relative element-L2 = 3.258e-01 (>1e-03) over 64 elements
```

Use the magnitude to localize the problem:

| relative L2     | likely cause                                              |
| --------------- | --------------------------------------------------------- |
| ~ 1e-12 to 1e-7 | passes; FP-order drift across platforms                   |
| ~ 1e-3 to 1e-2  | borderline; check `NewtonTol` and `GMREStol` in pdeapp    |
| > 1e-1          | real divergence — check partition / boundary handling     |

If both `codegen:X:V` and `facade:X:V` produce **identical** L2
magnitudes (lockstep), the divergence isn't facade-specific — it's
something both paths share (most likely a partition-handling bug
under MPI, or recorded baseline drift).

If only the `facade` stem fails, suspect a facade-specific
regression — start by comparing the per-iteration Newton trajectory
against `codegen` (capture stdout from both runs and `diff`).

The `doc/known_test_divergences.md` file tracks open known-failures
that aren't regressions — currently `orion:mpi`, `orion:mpi_gpu`,
`periodic:mpi`, `periodic:mpi_gpu` (all partition-handling bugs in
the legacy code, both stems lockstep).

---

## Where to go next

- **Authoring the math**: `doc/getting_started.md` (hand-written),
  `apps/library_example/README.md` (codegen via `pdemodel.txt`).
- **Model contract reference**: `doc/header_only_api.md`.
- **Backend matrix details**: `doc/tutorial_cpu_gpu_mpi.md`.
- **Test harness layers**: `doc/testing.md`.
- **App catalog (12 example PDEs)**: `doc/running_apps.md`.
