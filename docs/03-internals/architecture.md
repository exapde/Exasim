# Architecture

This file explains the internal data flow from a user's `Model`
struct + mesh + parameters to a converged solution.

## Layers

```
                        ┌── User code ──┐
                        │ Model struct  │   (math: Flux, Source, Ubou, Initu, …)
                        │ pdeapp.txt    │   (config: porder, tau, NewtonTol, …)
                        │ grid.bin /    │   (mesh: vertices + connectivity)
                        │ flat C arrays │
                        └───────┬───────┘
                                ↓
              ┌─── Driver — pick one ───┐
              │ <exasim/run.hpp>        │   ← legacy CLI driver
              │ ExasimSolver<M>          │   ← embedded library API
              └─────────────┬───────────┘
                            ↓
       ┌────────────── Preprocessing ──────────────┐
       │ Mesh ingest (initializeMesh /            │
       │   meshFromArrays / meshFromArraysDistributed) │
       │ Master element (basis, quadrature)        │
       │ ParMETIS partitioning (MPI only)          │
       │ DMD (distributed mesh decomposition)      │
       │ Connectivity (face/element neighbors)     │
       │ Periodic node merging (MPI only)          │
       └─────────────────┬─────────────────────────┘
                         ↓
       ┌─────────────── Runtime ───────────────────┐
       │ CDiscretization<M>: residual, Jacobian    │
       │ CPreconditioner<M>: ILU / RAS / RB        │
       │ CSolver<M>: Newton + GMRES                │
       │ DIRK time stepping (if tdep)              │
       │ SaveSolutions / SaveQoI                   │
       └───────────────────────────────────────────┘
```

The user code at the top is everything *outside* of `backend/`. The
preprocessing layer is `backend/Preprocessing/`. The runtime is
`backend/Discretization/` + `backend/Solution/`.

## Authoring paths and where they meet

The two authoring paths produce the same `Model` struct shape:

```
Hand-written:                   Codegen:
my_model.hpp (you write)        pdemodel.txt (you write)
                                       │
                                       ↓ text2code
                                my_model.hpp (auto-generated)
                                       │
                                       ↓
                          struct MyModel : ModelDefaults<MyModel> {
                              KOKKOS_INLINE_FUNCTION static
                              auto Flux(…)   { … }
                              auto Source(…) { … }
                              auto Ubou(…)   { … }
                              …
                          };
```

The runtime instantiates `CSolution<MyModel>`. The full Model
contract (every method, every signature, default semantics, and
indexing layouts) is documented in
[`02-model-contract.md`](../02-model-contract.md). For a worked
example see
[`tutorial/04-handwritten-cli/`](../../tutorial/04-handwritten-cli/README.md).

## Driving paths

### Legacy CLI (`<exasim/run.hpp>`)

```
int main(int argc, char** argv) {
    return exasim::run<MyModel>(argc, argv);
}
```

`exasim::run` reads `argv[1]` as `pdeapp.txt`, parses it, runs
preprocessing (writing per-rank `datain/*.bin` in MPI mode),
constructs `CSolution<M>(filein, fileout, ...)` which reads the
bins back, solves, and writes `dataout/*.bin`. Self-contained
binary; configured by text file.

### Embedded library (`ExasimSolver<M>`)

```
exasim::ExasimSolver<MyModel> solver;
solver.set_mesh(p.data(), t.data(), nv, ne, nve);
solver.add_boundary(tag, [](const double* x){ … });
solver.set_polynomial_order(3);
solver.set_physics_params({1.0});
solver.solve();
const double* udg = solver.udg();
```

The mesh, boundary classifiers, polynomial order, and physics
parameters all come from in-memory C++ values. The solve runs the
same preprocessing pipeline (plus ParMETIS in MPI mode), but stays
in memory throughout — no `datain/*.bin` round-trip.

Both drivers ultimately call `CSolution<M>::SolveProblem` after
preprocessing. The runtime kernels are identical.

## Preprocessing pipeline

### Single-rank (CPU or GPU)

```
SerialPreprocessing()
├─ initializeMesh                  Mesh ← grid.bin / set_mesh arrays
├─ initializeMaster                Master ← (porder, pgauss, elemtype)
└─ writeBinaryFiles                writes datain/{app,master,mesh,sol}.bin
                                   (legacy path) OR
                                   builds Preprocessed bundle in memory
                                   (facade path via CPreprocessing::take())

CSolution<M>(filein, …) or         reads datain/*.bin or consumes
CSolution<M>(struct, …)            in-memory Preprocessed bundle
   ↓
postInit: compGeometry,            face geometry, mass matrix,
          compMassInverse,         residual scratch, qpoints
          allocResidual…
```

### Multi-rank (MPI)

```
ParallelPreprocessing(comm) or     same dispatch with MPI side-effects
takeParallel(comm)
├─ initializeParMesh / meshFromArraysDistributed
│                                  Mesh ← per-rank slice (global IDs)
├─ initializeMaster
├─ callParMetis(mesh, pde, comm)   ParMETIS partitions for load balance,
│                                  populates dmd.elempart_local
├─ initializeDMD                   builds DMD (distributed mesh decomposition):
│                                  elempart, elempartpts, elemsend / elemrecv,
│                                  nbsd, sendrecvpts, t2t with global IDs
│                                  (cross-rank neighbors via nbinfo)
│  ├─ mke2e_fill_first_neighbors
│  ├─ setboundaryfaces             marks t2t with -k for boundary k
│  ├─ compute_dgnodes               builds xdg from p
│  ├─ mergePeriodicNodeIDs          ★ unifies global node IDs across periodic faces
│  └─ setperiodicfaces              updates t2t with periodic neighbors
└─ writemesh + writesol            (legacy path) writes per-rank
                                   datain/{mesh,sol}<r+1>.bin OR
                                   in-memory build via takeParallel()

CSolution<M>(…) → SolveProblem
```

★ The `mergePeriodicNodeIDs` step is the fix for the periodic-MPI
divergence: serial `setperiodicfaces` rewrites mesh connectivity `t`
so periodic-paired vertices share IDs (so `mkf2e_hash` later finds
the periodic faces as interior faces). The parallel version in
`parmetisexasim.hpp` has to do the equivalent across rank
boundaries: it Allgathers boundary node coords, q-coord-matches
them via `xiny()`, and rewrites `mesh.nodeGlobalID` to canonical IDs
(smaller of each pair). Without this, the runtime's `buildConn`
treats periodic faces as physical boundaries.

## Runtime pipeline

`CSolution<M>::SolveProblem` dispatches to:

- `SteadyProblem` for `tdep == 0`
- `DIRK` for `tdep == 1` (multi-stage time stepping)
- `SteadyProblem_PTC` for runmode 10/11 (pseudo-transient continuation)

Each calls `NewtonSolver` → `hdgAssembleResidual` /
`hdgAssembleLinearSystem` → GMRES + preconditioner → `UpdateSolution`.

### Time-dependent + DIRK

```
for (istep = 0 .. tsteps):
    PreviousSolutions(sol, sys, common)        # capture u_prev
    for (j = 0 .. tstages):
        common.currentstage = j
        common.time = time + dt[istep] * DIRKcoeff_t[j]
        UpdateSource(sol, sys, app, res, common)
        SteadyProblem(out, backend)
        UpdateSolution<M>(sol, sys, app, res, tmp, common, backend)
                                                # 7-arg overload —
                                                # uses res.Rq as accumulator,
                                                # leaves sys.u untouched.
                                                # The 4-arg overload's
                                                # UpdateSolutionDIRK extracts
                                                # into sys.u and produces
                                                # divergent trajectories on
                                                # naca0012unsteady — don't use it.
    SaveSolutions / SaveQoI / SaveParaview / …
    time += dt[istep]
```

## Key data structures

| Struct | Lives in | Holds |
|---|---|---|
| `PDE` | `backend/Preprocessing/structs.hpp` | runtime config (porder, tau, NewtonTol, …) parsed from `pdeapp.txt` or set via facade |
| `Mesh` | same | vertices `p`, connectivity `t`, boundary expressions, periodic mappings |
| `Master` | same | basis functions, quadrature points, shape derivatives at element nodes |
| `DMD` | same | distributed mesh decomposition: `elempart`, `nbinfo`, `elemsend`, `elemrecv` |
| `commonstruct` / `appstruct` / `meshstruct` / `solstruct` | `backend/Common/common.h` | runtime structs (raw pointers, GPU-friendly); built from the preprocessing structs by `CPreprocessing::take()` |

The conversion from preprocessing structs (`vector<int>`, `vector<double>`)
to runtime structs (`int*`, `dstype*`) happens in `buildAppStruct`,
`buildMasterStruct`, `buildMeshStruct`, `buildSolStruct` (file
`backend/Preprocessing/buildstructs.hpp`).

## In-memory MPI path

By default the embedded `ExasimSolver<M>::solve(mpiprocs, mpirank)`
runs `CPreprocessing::takeParallel(comm)` — same logic as
`ParallelPreprocessing` but it builds the runtime structs directly
in memory rather than via per-rank `datain/*.bin` files. The
`CSolution<M>(Preprocessed&&, …)` constructor consumes the bundle
and forwards to `CDiscretization<M>(P&&, …)`.

The escape hatch `EXASIM_FACADE_INMEMORY_MPI=0` reverts to the
file ABI for debugging — it writes per-rank
`datain/{app,master,mesh,sol}.bin` via `ParallelPreprocessing`,
then constructs `CSolution<M>(filein, …)` which reads them back.
Same final answers (validated up to 1e-12 relative L2).

## Per-target text2code output

`text2code` emits `my_model.hpp` and compiles
`libpdemodel{serial,cuda,hip}.so` directly into the consumer's
codegen example directory (`apps/library_example/<name>_codegen/`),
not a global `backend/Model/` shared singleton. This decouples
parallel ctest gates: different examples never race on the same
generated library, and the four backend ctests can run independently.

The `--out-dir` flag on `text2code` controls where the generated
files land; `apps/library_example/regenerate.sh` passes the example's
own directory.

## Where to look in the source

| Topic | File |
|---|---|
| Embedded library API | `include/exasim/solver_facade.hpp` |
| Legacy CLI driver | `include/exasim/run.hpp` |
| Single-rank preprocessing | `backend/Preprocessing/preprocessing.hpp` |
| MPI preprocessing | `backend/Preprocessing/preprocessing.hpp` (takeParallel) + `backend/Preprocessing/parmetisexasim.hpp` |
| Runtime kernels | `backend/Discretization/discretization.hpp` |
| Solver | `backend/Solution/solution.hpp` |
| Output | `backend/Solution/postsolution.hpp` |
| text2code | `text2code/text2code/` |
| In-memory struct builders | `backend/Preprocessing/buildstructs.hpp` |
