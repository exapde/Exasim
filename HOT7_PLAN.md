# HOT.7 — In-memory programmatic API: drop pdeapp.txt, datain/, dataout/

## Why

HOT.6.7 added a programmatic `CPreprocessing` constructor that takes
pre-populated `PDE` / `InputParams` / `ParsedSpec` structs instead of
parsing `pdeapp.txt`. `apps/library_example/solve_square` exercises
it: one `main()`, no `pdeapp.txt`, no `pdemodel.txt`, hand-written
model. It works.

But the program still:

1. Writes `grid.txt` to disk so `initializeMesh()` can read it back.
2. Has `CPreprocessing` write `datain/{app,master,mesh,sol}.bin` so
   `CSolution<M>` can read it back.
3. Has `CSolution<M>` write `dataout/outudg_np<r>.bin` for the user
   to read post-solve, instead of returning the result struct.
4. Forces boundary classification through string expressions parsed
   by tinyexpr at runtime, instead of taking C++ predicates.

Each of these is a file-IO hop or a string-parsing detour where the
data is already in C++. HOT.7 removes them.

## Goal

By end of HOT.7, an external app can:

```cpp
ExasimSolver<Poisson2D> solver;
solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4, /*nd=*/2);
solver.set_boundary([](const double* x){ return std::abs(x[1])<1e-8; }, /*tag=*/1);
solver.set_physics_params({1.0});
solver.set_polynomial_order(3);
solver.solve();
auto& udg = solver.solution_udg();   // host- or device-side view
```

— no files written or read in the hot path.

## Sequencing

Done in slices, each landable independently. Old paths (`pdeapp.txt`
+ `text2code` + `datain/`) keep passing the existing validate sweep
after every slice.

### Slice 7.0 (done in HOT.6.7)
- `CPreprocessing(PDE, InputParams, ParsedSpec, rank, comm)` — programmatic constructor
- `pdeFinalizeDerived(PDE&)` — exported helper, runs the
  flag/problem/factor/solversparam packing and the hybrid/tdep
  derivations
- `apps/library_example/solve_square/` — uses both
- Old paths verified byte-identical (poisson2d_codegen,
  poisson2d_template)

### Slice 7.1 — programmatic mesh
**Goal:** replace `initializeMesh(params, pde)` reading from
`pde.meshfile` on disk with an overload that accepts the mesh data
in memory.

- New header `backend/Preprocessing/inmemory_mesh.hpp` exposing:
  ```cpp
  inline Mesh meshFromArrays(const double* p, const int* t,
                             int np, int ne, int nve, int nd,
                             const InputParams& params, const PDE& pde);
  ```
- Move the disk-reading half of `initializeMesh` into a thin wrapper
  over `meshFromArrays`.
- `solve_square` switches from `write_square_mesh()` + file → array
  build + direct call.

**Old path:** `initializeMesh(params, pde)` keeps reading
`grid.bin/.txt/.msh` exactly as it does today. CI gates unchanged.

### Slice 7.2 — boundary expressions as predicates
**Goal:** boundary tagging today goes through tinyexpr-parsed strings
(`"abs(y)<1e-8"`). Add a parallel API that takes
`std::function<bool(const double*)>`.

- `setboundaryfaces` gains an overload:
  ```cpp
  int setboundaryfaces(int* f, int* t2lf, int* localfaces,
                       const double* p, const int* t,
                       const std::vector<std::function<bool(const double*)>>& exprs,
                       int nd, int elemtype, int ne, int nbndexpr);
  ```
- `Mesh` struct gains an optional
  `std::vector<std::function<bool(const double*)>> boundaryPreds`.
- When `params.boundaryPreds.size() > 0`, the predicate path is
  used; when empty, the legacy string path runs (unchanged).

**Old path:** strings still work via tinyexpr (no functional
regression). Codegen examples + hand-written `poisson2d_template`
keep passing.

### Slice 7.3 — in-memory ABI between CPreprocessing and CSolution<M>
**Goal:** today CPreprocessing serializes app/master/mesh/sol to
binary files, then CSolution<M> deserializes them. Replace this with
a typed handoff struct.

- New header `<exasim/preprocessed.hpp>` defining
  ```cpp
  struct Preprocessed {
      appstruct app;
      masterstruct master;
      meshstruct mesh;
      solstruct sol;
  };
  ```
- `CPreprocessing` exposes `Preprocessed take()` that builds the
  struct in memory (no `writepde`/`writemaster`/`writemesh`/`writesol`
  calls).
- `CSolution<M>` gains a constructor:
  ```cpp
  CSolution(Preprocessed&& pre, std::string fileout, ...);
  ```
  that consumes `Preprocessed` directly instead of reading from
  `filein/`.
- `pde.writemeshsol = 0` skips the file writes in the old path; the
  file writes become opt-in for users who want to inspect / restart.

**Old path:** legacy main flow keeps `writemeshsol=1`; the existing
file-driven `CSolution(filein, ...)` constructor still reads from
disk and reproduces today's behavior bit-for-bit.

### Slice 7.4 — in-memory output
**Goal:** today `CSolution<M>::SolveProblem` writes outudg/outuhat as
side effects. Add accessors that return the host views instead.

- `CSolution<M>::host_udg()`, `host_uhat()`, `host_qoi_volume()`,
  `host_qoi_boundary()` returning `std::span<const dstype>` on
  host memory (with a device→host copy if running on GPU/HIP).
- `pde.saveSolOpt = 0` makes the file writes opt-out; they remain
  the default for backward compatibility.

**Old path:** legacy bins are still written; `validate_codegen.sh`
still diffs them.

### Slice 7.5 — the `ExasimSolver<M>` façade
**Goal:** wrap 7.0–7.4 in the high-level fluent API shown in the
"Goal" section above. Lives in `<exasim/solver_facade.hpp>`.

### Slice 7.6 — cross-arch validation
- Run all four backends (CPU / GPU / MPI / MPI+GPU) against the new
  `solve_square` with each slice.
- Add CI workflows for GPU + MPI variants once self-hosted runners
  are available; CPU-only CI runs `solve_square` after slice 7.0.

### Slice 7.7 — tutorial rewrite
- Move §5 of `doc/tutorial_cpu_gpu_mpi.md` to §1; demote the codegen
  / `pdeapp.txt` flow to "convenience: codegen path".
- Show `ExasimSolver<M>` in the §3 "Run a model" block; the legacy
  `mpirun ./poisson2d_codegen pdeapp.txt` invocation moves to an
  appendix.

## Backward-compatibility invariant

After every slice, the following must pass byte-identical:

```
bash apps/library_example/validate_codegen.sh                      # CPU
bash apps/library_example/validate_codegen.sh --variant gpu        # GPU
bash apps/library_example/validate_codegen.sh --variant mpi --np 2 # MPI
bash apps/library_example/validate_codegen.sh --variant mpi_gpu --np 2  # MPI+GPU
```

The `_codegen`, `_codegen_gpu`, `_codegen_mpi`, `_codegen_mpi_gpu`
binaries from HOT.6.4 read `pdeapp.txt` and write `datain/`/`dataout/`
as they do today. Their CMake targets, link lines, and command-line
contracts are frozen.

## Effort estimate

- 7.1 mesh-from-arrays: ~0.5d
- 7.2 boundary predicates: ~0.5d
- 7.3 in-memory ABI: ~1d (largest — touches CSolution + CDiscretization)
- 7.4 in-memory output: ~0.5d
- 7.5 façade: ~0.5d
- 7.6 cross-arch: ~0.5d
- 7.7 tutorial: ~0.25d

Total: ~4 days focused; landable in 4–5 sessions.

## Open questions

1. Boundary predicates on GPU: lambdas captured by-value run as
   host-side callbacks during preprocessing only (not in the solve
   loop), so GPU portability is not a constraint. Confirmed.
2. ParMETIS partitioning: still needs the global element list as a
   contiguous array. `meshFromArrays` is fine for that — no in-memory
   API needed in the partitioning step itself.
3. CUDA-aware MPI dependency: stays as the multi-rank GPU
   prerequisite. Documented in the tutorial.
