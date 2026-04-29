# Phase 0.1 Inventory — backend/ as it stands

Snapshot of the runtime translation unit before any restructuring.
Used by Phase 1 to guide the split into a real static library.

## main.cpp's include chain

`backend/Main/main.cpp` is currently the **only TU compiled into every
executable**. It pulls everything else in via `#include "*.cpp"`.

### Public-shaped headers it consumes (already `.h` / `.hpp`)

```
backend/Common/common.h            # types, MPI/Kokkos macros, structs, memory utils
backend/Common/cpuimpl.h           # CPU-side numerical helpers
backend/Common/kokkosimpl.h        # Kokkos-side numerical helpers
backend/Common/pblas.h             # BLAS/MPI wrappers
```

### Implementation `.cpp` files included as if they were headers

```
backend/Discretization/discretization.cpp     # CDiscretization class
backend/Preconditioning/preconditioner.cpp    # CPreconditioner class
backend/Solver/solver.cpp                     # CSolver class
backend/Visualization/visualization.cpp       # CVisualization class
backend/PointLocator/pointlocator.cpp         # CPointLocator class
backend/Solution/solution.cpp                 # CSolution class
backend/Preprocessing/preprocessing.cpp       # CPreprocessing class (only when HAVE_SHARED_MODEL_LIB)
```

### Transitive `#include "*.cpp"` chains

| Top-level | Pulls in (relative to its directory) |
|---|---|
| `Discretization/discretization.cpp` | `discretization.h`, `gpuDeviceInfo.cpp`, `ioutilities.cpp`, `../Model/ModelDrivers.cpp`, `../Model/BuiltIn/BuiltinModelDrivers.cpp`, `../Model/KokkosDrivers.cpp`, `connectivity.cpp`, `readbinaryfiles.cpp`, `setstructs.cpp`, `residual.cpp`, `matvec.cpp`, `qoicalculation.cpp` |
| `Preconditioning/preconditioner.cpp` | `preconditioner.h`, `setprecondstruct.cpp`, `applymatrix.cpp` |
| `Solver/solver.cpp` | `solver.h`, `setsysstruct.cpp`, `getpoly.cpp`, `gmres.cpp`, `ptcsolver.cpp` |
| `Visualization/visualization.cpp` | (none — single file) |
| `PointLocator/pointlocator.cpp` | `pointlocator.h`, `pointlocation.cpp`, `pointwallmodel.cpp` |
| `Solution/solution.cpp` | `solution.h`, `previoussolutions.cpp`, `updatesolution.cpp`, `updatesource.cpp`, `timestepcoeff.cpp`, `avsolution.cpp` |
| `Preprocessing/preprocessing.cpp` | `preprocessing.h`, `TextParser.hpp`, `tinyexpr.cpp`, `helpersexasim.cpp`, `readpdeapp.cpp`, `readmesh.cpp`, `makemeshexasim.cpp`, `makemaster.cpp`, `makemasterexasim.cpp`, `domaindecomposition.cpp`, `writebinaryfilesexasim.cpp`, `parmetisexasim.cpp` |

## Subdir contents (full file listing)

```
backend/Common
  common.h  cpuimpl.h  kokkosimpl.h  pblas.h

backend/Discretization
  connectivity.cpp  discretization.cpp  discretization.h  geometry.cpp
  getuhat.cpp  gpuDeviceInfo.cpp  ioutilities.cpp  ismeshcurved.cpp
  KokkosDrivers.cpp  massinv.cpp  matvec.cpp  postdiscretization.cpp
  postdiscretization.h  qequation.cpp  qoicalculation.cpp  qresidual.cpp
  readbinaryfiles.cpp  residual.cpp  setstructs.cpp  uequation.cpp
  uresidual.cpp  wequation.cpp

backend/Preconditioning
  applymatrix.cpp  postpreconditioner.cpp  preconditioner.cpp
  preconditioner.h  setprecondstruct.cpp

backend/Solver
  getpoly.cpp  gmres.cpp  postsolver.cpp  ptcsolver.cpp
  setsysstruct.cpp  solver.cpp  solver.h

backend/Visualization
  visualization.cpp

backend/PointLocator
  pointcandidates.cpp  pointellipsoidgrid.cpp  pointinterpolation.cpp
  pointlocation.cpp  pointlocator.cpp  pointlocator.h
  pointwallmodel.cpp

backend/Solution
  avsolution.cpp  postsolution.cpp  postsolution.h  previoussolutions.cpp
  solution.cpp  solution.h  timestepcoeff.cpp  updatesolution.cpp
  updatesource.cpp

backend/Preprocessing
  connectivity.cpp  domaindecomposition.cpp  helpers.cpp  helpersexasim.cpp
  makemaster.cpp  makemasterexasim.cpp  makemesh.cpp  makemeshexasim.cpp
  parmetis.cpp  parmetisexasim.cpp  preprocessing.cpp  preprocessing.h
  readmesh.cpp  readpdeapp.cpp  structs.hpp  text2code.cpp  TextParser.hpp
  tinyexpr.cpp  tinyexpr.h  utils.cpp  wallmodelsampling*.cpp
  wallmodelsampling.h  writebinaryfiles.cpp  writebinaryfilesexasim.cpp
  (binary fixtures: gaussnodes.bin, grid.bin, masternodes.bin,
   partition.bin; text fixtures: pdeapp.txt, pdemodel.txt)
```

## ODR landmines (must fix before splitting into multiple TUs)

These work today because main.cpp is the only TU — splitting will surface them.

### External-linkage globals in `Common/common.h`

```
common.h:178:  dstype one = 1.0;
common.h:179:  dstype minusone = -1.0;
common.h:180:  dstype zero = 0.0;
common.h:193:  dstype cublasOne[1] = {one};
common.h:194:  dstype cublasMinusone[1] = {minusone};
common.h:195:  dstype cublasZero[1] = {zero};
```

**Fix**: declare `extern` in `common.h`, define once in a new `common.cpp`.

### Free-function bodies in `Common/{cpuimpl,kokkosimpl}.h`

- `cpuimpl.h`: ~52 non-template, non-static function definitions (1795 lines).
- `kokkosimpl.h`: ~122 non-template, non-static function definitions (3308 lines).
- `pblas.h`: 0 free-function bodies (all wrappers are `inline` / templated).

**Fix options** (per function):
1. Move definition to a new `.cpp` (preferred for non-templates).
2. Mark `inline` (small bodies).
3. Mark `static` (internal-linkage; only if the function is genuinely
   per-TU — unlikely here since the whole point is shared across modules).

For Phase 1 we expect most of these to move into `cpuimpl.cpp` and
`kokkosimpl.cpp` — **header keeps declarations, .cpp gets bodies**.

### `Model/ModelDrivers.cpp`, `Model/KokkosDrivers.cpp`

Currently `#include`d transitively via `Discretization/discretization.cpp`.
These dispatch to the user's PDE model (`KokkosFlux`, `HdgFlux`, …
declared in `Model/libpdemodel.hpp`). They are themselves implementation
files that need to become `.cpp` TUs of `exasim_core` — and they call
into the model symbols, which the user's `glue.cpp` (or text2code's
`libpdemodelserial.so`) provides at link time.

## What becomes public vs. private

### Public (installed under `include/exasim/`)

```
common.h, cpuimpl.h, kokkosimpl.h, pblas.h            (declarations only)
discretization.hpp                                    (CDiscretization + post)
preconditioner.hpp
solver.hpp
solution.hpp                                          (CSolution + post)
visualization.hpp                                     (DEFERRED to Phase 1.2 — class is in .cpp today)
pointlocator.hpp
preprocessing.hpp                                     (already library-shaped)
libpdemodel.hpp                                       (the model ABI)
run.hpp                                               (NEW — top-level driver)
model.hpp                                             (NEW — Phase 3, model contract docs)
kernels/*.hpp                                         (NEW — Phase 3, templated kernels)
```

#### Phase 1.1 status (this commit)

The shim files exist and forward to `backend/<dir>/*.h`. They let a future
TU `#include <exasim/foo.hpp>` and pick up the same declarations as today.

`visualization.hpp` is **not** shipped yet: `backend/Visualization/` has only
a `.cpp` file with the class body inline, so a forwarding shim would either
(a) include a `.cpp` and re-introduce the unity-build / ODR problem we are
trying to fix, or (b) need a brand-new `visualization.h` carved out — that
work belongs in Phase 1.2.

### Private (compiled into `exasim_core`, not installed)

Everything else under `backend/{Discretization,Solver,Solution,
Preconditioning,Visualization,PointLocator,Preprocessing}/*.cpp` — plus
the new `common.cpp`, `cpuimpl.cpp`, `kokkosimpl.cpp` carved out of the
ODR-landmine fix above.

## Build matrix to preserve

The following compile definitions currently flow from
`install/CMakeLists.txt` to `main.cpp`'s preprocessor:

```
_MPI         → HAVE_MPI
_CUDA        → HAVE_GPU + HAVE_CUDA
_HIP         → HAVE_GPU + HAVE_HIP
_TEXT2CODE   → HAVE_TEXT2CODE  → HAVE_SHARED_MODEL_LIB
_BUILTINMODEL→ HAVE_BUILTINMODEL → HAVE_SHARED_MODEL_LIB
_OPENMP      → HAVE_OPENMP (via -fopenmp)
_ENZYME      → HAVE_ENZYME
_MUTATIONPP  → HAVE_MPP
```

After the split, these become **PUBLIC interface compile definitions**
on `exasim_core_<backend>` so that any TU including `<exasim/*.h>` sees
the same `HAVE_*` switches. Expected library set:

- `exasim_core_serial` — no `_MPI` / `_CUDA` / `_HIP` (pure CPU, single rank)
- `exasim_core_serial_mpi` — `_MPI`
- `exasim_core_cuda`, `exasim_core_cuda_mpi` — `_CUDA` (+ `_MPI`)
- `exasim_core_hip`,  `exasim_core_hip_mpi`  — `_HIP`  (+ `_MPI`)

Each picks up the matching pre-built Kokkos under `kokkos/build{serial,cuda,hip}/`.

## Phase 1 execution order

1. **Common/** — fix ODR landmines, produce `common.cpp` /
   `cpuimpl.cpp` / `kokkosimpl.cpp`, install headers under
   `include/exasim/`. Lowest dependency, surfaces most of the build risk.
2. **Visualization/** — single file, no chain. Easy second.
3. **PointLocator/** — three .cpp files, no cross-deps to other backend modules.
4. **Discretization/** — biggest chain (12 files + Model drivers). The
   `Model/{Kokkos,Model}Drivers.cpp` move with it.
5. **Preconditioning/**, **Solver/**, **Solution/** — depend on Discretization.
6. **Preprocessing/** — depends on text2code-shared headers (`TextParser.hpp`,
   `structs.hpp`); compile-gated by `_TEXT2CODE` / `_BUILTINMODEL`.
7. **`backend/Main/main.cpp`** body moves to `<exasim/run.hpp>` +
   `src/run.cpp`; main becomes a 3-liner.

After each step: rebuild `cput2cEXASIM` / `cpumpit2cEXASIM`, run the
Phase 0.2 baseline, diff numerical outputs.
