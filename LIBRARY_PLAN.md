# Plan: Exasim as a Header-Only C++ Template Library

## Goal

Exasim ships as a tree of C++ headers. The user's own project does the work:

```cpp
// my_project/main.cpp
#include <exasim/discretization.hpp>
#include <exasim/solver.hpp>
#include <exasim/solution.hpp>
#include "my_model.hpp"   // user's PDE: pointwise math + hand-written Jacobians

int main(int argc, char** argv) {
    /* MPI/Kokkos init */
    MyModel model{...};
    exasim::CSolution<MyModel> sol(/*…*/, model);
    sol.solve();
    /* finalize */
}
```

Concretely:

- **No `libexasim_core.a`.** No static or shared library installed by Exasim.
  `find_package(Exasim)` gives an `INTERFACE` target that adds the include
  directory and propagates Kokkos/MPI/BLAS dependencies. Eigen-style.
- **`CDiscretization`, `CSolver`, `CSolution`, `CPreconditioner` are class
  templates parameterized on `Model`.** Internals call `M::flux(...)`,
  `M::flux_jac_uq(...)`, etc. directly — no `extern "C"` indirection.
- **Models are plain C++ structs** with `KOKKOS_INLINE_FUNCTION static`
  pointwise functions and **hand-written Jacobians** (no autodiff, no DSL).
- **The `libpdemodel.hpp` extern-C ABI goes away.** `text2code` keeps
  working: it now emits a `my_model.hpp` header (not a `.so`), which the
  user can use directly in place of writing one by hand.

## Why this is a real rewrite, not packaging

The current `backend/Main/main.cpp` is a unity build: it `#include`s
~30 `.cpp` files which provide free function bodies, out-of-class member
function definitions, and external-linkage globals. That works only
because main.cpp is the single TU. To put Exasim in *user* code, every
one of those bodies has to be either:

- moved inside a class body (auto-inline), or
- marked `inline`, or
- moved into a class template (auto-inline), or
- converted to an `inline` variable (C++17, for namespace-scope mutables).

And the FEM internals — the load-bearing classes — currently call into
`KokkosFlux`/`HdgFlux`/… extern symbols. Those calls have to become
`M::flux(...)` template calls so the user's struct drives the math.

Honest estimate: **15–25 working days** for v1 (Poisson 2D
end-to-end), more for full Navier–Stokes coverage.

## State carried forward from the previous (static-library) attempt

These commits remain valid; the new plan builds on them. They're
already on `library-port`:

| commit | what |
|---|---|
| 839506da | Phase 0 inventory + numerical baseline + `include/exasim/` shims |
| d516063b | `common.h` ODR-safe (BLAS sentinels → extern; utilities → inline; self-contained includes) |
| 1d1bb244 | `cpuimpl.h` free functions marked inline |
| 7501cecf | `kokkosimpl.h` free functions marked inline + MPI baseline refresh |
| 8e2ed9e7 | `Visualization/visualization.cpp` → `.hpp` (single header-only class) |
| e545952a | inventory order revision (docs only) |
| 920b90fb | `HAVE_*` macros + `EXASIM_COMM_*` extern globals in `common.h`/`common.cpp` |

`baseline/verify.sh` still gates every phase: serial bit-identical to
master, MPI matches the post-1.2b/2 baseline (QoI invariant).

## Phases

### HOT.1 — Headerize all of `backend/`

Goal: `backend/<subdir>/*.cpp` → `*.hpp`, every body ODR-safe across
multiple TUs. The user can `#include <exasim/*.hpp>` from their own TU
and link nothing.

**Substeps**:

1. **Common**: convert `backend/Common/common.cpp` (BLAS sentinels +
   `EXASIM_COMM_*`) to C++17 `inline` variables in `common.h`; delete
   `common.cpp`. Drop the entry in `install/CMakeLists.txt`.

2. **Visualization**: already done.

3. **PointLocator / Discretization / Solver / Preconditioning /
   Solution**: per subdir:
   - rename every `.cpp` → `.hpp`;
   - update `#include "..*.cpp"` chains to `.hpp`;
   - mark every free function definition `inline`;
   - for out-of-class member function definitions
     (`void CDiscretization::method(...)`), choose per case:
     a) move inside the class body when small,
     b) mark `inline` at the out-of-class definition when large.
   - convert any external-linkage namespace-scope globals to `inline`.

4. **Preprocessing**: same shape as the others. text2code-related
   binary fixtures stay where they are; only the source files change.

5. **Model/**: handled in HOT.2.

After HOT.1, `main.cpp`'s `#include` chain works unchanged but every
file in it is a header. The user can `#include <exasim/run.hpp>` and
get the same orchestration code into their own TU. The
`libpdemodel.hpp` extern-C ABI is still in place — for now, the user
provides those symbols (or text2code does) just like today.

**Gate**: `baseline/verify.sh` after each subdir flip.

### HOT.2 — Template on `Model`, drop the libpdemodel ABI

The headline change. Replace the extern-C boundary with a C++ type
boundary.

**Substeps**:

1. **Define the model contract** in `<exasim/model.hpp>`. Required:
   - `static constexpr int nd, ncu, ncw, nparam`;
   - `KOKKOS_INLINE_FUNCTION static void flux(double f[], …)`;
   - companion `flux_jac_uq`, `flux_jac_w` for HDG (hand-written);
   - same shape for `source`, `fbou`, `ubou`, `fhat`, `init*`,
     `tdfunc`, `vis_*`, `qoi_*`.
   - `exasim::ModelDefaults<Self>` CRTP base supplying no-op defaults
     for the optional methods; user's struct only overrides what its
     PDE actually has.

2. **Per-kernel templates** in `<exasim/kernels/*.hpp>` — one header
   per kernel family (~12: flux, source, fbou, ubou, fhat, uhat,
   stab, init, tdfunc, eos, vis, qoi). Each contains a templated
   `kernel<M>(...)` that performs the SoA gather/scatter and calls
   `M::flux(...)` (or `M::flux_jac_uq(...)` for HDG) pointwise inside
   `Kokkos::parallel_for`. No autodiff anywhere.

3. **Template `CDiscretization`/`CPreconditioner`/`CSolver`/`CSolution`
   on `Model`**. Every call site that currently invokes
   `KokkosFlux(...)` / `HdgFbou(...)` / etc. (~30 sites across the
   FEM internals) becomes `exasim::*_kernel<M>(...)`. Member function
   bodies live inside the class template body so they auto-inline.

4. **Drop `libpdemodel.hpp`**. The model ABI no longer exists in the
   header-only world. text2code is updated in HOT.3.

5. **`<exasim/run.hpp>`** is now a small templated façade
   `exasim::run<Model>(int argc, char** argv)` for users who want the
   default driver; users with custom orchestration call into
   `exasim::CSolution<Model>` directly.

**Gate**: `baseline/verify.sh` after each step that completes a kernel
family. Math is preserved by hand-coded Jacobians matching the
SymEngine-generated ones.

### HOT.3 — User-facing example

`apps/library_example/poisson2d/`:

```
my_model.hpp       # struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
                   #     static constexpr int nd = 2, ncu = 1, ncw = 0, nparam = 1;
                   #     KOKKOS_INLINE_FUNCTION static void flux(...) {...}
                   #     KOKKOS_INLINE_FUNCTION static void flux_jac_uq(...) {...}
                   #     KOKKOS_INLINE_FUNCTION static void source(...) {...}
                   #     KOKKOS_INLINE_FUNCTION static void ubou(...) {...}
                   # };
main.cpp           # MPI/Kokkos init + exasim::run<Poisson2D>(argc, argv);
CMakeLists.txt     # find_package(Exasim REQUIRED)
                   # find_package(Kokkos REQUIRED)
                   # find_package(MPI REQUIRED)
                   # add_executable(poisson_solver main.cpp)
                   # target_link_libraries(poisson_solver PRIVATE
                   #     Exasim::headers Kokkos::kokkos MPI::MPI_CXX)
pdeapp.txt         # runtime config: porder, GMRES tol, BCs (unchanged from today)
```

**Gate**: numerical match against `baseline/poisson2d_serial/` and
`baseline/poisson2d_mpi2/`.

### HOT.4 — text2code retargeting

text2code currently emits `Code2Cpp.cpp` + per-kernel `Kokkos*1.cpp`
+ links them into `libpdemodelserial.so`. After HOT.2 the `.so` has
nowhere to link to. Retarget text2code:

- emit `my_model.hpp` (one header) with the user's `Model` struct.
- the per-kernel files collapse to pointwise method bodies inside that
  struct: SymEngine's symbolic differentiation produces the
  `flux_jac_uq` etc. bodies that the user would otherwise write by
  hand. Same code shape, different file layout, no `.so`.
- drop the `executeCppCode` / `buildDynamicLibraries` steps in
  `text2code/text2code/CodeCompiler.cpp` — the user's CMake build
  now compiles the model header into their own TU.

After HOT.4, both authoring paths converge:

- hand-written model: user writes `my_model.hpp` themselves.
- text2code-generated model: user writes `pdemodel.txt`, runs
  `text2code pdemodel.txt`, gets `my_model.hpp` to `#include`.

In both cases the build is the same: user's project compiles
`my_model.hpp` along with the `<exasim/...>` templates, no `.so`.

### HOT.5 — Documentation + install rules

- `cmake/ExasimConfig.cmake.in` exporting `Exasim::headers` (INTERFACE
  target with include dir + Kokkos/MPI/BLAS find_dependency calls).
- `install(DIRECTORY include/exasim DESTINATION …)`.
- `doc/header_only_api.md`: model contract, kernel template shape,
  the ~30-method surface that `Model` must implement (or inherit
  no-ops for from `ModelDefaults`).

## Risks and how we'll address them

| Risk | Mitigation |
|---|---|
| Compile time blow-up for the user's TU (~10k lines of Exasim headers parsed per TU) | Start with PCH support; later, separate compilation of large pieces via "extern templates" if needed |
| Numerical drift when replacing SymEngine-generated Jacobians with hand-coded ones | Keep text2code as one source of generated Jacobians; baseline gate after every kernel family |
| MPI globals become `inline` variables — does ParMETIS / link order still work? | Phase 0.2 baseline catches it. C++17 inline variables are a standardized feature for exactly this case |
| Member function templates explode binary size | Acceptable for a header-only library. If a user complains, ship explicit instantiation TUs as a compile-time toggle |
| `Discretization/ioutilities.cpp`'s `writearray`/`NumberToString` templates and `connectivity.cpp`'s `xiny2` are reachable from multiple subdirs via include order | After header-ification, all are `inline` template definitions in `.hpp` files; multiple inclusion is fine |

## Definition of done

- `apps/library_example/poisson2d/` builds against installed Exasim
  with `find_package(Exasim)`, runs serial and MPI, matches numerical
  baseline.
- No `.a` produced or installed.
- No `libpdemodel.hpp` referenced.
- `text2code` runs end-to-end and emits a header-only model that
  drops into the same workflow.

## What restarts now

The first concrete step of HOT.1 is to convert
`backend/Common/common.cpp` (the only `.cpp` we ourselves added in the
static-library track) into C++17 `inline` variables in `common.h`,
removing the file. After that, work proceeds subdir by subdir through
HOT.1.

## Status update — landed through HOT.2 step 3.3

Completed (all green against `baseline/verify.sh`):

- HOT.1: `backend/` headerized (every `.cpp` → `.hpp`, all bodies
  marked `inline`).
- HOT.2 step 1: `<exasim/model.hpp>` Model contract + ModelDefaults<Self>.
- HOT.2 step 2: 29 kernel templates in `<exasim/kernels/*.hpp>`,
  cross-example contract fixes (separate `fbou_hdg`, `nco`, naming).
- HOT.2 step 3.1: `<exasim/drivers.hpp>` — templated `*Driver<M>` wrappers.
- HOT.2 step 3.2: `<exasim/detail/abi_adapter.hpp>` — AbiAdapter marker.
- HOT.2 step 3.3: **the four FEM classes are templated on M** with default
  `M = AbiAdapter`. CDiscretization<M>, CSolver<M>, CPreconditioner<M>,
  CSolution<M> all exist. main.cpp uses `CSolution<>` (default arg)
  and the existing build is bit-identical to baseline.

## HOT.2 step 3.4 — remaining work, scoped precisely

The four FEM classes are templated, but `M` is currently a **phantom**
type parameter inside their member functions: bodies still call the
old non-templated `FluxDriver(...)`, `SourceDriver(...)`, etc.
defined in `backend/Model/{KokkosDrivers,ModelDrivers}.cpp`. So
`CSolution<MyModel>` *compiles* but doesn't actually invoke the
user's pointwise math — it goes through the libpdemodel.hpp ABI.

To make M actually flow through to the kernels, the FEM-internal
**chain free functions** also need templating, because they're the
ones that call the driver functions from inside a context where M
must be selected. Audit (after step 3.3 commit `e8c261a4`):

| File                                                | inline funcs | driver calls |
|-----------------------------------------------------|--------------|--------------|
| `Discretization/discretization.hpp`                 | (member fns) | 7            |
| `Discretization/postdiscretization.hpp`             | (member fns) | 7            |
| `Discretization/residual.hpp`                       | 13           | 5            |
| `Discretization/uequation.hpp`                      | 7            | 15           |
| `Discretization/uresidual.hpp`                      | 8            | 12           |
| `Discretization/qoicalculation.hpp`                 | 4            | 2            |
| `Discretization/getuhat.hpp`                        | 5            | 2            |
| `Discretization/qresidual.hpp`                      | 8            | 0            |
| `Discretization/wequation.hpp`                      | 3            | 0 (calls Hdg*) |
| `Discretization/qequation.hpp`                      | 7            | 0 (calls Hdg*) |
| `Solution/solution.hpp`                             | (member fns) | 3            |
| `Solution/postsolution.hpp`                         | (member fns) | 3            |
| **Total**                                           | **55**       | **56**       |

The 56 driver call sites are the leaf operations. The 55 free
functions (and the ~70 CDiscretization/CSolution member functions
already templated in 3.3) form the path between
`CSolution<M>::SolveProblem()` and a leaf driver call.

### Concrete plan for step 3.4

1. **Add a dispatch helper** `<exasim/detail/driver_dispatch.hpp>`
   defining `EXASIM_DRIVER_CALL(Name, ...)` — an `if constexpr` macro
   that routes to `::Name(...)` when `M = AbiAdapter` and to
   `::exasim::Name<M>(...)` otherwise. This was prototyped in this
   session and reverted; design is settled.

2. **Template the 55 chain free functions** on `M`:
   `inline void uEquationElemBlock(...)` →
   `template <class M> inline void uEquationElemBlock(...)`.
   Mechanical sed pass over the 8 chain files.

3. **Update the call sites of the chain functions** inside
   `CDiscretization<M>::method()` bodies to specify `<M>`:
   `uEquationElemBlock(args)` → `uEquationElemBlock<M>(args)`.
   Several hundred sites; mechanical sed pass.

4. **Replace the 56 leaf driver calls** with `EXASIM_DRIVER_CALL`
   invocations.

5. **Drop `libpdemodel.hpp`** in step 3.5 once nothing references it.

Estimated effort: 1–2 focused sessions (~3–6 hours). All edits are
mechanical. Risk: low (sed sweeps + compile error iteration), but
volume is high. Numerical baseline gates each step.

### What this session did not commit

A first attempt at step 3.4 — adding the dispatch header and macro,
sweeping driver call sites in 9 files — was reverted because the chain
free functions were not yet templated, leaving the macro expansion
referring to an undeclared `M`. The reverted edits are mechanically
correct in shape; they just need the chain templating to land first.
The branch is at the stable `e8c261a4` state.
