# Mixed-Precision Mode — Implementation Plan

This plan is grounded in `research/mixed-precision-codebase-survey.md`.
Read that document first; the rationale for the choices below lives
there.

The plan is split into phases that ship behavior-preserving by
default and let the user opt into mixed precision via a single
template parameter. Each phase is independently mergeable.

---

## 0. Summary

- Introduce a **precision policy** struct that names the precision of
  each numerically-distinct data layer (geometry, solution,
  residual/Jacobian, Krylov, preconditioner, halo).
- Thread the policy as a defaulted second template parameter on
  `CDiscretization`, `CSolver`, `CPreconditioner` — alongside the
  existing `M` parameter.
- Replace the global `typedef dstype` with policy-driven type
  aliases inside the templates. Provide a `dstype` back-compat alias
  for the default policy so existing code paths keep working.
- Provide `ArrayCast<TFrom,TTo>` Kokkos kernels and `mpi_type<T>()`
  traits to handle conversion boundaries.
- Templatize the BLAS dispatch (pblas.h) so the same call site can
  run S* or D* based on the *type of the pointer*, not a global
  `#ifdef USE_FLOAT`.
- Roll out by lowering one layer at a time, gated by the policy
  field, with regression on tutorial cases at each step.

---

## 1. Design

### 1.1 Precision policy struct

```cpp
namespace exasim {

// Empty marker for "leave at default precision".
struct DefaultPrecision {
    using TGeom    = double;   // shape funcs, jac, Xx, mesh coords
    using TState   = double;   // udg, wdg, uh, xdg, sdg, odg
    using TResid   = double;   // Ru, Rh, Rq
    using TJac     = double;   // H, K, F, Mass, Minv, C, D, B, G, E
    using TKrylov  = double;   // sys.v, sys.r, sys.b, sys.x, sys.u
    using THess    = double;   // s, cs, sn, y, H (Hessenberg), lam
    using TPrec    = double;   // precond.W, precond.U
    using THalo    = double;   // MPI halo buffers
    using TModel   = double;   // pointwise math inside flux/source/etc.
};

// Classic "FP32 preconditioner + FP64 Krylov": one knob change.
struct LowPrecPreconditioner : DefaultPrecision {
    using TPrec   = float;
    using TJac    = float;    // implies Minv, H, K stored in float
};

// Mixed-precision iterative refinement: lower everything inside the
// linear solve, keep outer Newton in FP64.
struct MixedIR : DefaultPrecision {
    using TResid  = float;
    using TJac    = float;
    using TKrylov = float;
    using TPrec   = float;
    using THalo   = float;
    using TModel  = float;
    // THess and TState stay FP64 — outer Newton update is FP64.
};

} // namespace exasim
```

The policy is a *trait struct*. Users can write their own policies by
inheriting and overriding individual `using`s.

### 1.2 Class-template surgery

Today:
```cpp
template <class M = exasim::detail::AbiAdapter> class CDiscretization;
template <class M = exasim::detail::AbiAdapter> class CSolver;
template <class M = exasim::detail::AbiAdapter> class CPreconditioner;
```

After:
```cpp
template <class M = exasim::detail::AbiAdapter,
          class Prec = exasim::DefaultPrecision>
class CDiscretization;
// same for CSolver, CPreconditioner.
```

Existing call sites that name only `M` keep working because `Prec`
defaults. Existing `main.cpp` (`backend/Main/main.cpp`) and the
library examples need no source change.

### 1.3 Storage struct surgery

Two options:

**Option A — templatize the structs:** `solstruct<Prec>`,
`resstruct<Prec>`, … and pick the right one inside the FEM class.
- Pro: cleanest separation; type system enforces correctness.
- Con: very wide change; every helper in
  `setstructs.hpp`/`ioutilities.hpp` must take a `Prec`.

**Option B — keep the structs with a single precision per struct
field, but selected at compile time from the policy:** the struct
becomes a thin wrapper that holds the right pointer types via
`using`:
```cpp
template <class Prec> struct solstruct_t {
    typename Prec::TState* udg = nullptr;
    typename Prec::TState* wdg = nullptr;
    typename Prec::TGeom*  xdg = nullptr;
    …
};
using solstruct = solstruct_t<DefaultPrecision>;   // back-compat alias
```
- Pro: surgical, behavior-preserving at default policy, no breakage
  in `main.cpp`.
- Con: still touches every struct, but only at the type-alias level.

**Recommendation: Option B.** Adopt template structs, but expose
`dstype`-typed aliases for the default policy so the legacy
ABI/codegen path keeps compiling. The chosen precision is then a
single template parameter threaded down from the top.

### 1.4 BLAS dispatch

Replace the `#ifdef USE_FLOAT` macro path in `pblas.h` with
trait-driven overloads:

```cpp
template <class T> struct blas_traits;
template <> struct blas_traits<float> {
    static constexpr auto nrm2  = SNRM2;
    static constexpr auto dot   = SDOT;
    static constexpr auto gemm  = SGEMM;
    static constexpr auto cublas_gemm = cublasSgemm;
    static constexpr auto mpi_type    = MPI_FLOAT;
    …
};
template <> struct blas_traits<double> { /* DGEMM, MPI_DOUBLE, … */ };

// Templated free function replaces the `#ifdef USE_FLOAT` body:
template <class T>
inline void PDOT(cublasHandle_t handle, Int m, T* x, Int incx, T* y,
                 Int incy, T* global_dot, Int backend) { … }
```

The existing macro names (`cublasGEMM`, etc.) can remain as aliases
to `blas_traits<dstype>::gemm` for the legacy build.

### 1.5 Cast kernels

```cpp
template <class TFrom, class TTo>
inline void ArrayCast(TTo* dst, const TFrom* src, int N) {
    Kokkos::parallel_for("ArrayCast", N, KOKKOS_LAMBDA(const size_t i) {
        dst[i] = static_cast<TTo>(src[i]);
    });
}
```
Lives next to `ArrayCopy` in `backend/Common/kokkosimpl.h`. Use these
at all boundaries identified in research §4.

### 1.6 Codegen / model boundary

- Short term: do **not** change codegen. The codegen ABI stays
  `dstype = TState`. All "low precision" work happens at the
  *FEM-internal* level, with casts at the MatVec boundary.
- Medium term: extend `text2code/CodeGenerator.cpp` to emit
  `template <class T>` versions of the generated kernels. The
  templated `M::flux` path (`kernels/flux.hpp`) only needs the kernel
  body to accept a `T*` instead of `dstype*` and use `T` for the
  stack-local arrays. Models opting into FP32 evaluation declare it
  on their `M` type:
  ```cpp
  struct MyModel { using compute_t = float; … };
  ```

### 1.7 Where conversions live

| Boundary                          | Cast direction           | New site                                                |
|-----------------------------------|--------------------------|----------------------------------------------------------|
| MatVec entry                      | `TKrylov → TJac/TModel`  | `evalMatVec` start (`discretization.hpp:798-813`)        |
| MatVec exit                       | `TJac → TKrylov`         | `evalMatVec` end                                         |
| ApplyPreconditioner entry         | `TKrylov → TPrec`        | `CPreconditioner::ApplyPreconditioner`                   |
| ApplyPreconditioner exit          | `TPrec → TKrylov`        | same                                                     |
| MPI halo pack                     | `TState/TJac → THalo`    | `GetCollumnAtIndex` → buffsend (`matvec.hpp:268-300`)    |
| MPI halo unpack                   | `THalo → TState/TJac`    | buffrecv → `PutCollumnAtIndex`                           |
| Dot product accumulator           | `T → double`             | `PDOT` reduction loop in pblas (`pblas.h:278-321`)       |
| Newton update                     | `TKrylov dx → TState u`  | `PTCsolver`/`NonlinearSolver` update step                |

---

## 2. Phased rollout

Each phase is independently testable and can be merged separately.

### Phase 0 — scaffolding (no behavior change)

**Deliverables:**
- `include/exasim/precision.hpp` defining `DefaultPrecision`,
  `LowPrecPreconditioner`, `MixedIR`, and the `blas_traits<T>`
  primary template + `float`/`double` specializations.
- `ArrayCast<TFrom,TTo>` kernel in `kokkosimpl.h`.
- `mpi_type<T>()` helper.
- Templatize `CDiscretization<M, Prec=DefaultPrecision>`,
  `CSolver<M, Prec=DefaultPrecision>`,
  `CPreconditioner<M, Prec=DefaultPrecision>`. Defaults preserve
  current behavior; *no FEM-internal code changes yet.*

**Verification:**
- Build the whole repo; the `dstype` typedef still wins and every
  test passes bit-identically. (Run the tutorial poisson/naca tests;
  they assert binary-equal residuals.)
- `tests/ctest` parity.

### Phase 1 — type-aliased structs (no behavior change)

**Deliverables:**
- Convert `solstruct`, `resstruct`, `sysstruct`, `precondstruct`,
  `masterstruct`, `meshstruct`, `tempstruct` to
  `*_t<Prec>` templates with `*_t<DefaultPrecision>` typedef'd back
  to the legacy name.
- `CDiscretization<M,Prec>::sol` becomes `solstruct_t<Prec>`,
  similarly for the rest.
- `setstructs.hpp`, `ioutilities.hpp`, `setsysstruct.hpp`,
  `setprecondstruct.hpp` get templated on `Prec`. The bodies don't
  change; only the field-pointer types reflect the policy.

**Verification:**
- Bit-identical residuals at `DefaultPrecision`.
- Instantiate `CDiscretization<AbiAdapter, DefaultPrecision>` and
  `CDiscretization<AbiAdapter, LowPrecPreconditioner>` both compile.
  At this phase the second one only *allocates* `res.K`/`res.Minv` in
  FP32 — no FP32 *kernel* sees them yet, so it crashes at first use.
  That's expected; phase 2 wires the kernels.

### Phase 2 — templated BLAS dispatch

**Deliverables:**
- Replace `#ifdef USE_FLOAT` macro bodies in `pblas.h` with templated
  function bodies via `blas_traits<T>`.
- Replace `MPI_DOUBLE` literals with `mpi_type<THalo>()` calls.
- Kokkos kernels in `kokkosimpl.h` become templates:
  `template <class T> inline void ArrayAXPBY(T* y, T* x, const T* z, T a, T b, int N)`.
  Each call site picks its own `T`. The Kokkos lambda body is
  unchanged (already type-generic).

**Verification:**
- Run with `Prec = DefaultPrecision` (everything `double`) and the
  global `USE_FLOAT` build (everything `float`). Both must produce
  the same answers as before phase 2.
- Add a unit test: a hand-rolled `PDOT<float>` and `PDOT<double>` on
  the same data, comparing against a reference reduction.

### Phase 3 — FP32 preconditioner (the first real mixed-precision config)

**Deliverables:**
- `LowPrecPreconditioner` policy uses `TPrec = TJac = float`,
  everything else `double`.
- In `CPreconditioner::ApplyPreconditioner`:
  - Cast input vector `x` from `TKrylov` to `TPrec` (using
    `ArrayCast`) into a temporary FP32 buffer.
  - Run the existing block-Jacobi / Schwarz / Block-ILU0 path on
    FP32 storage (`res.K`, `res.Minv`).
  - Cast the result back to `TKrylov`.
- Build the preconditioner data (`hdgBlockILU0`,
  `hdgBlockJacobi`, `compMassInverse`) directly in FP32, or build in
  FP64 and down-cast — pick whichever is cheaper to implement
  (probably down-cast first, optimize to native-build later).
- Polynomial preconditioner: compute `getPoly` (Ritz values, Leja
  sort) in FP64; apply `ApplyPoly` in `TPrec`.

**Verification:**
- Re-run `tests/`: GMRES iteration count should change by ≤ ~10% on
  the tutorial cases. Final solution should match FP64 to about FP64
  tolerance (because the *outer* solver is FP64 — the preconditioner
  is just an approximation).
- Profile: HDG `hdgBlockJacobi`/`Block-ILU0` apply time should
  roughly halve.

### Phase 4 — FP32 Krylov (TKrylov = float)

**Deliverables:**
- `sys.v` (Krylov basis), `sys.r`, `sys.b`, `sys.x`, `sys.u` become
  FP32 storage.
- `s`, `cs`, `sn`, `y`, `H` (Hessenberg), `lam`, the polynomial-
  preconditioner Ritz values stay FP64 (`THess`).
- `MGS`/`CGS` dot products accumulate in FP64 even on FP32 inputs.
  Implement either via custom Kokkos reductions or via the
  `cublasSdot`-with-`double`-result extension (cuBLAS supports mixed
  precision dot via `cublasSdotEx`).
- `evalMatVec` casts `v` and the residual workspace in/out at its
  boundary.
- Switch default orthogonalization to CGS with one reorthogonalization
  pass when `is_same_v<TKrylov, float>` — this trades a single extra
  PGEMTV per step for stability.

**Verification:**
- Iteration counts should rise modestly. Final solution error
  (computed against the FP64-everywhere baseline) should be within a
  few ULPs of `eps(TKrylov)` for well-conditioned cases.
- New regression case: a stiff Reynolds=10⁶ NACA run. Verify the
  Newton outer loop still converges — if not, that's the signal to
  bump `TKrylov` back up for that problem.

### Phase 5 — Mixed-precision Newton iterative refinement

**Deliverables:**
- `MixedIR` policy: `TState = THess = double`, everything else `float`.
- `PTCsolver`/`NonlinearSolver` outer loop:
  - Evaluate `b = R(u)` in `TResid = float`, then up-cast to FP64.
    *Alternative:* evaluate in FP64 (keep `TResid = double`) once per
    Newton step. The latter is safer; pick it first.
  - Build the FP32 preconditioner and FP32 Krylov state.
  - Run GMRES → produces `dx` in `TKrylov = float`.
  - Up-cast `dx` to FP64.
  - Apply `u_FP64 += dx_FP64` (the update is in FP64 to recover
    Newton's quadratic convergence).
- Document the divergence guard: if `nrmr > 1.05 * old_nrmr` (the
  existing check in `LinearSolver` at `ptcsolver.hpp:89`), the RB
  preconditioner is reset; for mixed-precision we add a secondary
  guard that increases `TKrylov` to FP64 for the next outer step.

**Verification:**
- Compare nonlinear iteration count (`PTCsolver` outer iterations)
  vs. the FP64-everywhere baseline on tutorial poisson/naca:
  expect ≤ 2× outer iterations.
- Wall-clock: expect ~1.5-2× speedup on memory-bandwidth-bound HDG
  cases. CPU-only runs likely gain less.

### Phase 6 — Model evaluation in low precision (optional)

**Deliverables:**
- Add `using compute_t = T;` to the Model concept (`include/exasim/model.hpp`).
- Make the kernel templates in `include/exasim/kernels/*.hpp`
  parameterize on `T = M::compute_t` (currently they hard-code
  `double` for the stack arrays). The Kokkos lambda body becomes
  `T x[nd]; T uq[Nq]; …`.
- Codegen update (`text2code/CodeGenerator.cpp`): emit
  `template <class T>` versions of each generated kernel, with
  `T* f`, `const T* x`, etc. Existing codegen string literals (line
  1509-2240) replace `dstype` with the type parameter.
- Selecting model precision is then per-Model, decoupled from the
  global `Prec` policy.

**Verification:**
- For each tutorial: run with `M::compute_t = float` and verify the
  difference vs `M::compute_t = double` is within FP32 machine
  precision on the residual.
- Flag problematic models where this fails (typically: cv compressible
  flow with pressure-from-energy recovery).

### Phase 7 — MPI halo in low precision (optional)

**Deliverables:**
- Halo packing/unpacking routes through `ArrayCast` between the live
  field precision (`TState`, `TJac`) and the halo buffer precision
  (`THalo`).
- MPI calls use `blas_traits<THalo>::mpi_type` instead of the
  hard-coded `MPI_DOUBLE` literals.

**Verification:**
- Run 4-rank tutorial cases (poisson2d, naca) and compare against
  the FP64-halo baseline. Difference should be within `eps(THalo) ×
  N_halo_messages` of the baseline.
- Test interface flux maps with `THalo = float`: these are the most
  sensitive halo (carry per-face coupling data), so they may need
  to stay in `TState` precision via a `THaloInterface` distinct
  policy field.

---

## 3. File-level work breakdown

The phases above touch these files (non-exhaustive; verify by
`grep -l dstype`):

**Phase 0 (new files):**
- `include/exasim/precision.hpp` (new)
- `include/exasim/detail/blas_traits.hpp` (new)
- `include/exasim/detail/mpi_traits.hpp` (new)

**Phase 0 (edits, defaulted-template):**
- `backend/Discretization/discretization.h`
- `backend/Solver/solver.h`
- `backend/Preconditioning/preconditioner.h`

**Phase 1 (template-ize the structs in `common.h`):**
- `backend/Common/common.h` — `solstruct`, `resstruct`, `sysstruct`,
  `precondstruct`, `masterstruct`, `meshstruct`, `tempstruct`,
  `appstruct`, `commonstruct`, `wallmodelstruct` (and the
  `freememory`/`printinfo`/`sizeoffloat` members of each)
- `backend/Discretization/setstructs.hpp`,
  `backend/Discretization/ioutilities.hpp`
- `backend/Solver/setsysstruct.hpp`
- `backend/Preconditioning/setprecondstruct.hpp`

**Phase 2 (BLAS/Kokkos templating):**
- `backend/Common/pblas.h` — convert every macro-dispatched function
  body to a function template.
- `backend/Common/kokkosimpl.h` — convert each `inline void Name(dstype*, …)`
  to `template <class T> inline void Name(T*, …)`.
- `backend/Common/cpuimpl.h` — same.

**Phase 3 (preconditioner):**
- `backend/Preconditioning/preconditioner.hpp`,
  `backend/Preconditioning/applymatrix.hpp`,
  `backend/Preconditioning/setprecondstruct.hpp`,
  `backend/Preconditioning/postpreconditioner.hpp`
- `backend/Discretization/massinv.hpp`

**Phase 4 (Krylov):**
- `backend/Solver/gmres.hpp`, `backend/Solver/getpoly.hpp`,
  `backend/Solver/setsysstruct.hpp`

**Phase 5 (mixed-IR Newton):**
- `backend/Solver/ptcsolver.hpp`
- `backend/Solver/solver.hpp`

**Phase 6 (model templating):**
- `include/exasim/drivers.hpp`, `include/exasim/kernels/*.hpp`
- `text2code/text2code/CodeGenerator.cpp`
- `backend/Model/libpdemodel.hpp` (extend ABI to add `_float`
  variants if needed)

**Phase 7 (halo):**
- `backend/Discretization/matvec.hpp`,
  `backend/Discretization/residual.hpp`,
  `backend/Discretization/qequation.hpp`,
  `backend/Discretization/wequation.hpp`,
  `backend/Solver/setsysstruct.hpp` — all `MPI_DOUBLE` sites.

---

## 4. Verification plan

A standing verification harness across all phases:

1. **Bit-identity test** (phases 0-2): run tutorial poisson2d at
   `Prec=DefaultPrecision`, assert per-iteration residual norm
   matches the pre-change build byte-for-byte.
2. **Approximate-equivalence test** (phases 3+): for each new policy,
   measure
   - Newton outer iteration count vs. baseline (warning > 50%
     increase)
   - GMRES inner iteration count vs. baseline (warning > 50%
     increase)
   - Final solution `‖u_mixed − u_FP64‖∞ / ‖u_FP64‖∞` (fail if it
     exceeds `100·eps(TKrylov)`)
3. **Performance regression** test: track wall-clock for tutorial
   naca at fixed mesh + degree. Fail-on-regression for the default
   policy; expected-improvement for `LowPrecPreconditioner` and
   `MixedIR`.
4. **Stress test** at end of phase 5: re-run the test cases pinned in
   `tests/` (existing CTest matrix). Run each at all three policies
   (`DefaultPrecision`, `LowPrecPreconditioner`, `MixedIR`).
5. **GPU + MPI cross-check** (phases 5+): run on `lanka-dgx0` with 2
   MPI ranks per node × 2 nodes; compare residual norms rank-by-rank
   to verify no precision-dependent divergence.

---

## 5. Risks & mitigations

- **Compile-time fan-out from templating structs.** Mitigation: keep
  the policy a *small* finite set (`DefaultPrecision`,
  `LowPrecPreconditioner`, `MixedIR`) and explicitly-instantiate
  those in the same `.cpp` files that today instantiate
  `CDiscretization<AbiAdapter>`. The user-Model path then triggers
  template-instantiation cost only for their chosen policy.
- **Stiff problems lose digits in FP32 preconditioner.** Mitigation:
  policy is per-instance, not per-build — the user can fall back to
  `DefaultPrecision` for stiff PDE systems without rebuilding.
- **Codegen path stays FP64 in phases 0-5.** That's fine — phase 3
  and phase 4 only need the *FEM-internal* matrices in FP32, not the
  generated model kernels. Phase 6 is optional.
- **MPI halo precision divergence between ranks.** Mitigation: the
  `mpi_type<THalo>()` trait makes the type explicit; every MPI call
  uses the same trait, so all ranks agree by construction.
- **`MPI_Allreduce` of FP32 dot products has worse accuracy than
  FP64.** Mitigation: phase 4 spec already requires reduction
  accumulator in FP64; the `MPI_Allreduce` uses FP64 even if `V[i]`
  is FP32.
- **Enzyme path under FP32.** Mitigation: keep
  `HAVE_ENZYME ⇒ TResid=double` as a guard in phase 5 until Enzyme
  FP32 differentiation is verified end-to-end.

---

## 6. Deliverable order (a single PR could ship phases 0-1)

Suggested PR cadence:

- **PR-A**: Phase 0 + Phase 1 (scaffolding + struct templating, no
  behavior change). Smallest review surface for the riskiest
  structural change.
- **PR-B**: Phase 2 (BLAS / Kokkos templating). Stand-alone, fully
  testable.
- **PR-C**: Phase 3 (FP32 preconditioner). First user-visible win.
- **PR-D**: Phase 4 (FP32 Krylov). Adds the reduction-precision
  machinery; testable with phase 3 disabled.
- **PR-E**: Phase 5 (mixed-IR Newton). End-state for the
  "all-defaults" mixed-precision config.
- **PR-F**: Phase 6 + Phase 7 (model + halo). Optional finishers.

Each PR shifts a single coordinate and ships its own regression test
update.
