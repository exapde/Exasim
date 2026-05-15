# Mixed-Precision Mode — Codebase Survey

Goal: identify where in Exasim a templated mixed-precision mode would
have the highest payoff, what data and code surfaces currently bake in
a single global precision, and where conversions would need to live.

This document is *research only*. The implementation plan lives in
`plans/mixed-precision-plan.md`.

---

## 1. How precision is wired today

### 1.1 A single global typedef
- `backend/Common/common.h:138-142`:
  ```cpp
  #ifdef USE_FLOAT
  typedef float dstype;
  #else
  typedef double dstype;   // double is default precision
  #endif
  ```
- `dstype` is referenced in **177 backend files** (`grep -l dstype`),
  including every struct in the FEM data model (`solstruct`,
  `resstruct`, `sysstruct`, `precondstruct`, `appstruct`,
  `masterstruct`, `meshstruct`, `tempstruct`, `commonstruct`,
  `wallmodelstruct`).
- Every BLAS wrapper, every Kokkos kernel, every codegen-emitted
  `Kokkos*` / `Hdg*` model routine takes `dstype*`. Switching between
  float and double is a *whole-program* compile-time flag.

### 1.2 BLAS / cuBLAS / hipBLAS dispatch
- `backend/Common/common.h:318-350` defines `cublasGEMM`, `cublasDOT`,
  … as preprocessor aliases:
  ```cpp
  #ifdef USE_FLOAT
  #define cublasGEMM cublasSgemm
  #else
  #define cublasGEMM cublasDgemm
  #endif
  ```
- `backend/Common/pblas.h` mirrors the same structure for every
  routine: `cpuNode2Gauss`, `Node2Gauss`, `Gauss2Node`,
  `gpuComputeInverse`, `cpuComputeInverse`, `Inverse`, `PDOT`,
  `PNORM`, `DOT`, `ArrayCopy`, `ArrayMultiplyScalar`, `ArrayAXPY`,
  `ArrayAXPBY`, `PGEMNV`, `PGEMTV`, `PGEMTM`, `PGEMNMStridedBached`.
  Each is a free function taking `dstype*`, dispatching to S/D BLAS by
  `#ifdef USE_FLOAT`.
- The handles (`cublasHandle_t`, `hipblasHandle_t`) are not typed by
  precision, so a single handle works in both modes.

### 1.3 MPI types
- 75 `MPI_DOUBLE` references across the backend. `pblas.h:311-316`,
  `pblas.h:557-572`, `pblas.h:678-693` use the same `#ifdef USE_FLOAT`
  pattern (`MPI_FLOAT` vs `MPI_DOUBLE`).
- Halo / Schur communication in `matvec.hpp`, `residual.hpp`,
  `setsysstruct.hpp`, `qequation.hpp`, `wequation.hpp` calls
  `MPI_Isend(... MPI_DOUBLE ...)` directly — a hardcoded
  `MPI_DOUBLE`. These will need to become traits-driven.

### 1.4 Codegen / text2code
- `text2code/text2code/CodeGenerator.cpp` emits all model files
  (`KokkosFlux.cpp`, `HdgFlux.cpp`, …) with `dstype*` signatures
  (e.g. line 1590: `"void " << funcname << "(dstype* f, …, const dstype time, …)"`).
- `backend/Model/libpdemodel.hpp` declares the **ABI boundary** for
  the generated code as `dstype*` (all 40+ `Kokkos*` / `Hdg*` /
  `cpuInit*` entries).
- The new templated model path (`include/exasim/drivers.hpp`,
  `include/exasim/kernels/*.hpp`) is still `dstype`-typed at the
  storage layer but inlines `M::flux(...)` on a stack-local `double[]`
  (see `kernels/flux.hpp:75-100`). So the *math* inside the kernel is
  on stack-local `double`, while the array storage matches `dstype`.
  This is significant — see §4.

### 1.5 Hand-written model path is first-class
- Per memory `feedback_handwritten_path`: the templated `M::flux` path
  is a peer of codegen. Any plan must support both: the codegen ABI
  (`libpdemodel.hpp`) and the templated `M::flux` kernels.

---

## 2. The solve, end to end (what mixed precision needs to know)

The HDG / LDG nonlinear solve runs the following chain (see
`backend/Solver/ptcsolver.hpp`, `gmres.hpp`, `matvec.hpp`):

```
PTCsolver / NonlinearSolver        (Newton outer loop)
  └─ LinearSolver                  (per iteration)
       ├─ evalResidual / hdgAssembleLinearSystem   →  fill sys.b, res.K/H/F
       ├─ (optional) ComputeInitialGuessAndPreconditioner   →  W, U, RB
       └─ GMRES
            ├─ evalMatVec     →  J·v  (LDG: finite-diff via Residual;
            │                          HDG: dense GEMM against res.H)
            ├─ ApplyPreconditioner
            │     ├─ LDG: ApplyMassInv on res.Minv
            │     └─ HDG: Block-Jacobi / Schwarz / Block-ILU0 on res.K
            ├─ (optional) ApplyPoly   →  polynomial preconditioner using sys.lam
            ├─ MGS / CGS              →  orthogonalize sys.v[m]
            ├─ Givens rotations       →  s, cs, sn
            └─ UpdateSolution          →  sys.x += sys.v · y
```

### 2.1 Per-stage data and where it lives
The structs in `backend/Common/common.h` group the per-stage data:

| Stage              | Struct          | Key arrays                                                                     |
|--------------------|-----------------|---------------------------------------------------------------------------------|
| Mesh / geometry    | `masterstruct`  | `shapegt`, `shapfgw`, `Xx`, `jac`, `gpe`, `gwe`, `gpf`, `gwf` (static)         |
| Solution state     | `solstruct`     | `udg`, `wdg`, `odg`, `uh`, `xdg`, `udgprev*`, `wprev*`                          |
| Residual / Jacobian| `resstruct`     | `Ru`, `Rh`, `Rq`, `Mass`, `Minv`, `K`, `H`, `F`, `D`, `B`, `C`, `G`, `E`, `Ri`, `Hi`, `Ki`, `Gi` |
| Krylov / outer     | `sysstruct`     | `x`, `u`, `r`, `b`, `v` (basis), `q`, `p`, `tempmem` (s, cs, sn, y, H), `lam`, `randvect`, prev solutions for BDF/DIRK |
| Preconditioner     | `precondstruct` | `W`, `U` (reduced basis), `ipiv`                                                |
| Halo buffers       | `tempstruct`    | `tempn`, `tempg`, `buffrecv`, `buffsend`, `bufffacerecv`, `bufffacesend`        |

Every field is `dstype*`. There is no per-array precision today.

### 2.2 Matrix-vector kernels — the hot spot
- **LDG** (`backend/Discretization/matvec.hpp:40-101`,
  `MatVec<M>`): builds `J·v` by finite difference,
  `(R(u+εv) − R(u))/ε`. Two full residual evaluations per outer
  iteration call. The `Residual` call is the most expensive step in
  the solve (it walks the entire element/face loop, builds geometry,
  calls every kernel in the model).
- **HDG** (`backend/Discretization/matvec.hpp:246-395`,
  `hdgMatVec<M>`): pre-assembled global Jacobian `res.H` (block
  matrix, size `(ncu·npf·nfe)²·ne`) multiplied against `v` via
  `PGEMNMStridedBached`. The matrix `H` and intermediate buffers are
  the largest GPU allocations in HDG runs.

### 2.3 Preconditioner kernels
- LDG: `ApplyMassInv` (`backend/Preconditioning/applymatrix.hpp:66-85`)
  is a single `Gauss2Node` GEMM (curved) or `ApplyJacInv + Gauss2Node`
  (straight). Bandwidth-bound.
- HDG: three modes (`backend/Preconditioning/preconditioner.hpp:272-325`):
  - Block-Jacobi: invert per-face `res.K` and apply.
  - Elemental Schwarz: `hdgMatVec` on a precomputed
    `(ncu·npf·nfe)²·ne` block.
  - Super-element Schwarz with Block-ILU0: forward/backward solve via
    `PGEMNMStridedBached` against `res.K`.
- The polynomial preconditioner (`backend/Solver/gmres.hpp:91-130`,
  `ApplyPoly`) wraps the base preconditioner with a Chebyshev / Leja
  iteration. Each step calls `evalMatVec` + `ApplyPreconditioner`.

### 2.4 Krylov accumulators
- `MGS` / `CGS` (`backend/Solver/getpoly.hpp:33-64`,
  `backend/Solver/gmres.hpp:81-88`): the `PDOT`s and `ArrayAXPY`s
  accumulating Hessenberg entries `H[k]` and renormalizing `V[m]`.
- `cpuApplyGivensRotation`, `cpuBackSolve`, `UpdateSolution`: small
  serial work on `tempmem` (host-pinned). Always cheap and almost
  always should stay FP64 — these are the entries that determine
  GMRES stability.

### 2.5 Reduced-basis preconditioner state
- `precond.W`, `precond.U` (`backend/Common/common.h:1603-1633`,
  `backend/Preconditioning/preconditioner.hpp:78-127`): basis of past
  Newton increments, size `ndof × RBdim`. Used to seed the initial
  guess and as a low-rank preconditioner.

---

## 3. Candidate components for mixed precision

These are sorted by *expected payoff / risk*. The exact numbers below
are upper bounds — actual speedups depend on roofline position.

### 3.1 Inner Krylov data in low precision (high payoff, low risk)
- Holding `sys.v` (the Krylov basis) and the residual workspace `Ru`
  in FP32 halves memory traffic for the `MatVec` and orthogonalization
  loops. These are the dominant GPU memory-traffic consumers in a
  GMRES restart cycle.
- Hessenberg entries `H`, Givens parameters `cs`, `sn`, `s`, `y` stay
  FP64 — they are tiny (`(nrest+1)²` total) and they govern
  convergence stability.
- Conversions: down-cast `sys.r → sys.v[0]` once per restart;
  `evalMatVec` output is cast back to FP32 before MGS.
- *Numerical concern*: MGS in FP32 is famously fragile beyond ~30
  iterations. Mitigations: keep `PDOT` accumulator in FP64 (the
  `MPI_Allreduce` already takes a scratch + global pair, easy to make
  FP64); switch from MGS to CGS-with-reorthog when Krylov is in FP32.

### 3.2 Preconditioner data in low precision (high payoff, low risk)
- `res.K`, `res.H`, `res.Minv`, `res.Mass2`, `res.Minv2`,
  `precond.W`, `precond.U` are touched once per linear iteration but
  not by the convergence test. The classic "FP32 preconditioner /
  FP64 GMRES" arrangement.
- Block-ILU0 / Block-Jacobi factorizations of `res.K` can be done and
  stored in FP32. The inversion of `npf·ncu`-sized diagonal blocks in
  `hdgBlockILU0` (`matvec.hpp:120-168`) is `O(ncf³)` per block and is
  precision-tolerant for well-conditioned local blocks.
- Conversions: down-cast `r → r_lo`, run preconditioner in FP32,
  up-cast result. A single `ArrayCast<TFrom,TTo>` kernel suffices.
- *Numerical concern*: stiff problems (strong AV, high Reynolds) may
  see GMRES iteration counts increase. Plan to keep precondition
  precision configurable per problem.

### 3.3 Jacobian-vector product (finite difference) in low precision (medium payoff, medium risk)
- The LDG `MatVec` evaluates `(R(u+εv) − R(u))/ε`. The
  `ε ≈ √machine_eps` cancellation already costs ~half precision; in
  FP64 you still get ~7-8 digits of `J·v`. In FP32 you lose to about 3
  digits, often enough for an *inner* iteration but not for *outer*
  Newton.
- Doable as: down-cast `u`, `v` to FP32, run `Residual` in FP32 (this
  means the *model evaluation* runs in FP32), demote `res.Ru` to
  FP32, do the subtraction and divide in FP32, up-cast the result.
- *Numerical concern*: with FP32 model evaluation, near-stagnation can
  drown in noise. Plan to use the `matvecOrder=2` central-difference
  form as the default for FP32 MatVec — error scales as `ε² ·
  ‖R‴‖∞`, much better than the FP32 `ε` truncation.

### 3.4 Model evaluation (Flux/Source/Fhat/Fbou) in low precision (medium payoff, model-dependent risk)
- Codegen and the templated `M::flux` kernels both already operate on
  small stack-local arrays (`double f_local[ncu*nd]` etc., see
  `kernels/flux.hpp`). Two ways to lower the *math* to FP32:
  1. Make the kernel template type parameter `T` instead of hard
     `double`. Stack arrays become `T f_local[…]`.
  2. Keep stack math in FP64 but make the *storage* I/O FP32 (load
     FP32, expand to FP64, compute, contract to FP32). This is the
     "TensorCore-style" path.
- Choice (2) is the safer first step (math is unaffected; you only
  pay the bandwidth dividend on udg/Ru I/O). Choice (1) needs
  per-equation sensitivity analysis — for example, Euler with
  pressure recovery from `(γ-1)(E - ½ρ|u|²)` is famously badly
  conditioned and *must* stay FP64.

### 3.5 MPI buffers in low precision (medium payoff, low risk)
- Halo packing/unpacking can be done in FP32 with negligible accuracy
  loss when the halo data is itself a Krylov basis vector or residual
  fragment (those buffers transfer per outer GMRES iteration).
- Switch `MPI_DOUBLE → MPI_FLOAT` everywhere depends on the buffer
  precision, not the global precision.

### 3.6 Geometry / static data in low precision (low payoff, low risk)
- `masterstruct` (shape functions, Gauss weights), curvature-related
  fields on `meshstruct` are static after construction. Holding them
  in FP32 saves cache footprint in geometry kernels (`Node2Gauss`,
  `Gauss2Node`, `ApplyXxJac`). Speedup is small for problems where
  geometry kernels aren't the bottleneck.

### 3.7 Mixed-precision iterative refinement at the Newton level (high payoff, classical pattern)
- Outer Newton: accumulate `b = -R(u)` in FP64, solve `J · dx = b` in
  FP32 (entire GMRES in FP32, including J·v, preconditioner, basis),
  apply update `u += dx` in FP64. This is the "GMRES-IR" / "TC-IR"
  pattern that gets ~1.7-2× speedup on the linear solve without
  touching outer Newton convergence.
- Maps cleanly onto Exasim: `sys.u`, `sys.b` stay FP64;
  `sys.x` (the Newton update vector), `sys.r`, `sys.v[…]`, and the
  whole `res.*` matrices are FP32.

---

## 4. Where conversions must live

A mixed-precision build must have well-defined precision boundaries.
Each boundary is a small set of cast operations.

### 4.1 Boundary: MatVec entry/exit
```
evalMatVec(Jv_lo, v_lo, u_hi, Ru_hi):    // lo = FP32, hi = FP64 (configurable)
    cast u_hi → u_lo
    cast v_hi → v_lo                   // (no-op if Krylov already lo)
    R(u_lo + ε·v_lo) − R(u_lo) → Jv_lo  // model in lo
    [Jv_lo stays lo; outer Givens promotes when needed]
```
Cast kernels needed: `ArrayCast<FP64,FP32>`, `ArrayCast<FP32,FP64>`
(both Kokkos `parallel_for` over `N`).

### 4.2 Boundary: ApplyPreconditioner entry/exit
- LDG (mass-inverse) and HDG (block solve) take `r` and produce
  `M⁻¹·r`. If `res.Minv`, `res.K`, etc. are stored in FP32 and `r` is
  FP64, down-cast `r` on entry; up-cast the result on exit.

### 4.3 Boundary: MPI halo
- Pack low-precision halo into low-precision buffer, send with
  matching `MPI_FLOAT`. Receiver up-casts as part of unpack if the
  consumer wants FP64.
- Sites: `matvec.hpp:282-300, 437-456, 585-602`,
  `residual.hpp` (similar pattern), `setsysstruct.hpp:111-201`,
  `qequation.hpp`, `wequation.hpp`.
- Plumbing: replace `MPI_DOUBLE` with a `mpi_type<T>()` helper.

### 4.4 Boundary: orthogonalization accumulators
- Inside `MGS` / `CGS` / `PDOT` (`pblas.h:278-321`,
  `getpoly.hpp:33-64`): even when `V[i]` is FP32, the inner product
  should be accumulated in FP64. Either use the FP64 BLAS even on
  FP32 vectors (Kahan-style; cuBLAS does not expose this directly) or
  write a small custom Kokkos reduction that holds the accumulator in
  FP64.

### 4.5 Boundary: codegen ABI
- Each generated `KokkosFlux(dstype*, …)` is bound to `dstype`. To
  support per-call precision selection we either:
  - Generate two variants (`KokkosFlux_float`, `KokkosFlux_double`)
    and dispatch at runtime, or
  - Template the generator so the model file emits `template <class T>
    void KokkosFlux(T* f, …)`, and the FEM side instantiates.
- The templated `M::flux` path (`kernels/flux.hpp`) is already poised
  for option 2 — the kernel template just needs to take a precision
  type parameter alongside `M`.

---

## 5. Constraints from the existing code

### 5.1 The class templates already exist
- `CDiscretization<M>`, `CSolver<M>`, `CPreconditioner<M>` already
  default to `M = AbiAdapter` (`include/exasim/detail/abi_adapter.hpp`).
  This is the same mechanism that today lets a user instantiate
  `CSolution<MyModel>` to bypass the codegen ABI.
- Adding a second template parameter `Precision = DefaultPrecision`
  is mechanically possible. The `EXASIM_DRIVER_CALL` macro
  (`include/exasim/detail/driver_dispatch.hpp`) only dispatches on
  `M`; we'd need to thread the precision policy through it too.

### 5.2 The struct fields are not templates
- `solstruct`, `resstruct`, etc. are *not* templates today. Every
  field is `dstype*`. To support multiple precisions in one process,
  we either:
  - Templatize each struct on a precision policy
    (`solstruct<Prec>`), which is the cleanest but a big change, or
  - Add parallel low-precision fields next to the existing ones
    (`dstype* udg`; `float* udg_lo`), with a flag for which is live.
    This is the surgical option.

### 5.3 The Inverse / Schur / batched-GEMM helpers cast inputs
- `pblas.h:163-203` (`gpuComputeInverse`) creates a `dstype**` array
  of pointers from a `dstype*` flat buffer. Type-correct
  templating works; the cast `(const dstype **)Ap_d` will need to
  become a `const T**`.
- `pblas.h:574-731` (commented-out `PGEMTV2`, `PGEMNM`) shows the
  authors already prototyped templated variants but kept the macro
  path. The templated path can be revived.

### 5.4 The `feedback_gpu_naming` memory rule
- Kokkos provides portability; CUDA + HIP share one CMake target. So
  any new cast/precision helpers must live inside Kokkos kernels (the
  same way `kokkosimpl.h` is structured), not as separate
  CUDA-or-HIP kernels.

### 5.5 The hand-written model is first-class
- Per `feedback_handwritten_path`: when writing the precision policy
  hooks, both the codegen ABI and the templated kernel path must be
  supported. The plan must give the user a way to declare *their*
  model's precision policy without regenerating codegen.

---

## 6. Open questions

- **Is the bottleneck arithmetic or bandwidth?** For HDG on GPU,
  `hdgMatVec` is dominated by `PGEMNMStridedBached` against `res.H`,
  which is memory-bandwidth-bound at typical problem sizes — that
  favors lowering `res.H` to FP32. For LDG, the bottleneck is the
  model evaluation in `Residual` — that favors lowering the model
  evaluation. Strategy should default differently per spatial scheme.
- **What does Enzyme do under FP32?** `matvec.hpp:54-58` has the
  Enzyme automatic-differentiation path under `HAVE_ENZYME`. Enzyme
  needs to be able to differentiate a FP32 residual; this is
  Enzyme-supported but adds a constraint.
- **Polynomial preconditioner under FP32?** `ApplyPoly` calls
  `evalMatVec` and `ApplyPreconditioner` inside a Chebyshev loop. If
  both are FP32, the Ritz value estimate from `getPoly` may lose
  enough digits to mis-place the interpolation nodes. Recommend
  computing Ritz values in FP64 even when applying the polynomial in
  FP32.
- **Does the codegen support templated emission today?** No — the
  generator hard-codes `dstype` in 30+ string literals. Touching this
  is a separable workstream.
