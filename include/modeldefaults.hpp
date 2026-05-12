// SPDX-License-Identifier: see LICENSE
//
// <exasim/model.hpp> — the Model contract.
//
// A user PDE model is a C++ struct that supplies pointwise math (and,
// for HDG, hand-written pointwise Jacobians) which the templated
// kernels in <exasim/kernels/*.hpp> call inside `Kokkos::parallel_for`
// to evaluate volume / surface integrals at quadrature points.
//
// There is no autodiff and no DSL. Users write `flux(f, x, uq, w, mu, t)`
// in plain C++ and, where the discretization needs derivatives, write
// `flux_jac_uq(f_uq, …)` themselves. exasim::ModelDefaults<Self> (CRTP
// base) supplies zero-fill no-op implementations of every optional
// method so users override only what their PDE actually has.
//
// Boundary identification is NOT part of the Model contract
// ---------------------------------------------------------
// The user's `M::fbou(...)`, `M::ubou(...)`, `M::qoi_boundary(...)`
// methods receive an integer `ib` boundary tag and dispatch on it:
//
//   void fbou(dstype fb[], int ib, …) {
//       if (ib == 1)      { /* wall */ }
//       else if (ib == 2) { /* inflow */ }
//       …
//   }
//
// The mapping from physical geometry to the integer `ib` happens at
// the runtime / preprocessing layer, not here. Users still write the
// usual `pdeapp.txt` config:
//
//   boundaryconditions  = [1, 2, 1, 3];
//   boundaryexpressions = ["abs(y)<1e-8", "x>9-1e-2", …];
//
// Preprocessing's `tinyexpr` evaluator parses the expression strings,
// tags each mesh face with the matching integer, and the kernels see
// only the resulting `int ib`. Curved-boundary expressions and
// periodic boundary matching are also preprocessing-layer concerns.
//
// Storage layout convention
// -------------------------
// The kernels see SoA buffers laid out as:
//   xdg [k*ng + i]   for k in [0, nd),       (spatial coords)
//   udg [k*ng + i]   for k in [0, Nq),       (u and ∇u, Nq = ncu*(1+nd))
//   wdg [k*ng + i]   for k in [0, ncw),      (auxiliary scalar fields)
//   param[k]         for k in [0, nparam),   (physics params; not per-i)
//
// At each quadrature point i, the kernel gathers small fixed-size
// arrays and passes them to the model's pointwise function:
//
//   dstype x [nd], uq[Nq], w[ncw];
//   for k in [0, nd):  x [k] = xdg[k*ng + i];
//   for k in [0, Nq):  uq[k] = udg[k*ng + i];
//   for k in [0, ncw): w [k] = wdg[k*ng + i];
//   dstype f_local[ncu*nd];
//   Model::flux(f_local, x, uq, w, param, t);
//   for k in [0, ncu*nd): f[k*ng + i] = f_local[k];
//
// Jacobian outputs use a row-major convention with the output index
// outer:
//   f_uq [ i * Nq  + j ] = ∂f_i / ∂uq_j
//   f_w  [ i * ncw + j ] = ∂f_i / ∂w_j
//   fb_uh[ i * ncu + j ] = ∂fb_i / ∂uh_j   (boundary trace Jacobian)
//
// Required compile-time configuration
// -----------------------------------
// Every Model must publish:
//
//   static constexpr int nd;       // spatial dimension (1, 2, or 3)
//   static constexpr int ncu;      // # primary unknowns
//   static constexpr int ncw;      // # auxiliary scalar fields (often 0)
//   static constexpr int nparam;   // # physics parameters
//
//   // Optional (default 0 from ModelDefaults<Self>):
//   static constexpr int nco;      // # "other DG" auxiliary fields
//                                  // (the `v` / `odg` arrays in
//                                  // libpdemodel and text2code's DSL)
//
//   // Derived (provided by ModelDefaults<Self>; users may override):
//   static constexpr int Nq = ncu * (1 + nd);
//
// Naming alignment with the text2code DSL (pdemodel.txt) and the
// libpdemodel.hpp ABI:
//
//   text2code DSL    libpdemodel ABI    Model contract (this file)
//   ─────────────────────────────────────────────────────────────
//   x                xdg                x[]    (spatial coords)
//   uq               udg                uq[]   ((u, ∇u) packed)
//   v                odg                — (only nco; not pointwise)
//   w                wdg                w[]    (auxiliary scalars)
//   uhat             uhg                uh[]   (HDG trace)
//   n                nlg                n[]    (face normal)
//   tau              tau                tau[]  (stabilization)
//   mu               param              mu[]   (physics params)
//   eta              uinf               uinf[] (free-stream / reference)
//   uext             uext               uext[] (Fext only; multi-domain)
//   t                time               t      (current time)
//
// Required pointwise functions
// ----------------------------
// The minimal set every Model must define (no defaults):
//
//   KOKKOS_INLINE_FUNCTION static
//   void flux  (dstype f [/*ncu*nd*/], const dstype x[/*nd*/],
//               const dstype uq[/*Nq*/], const dstype w[/*ncw*/],
//               const dstype mu[/*nparam*/], const dstype uinf[/*ncu*/],
//               dstype t);
//
//   KOKKOS_INLINE_FUNCTION static
//   void initu (dstype ui[/*ncu*/], const dstype x[/*nd*/],
//               const dstype uinf[/*ncu*/], const dstype mu[/*nparam*/]);
//
// Volume pointwise methods (`flux`, `source`, `sourcew`, `tdfunc`,
// `avfield`, `eos`, `eos_du`, `eos_dw`) all take args
//   (out, x, uq, w, mu, uinf, t)
// after the output buffer. `uinf` is pointer-passed and may be
// nullptr — methods that need free-stream values dereference at
// their own risk.
//
// Optional pointwise functions (default = zero-fill, via ModelDefaults):
//
//   source, sourcew, fbou, ubou, fhat, uhat, stab,
//   tdfunc, eos, eos_du, eos_dw, avfield,
//   init{q, udg, wdg, odg}, monitor, output,
//   vis_scalars, vis_vectors, vis_tensors,
//   qoi_volume, qoi_boundary
//
// HDG Jacobians (required iff Discretization::HDG is selected; default
// = zero-fill, so users get a compile error / zero-residual Jacobian
// only if they forget):
//
//   flux_jac_uq, flux_jac_w,
//   source_jac_uq, source_jac_w,
//   fbou_jac_uq, fbou_jac_uh, fbou_jac_w,
//   fhat_jac_uq, fhat_jac_uh, fhat_jac_w,
//   stab_jac_uq, stab_jac_uh,
//   eos_jac_uq, eos_jac_w
//
// All Jacobians are pointwise and `KOKKOS_INLINE_FUNCTION static`. See
// the per-kernel headers in <exasim/kernels/> for exact signatures.

#pragma once

#include "driver_abi.hpp"

template <class Self>
struct ModelDefaults {
    // Default count of "other DG" auxiliary fields (the `v` / `odg`
    // arrays in libpdemodel and the text2code DSL `vectors v(nco)`).
    // Most PDEs leave this at 0; only models with an external scalar
    // field (artificial viscosity, level-set, … plumbed through `odg`)
    // override.
    static constexpr int nco = 0;

    // Helper: zero-fill an output buffer of size N at compile time.
    template <int N>
    KOKKOS_INLINE_FUNCTION static void zero_fill_(dstype f[]) {
        for (int k = 0; k < N; ++k) f[k] = 0.0;
    }

    // ---- Volume terms (default: zero) ----
    //
    // Volume pointwise methods take 8 args after the output buffer:
    //   (out, x[nd], uq[Nq], v[nco], w[ncw], mu[nparam], uinf[], t)
    // `v` is the auxiliary "other DG" field (pdemodel.txt's `vectors v(nco)`,
    // libpdemodel's `odg`). `uinf` is pointer-passed and may be nullptr —
    // methods that need free-stream values dereference at their own risk.

    KOKKOS_INLINE_FUNCTION static
    void source(dstype s[],
                const dstype /*x*/[],  const dstype /*uq*/[],
                const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(s);
    }

    KOKKOS_INLINE_FUNCTION static
    void sourcew(dstype sw[],
                 const dstype /*x*/[],  const dstype /*uq*/[],
                 const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                 const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) zero_fill_<Self::ncw>(sw);
    }

    KOKKOS_INLINE_FUNCTION static
    void tdfunc(dstype m[],
                const dstype /*x*/[],  const dstype /*uq*/[],
                const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                const dstype /*uinf*/[], dstype /*t*/) {
        // Identity mass weighting by default.
        for (int k = 0; k < Self::ncu; ++k) m[k] = 1.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void avfield(dstype av[],
                 const dstype /*x*/[],  const dstype /*uq*/[],
                 const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                 const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(av);
    }

    // ---- Boundary terms (default: zero) ----

    // Boundary kernels: pointwise functions seen by every kernel call from
    // <exasim/kernels/boundary.hpp>. Args (after ub/fb output and `ib` tag):
    //   x[nd], uq[Nq], v[nco], w[ncw], uh[ncu], n[nd], tau[ncu], mu[nparam], uinf, t
    KOKKOS_INLINE_FUNCTION static
    void ubou(dstype ub[], int /*ib*/,
              const dstype /*x*/[],  const dstype /*uq*/[],
              const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
              const dstype /*n*/[],  const dstype /*tau*/[],
              const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(ub);
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou(dstype fb[], int /*ib*/,
              const dstype /*x*/[],  const dstype /*uq*/[],
              const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
              const dstype /*n*/[],  const dstype /*tau*/[],
              const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(fb);
    }

    // HDG boundary residual — distinct from `fbou`.
    //
    // In every example PDE in apps/, `Fbou` (LDG path) and `FbouHdg`
    // (HDG path) are *different math*: `Fbou` is the boundary flux
    // contribution in a volume residual; `FbouHdg` is the trace-side
    // residual that pins the hybrid trace `uhat`. Same signature
    // shape; the HDG kernels in <exasim/kernels/boundary.hpp> route
    // through `fbou_hdg` and its three companion Jacobians:
    //   fbou_hdg_jac_uq, fbou_hdg_jac_w, fbou_hdg_jac_uh
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(dstype fb[], int /*ib*/,
                  const dstype /*x*/[],  const dstype /*uq*/[],
                  const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                  const dstype /*n*/[],  const dstype /*tau*/[],
                  const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(fb);
    }

    // Interface kernels see both sides (uq1, uq2, v1, v2, w1, w2) plus the
    // trace uh, normal n, stabilization tau, and physics params mu. uinf is
    // pointer-passed and may be nullptr — model methods that need
    // free-stream values dereference at their own risk.
    KOKKOS_INLINE_FUNCTION static
    void fhat(dstype fh[],     const dstype /*x*/[],
              const dstype /*uq1*/[],  const dstype /*uq2*/[],
              const dstype /*v1*/[],   const dstype /*v2*/[],
              const dstype /*w1*/[],   const dstype /*w2*/[],
              const dstype /*uh*/[],   const dstype /*n*/[],   const dstype /*tau*/[],
              const dstype /*mu*/[],   const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(fh);
    }

    KOKKOS_INLINE_FUNCTION static
    void uhat(dstype uh[],     const dstype /*x*/[],
              const dstype /*uq1*/[],  const dstype /*uq2*/[],
              const dstype /*v1*/[],   const dstype /*v2*/[],
              const dstype /*w1*/[],   const dstype /*w2*/[],
              const dstype /*trace*/[],const dstype /*n*/[],   const dstype /*tau*/[],
              const dstype /*mu*/[],   const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(uh);
    }

    KOKKOS_INLINE_FUNCTION static
    void stab(dstype tau[],    const dstype /*x*/[],
              const dstype /*uq1*/[],  const dstype /*uq2*/[],
              const dstype /*v1*/[],   const dstype /*v2*/[],
              const dstype /*w1*/[],   const dstype /*w2*/[],
              const dstype /*uh*/[],   const dstype /*n*/[],   const dstype /*tau_in*/[],
              const dstype /*mu*/[],   const dstype /*uinf*/[], dstype /*t*/) {
        for (int k = 0; k < Self::ncu; ++k) tau[k] = 1.0;
    }

    // ---- Equation of state (default: identity / zero) ----

    KOKKOS_INLINE_FUNCTION static
    void eos(dstype e[],
             const dstype /*x*/[],  const dstype /*uq*/[],
             const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
             const dstype /*uinf*/[], dstype /*t*/) {
        zero_fill_<Self::ncu>(e);
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_du(dstype ed[],
                const dstype /*x*/[],  const dstype /*uq*/[],
                const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                const dstype /*uinf*/[], dstype /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) ed[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_dw(dstype ew[],
                const dstype /*x*/[],  const dstype /*uq*/[],
                const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) ew[k] = 0.0;
        }
    }

    // ---- Initial conditions (default: zero, except `initu` which is required) ----

    KOKKOS_INLINE_FUNCTION static
    void initq(dstype q[], const dstype[], const dstype[], const dstype[]) {
        zero_fill_<Self::ncu * Self::nd>(q);
    }

    KOKKOS_INLINE_FUNCTION static
    void initwdg(dstype w[], const dstype[], const dstype[], const dstype[]) {
        if constexpr (Self::ncw > 0) zero_fill_<Self::ncw>(w);
    }

    KOKKOS_INLINE_FUNCTION static
    void initudg(dstype udg[], const dstype[], const dstype[], const dstype[]) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        zero_fill_<Nq>(udg);
    }

    KOKKOS_INLINE_FUNCTION static
    void initodg(dstype odg[], const dstype[], const dstype[], const dstype[]) {
        // No default; `nco` is determined by the discretization, not by Self.
        // Concrete kernels pass the correct count in.
        (void)odg;
    }

    KOKKOS_INLINE_FUNCTION static
    void cpuinitu(dstype u[], const dstype x[], const dstype uinf[],
                  const dstype mu[]) {
        Self::initu(u, x, uinf, mu);
    }

    KOKKOS_INLINE_FUNCTION static
    void cpuinitq(dstype q[], const dstype x[], const dstype uinf[],
                  const dstype mu[]) {
        Self::initq(q, x, uinf, mu);
    }

    KOKKOS_INLINE_FUNCTION static
    void cpuinitudg(dstype udg[], const dstype x[], const dstype uinf[],
                    const dstype mu[]) {
        Self::initudg(udg, x, uinf, mu);
    }

    KOKKOS_INLINE_FUNCTION static
    void cpuinitwdg(dstype w[], const dstype x[], const dstype uinf[],
                    const dstype mu[]) {
        Self::initwdg(w, x, uinf, mu);
    }

    KOKKOS_INLINE_FUNCTION static
    void cpuinitodg(dstype odg[], const dstype x[], const dstype uinf[],
                    const dstype mu[]) {
        Self::initodg(odg, x, uinf, mu);
    }

    // ---- Visualization, QoI, monitor/output (default: empty / zero) ----
    //
    // Volume-shape methods (`vis_*`, `qoi_volume`, `monitor`, `output`)
    // take the same args as `flux`/`source`:
    //   (out, x, uq, w, mu, uinf, t)
    // `qoi_boundary` adds an `ib` tag and the boundary-trace args
    // (uh, n, tau).

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                     const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                     const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                     const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                     const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_tensors(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                     const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                     const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                    const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                    const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(dstype[], int /*ib*/,
                      const dstype /*x*/[],  const dstype /*uq*/[],
                      const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                      const dstype /*n*/[],  const dstype /*tau*/[],
                      const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void monitor(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                 const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                 const dstype /*uinf*/[], dstype /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void output(dstype[], const dstype /*x*/[], const dstype /*uq*/[],
                const dstype /*v*/[], const dstype /*w*/[], const dstype /*mu*/[],
                const dstype /*uinf*/[], dstype /*t*/) { }

    // ---- HDG Jacobians (default: zero) ----
    //
    // For HDG models these MUST be overridden by the user — the
    // defaults give zero Jacobians, which makes the linearized system
    // zero and Newton fail to converge. The defaults exist only so
    // that LDG models (which don't call into the *_jac_* paths) don't
    // need to write empty stubs.

    // Volume Jacobians: same arg shape as the value methods —
    //   (out, x, uq, v, w, mu, uinf, t)
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(dstype f_uq[],
                     const dstype /*x*/[],  const dstype /*uq*/[],
                     const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                     const dstype /*uinf*/[], dstype /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Self::nd * Nq; ++k) f_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_w(dstype f_w[],
                    const dstype /*x*/[],  const dstype /*uq*/[],
                    const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                    const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::nd * Self::ncw; ++k) f_w[k] = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_uq(dstype s_uq[],
                       const dstype /*x*/[],  const dstype /*uq*/[],
                       const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                       const dstype /*uinf*/[], dstype /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) s_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_w(dstype s_w[],
                      const dstype /*x*/[],  const dstype /*uq*/[],
                      const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*mu*/[],
                      const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) s_w[k] = 0.0;
        }
    }

    // LDG-path boundary Jacobians. These are not currently consumed by
    // any kernel — text2code's libpdemodel.hpp ABI has `KokkosFbou`
    // (LDG, value-only) and `HdgFbou` (HDG, value + Jacobians) but no
    // LDG-path Jacobians. Kept in the contract for completeness in
    // case a future numerical scheme needs them.
    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_uq(dstype fb_uq[], int /*ib*/,
                     const dstype /*x*/[],  const dstype /*uq*/[],
                     const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                     const dstype /*n*/[],  const dstype /*tau*/[],
                     const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) fb_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_uh(dstype fb_uh[], int /*ib*/,
                     const dstype /*x*/[],  const dstype /*uq*/[],
                     const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                     const dstype /*n*/[],  const dstype /*tau*/[],
                     const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        for (int k = 0; k < Self::ncu * Self::ncu; ++k) fb_uh[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_w(dstype fb_w[], int /*ib*/,
                    const dstype /*x*/[],  const dstype /*uq*/[],
                    const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                    const dstype /*n*/[],  const dstype /*tau*/[],
                    const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) fb_w[k] = 0.0;
        }
    }

    // HDG-path boundary Jacobians (companions to `fbou_hdg`). These
    // ARE consumed by `hdg_fbou_kernel<M>` in
    // <exasim/kernels/boundary.hpp>: the user MUST override them to
    // get a converging HDG Newton solve. Defaults zero-fill so that
    // LDG-only models compile without writing empty stubs.
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uq(dstype fb_uq[], int /*ib*/,
                         const dstype /*x*/[],  const dstype /*uq*/[],
                         const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                         const dstype /*n*/[],  const dstype /*tau*/[],
                         const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) fb_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(dstype fb_uh[], int /*ib*/,
                         const dstype /*x*/[],  const dstype /*uq*/[],
                         const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                         const dstype /*n*/[],  const dstype /*tau*/[],
                         const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        for (int k = 0; k < Self::ncu * Self::ncu; ++k) fb_uh[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_w(dstype fb_w[], int /*ib*/,
                        const dstype /*x*/[],  const dstype /*uq*/[],
                        const dstype /*v*/[],  const dstype /*w*/[],  const dstype /*uh*/[],
                        const dstype /*n*/[],  const dstype /*tau*/[],
                        const dstype /*mu*/[], const dstype /*uinf*/[], dstype /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) fb_w[k] = 0.0;
        }
    }
    // (fhat_jac_*, stab_jac_*, eos_jac_* follow the same shape; added
    // when the corresponding kernel template is wired up.)
};
