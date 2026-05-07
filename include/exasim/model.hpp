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
//   void fbou(double fb[], int ib, …) {
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
//   double x [nd], uq[Nq], w[ncw];
//   for k in [0, nd):  x [k] = xdg[k*ng + i];
//   for k in [0, Nq):  uq[k] = udg[k*ng + i];
//   for k in [0, ncw): w [k] = wdg[k*ng + i];
//   double f_local[ncu*nd];
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
//   void flux  (double f [/*ncu*nd*/], const double x[/*nd*/],
//               const double uq[/*Nq*/], const double w[/*ncw*/],
//               const double mu[/*nparam*/], const double uinf[/*ncu*/],
//               double t);
//
//   KOKKOS_INLINE_FUNCTION static
//   void initu (double ui[/*ncu*/], const double x[/*nd*/],
//               const double uinf[/*ncu*/], const double mu[/*nparam*/]);
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

#include <Kokkos_Core.hpp>

#include "common.h"   // dstype + Int

namespace exasim {

// Discretization tag — published as `Self::disc` by the user model so
// kernels can `if constexpr` on it. Defaults to LDG (Local DG) since
// it requires no Jacobians.
enum class Discretization {
    LDG,   // Local DG: only flux/source/initu values needed.
    HDG    // Hybridized DG: requires hand-written pointwise Jacobians.
};

// CRTP base supplying zero-fill no-op defaults for every optional
// method on the Model contract. Users inherit and override only what
// their PDE has.
//
//   struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
//       static constexpr int nd = 2, ncu = 1, ncw = 0, nparam = 1;
//       static constexpr auto disc = exasim::Discretization::HDG;
//
//       KOKKOS_INLINE_FUNCTION static
//       void flux(double f[], const double[], const double uq[],
//                 const double[], const double mu[], double) {
//           f[0] = mu[0]*uq[1];
//           f[1] = mu[0]*uq[2];
//       }
//
//       KOKKOS_INLINE_FUNCTION static
//       void flux_jac_uq(double f_uq[], const double[], const double[],
//                        const double[], const double mu[], double) {
//           constexpr int Nq = ncu*(1+nd);     // 3
//           for (int k = 0; k < ncu*nd*Nq; ++k) f_uq[k] = 0;
//           f_uq[0*Nq + 1] = mu[0];
//           f_uq[1*Nq + 2] = mu[0];
//       }
//
//       // ... source, ubou, initu — explicit overrides for what's used.
//   };
template <class Self>
struct ModelDefaults {
    // Default discretization tag. Override by `static constexpr auto
    // disc = exasim::Discretization::HDG;` in the derived struct.
    static constexpr Discretization disc = Discretization::LDG;

    // Default count of "other DG" auxiliary fields (the `v` / `odg`
    // arrays in libpdemodel and the text2code DSL `vectors v(nco)`).
    // Most PDEs leave this at 0; only models with an external scalar
    // field (artificial viscosity, level-set, … plumbed through `odg`)
    // override.
    static constexpr int nco = 0;

    // Helper: zero-fill an output buffer of size N at compile time.
    template <int N>
    KOKKOS_INLINE_FUNCTION static void zero_fill_(double f[]) {
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
    void source(double s[],
                const double /*x*/[],  const double /*uq*/[],
                const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(s);
    }

    KOKKOS_INLINE_FUNCTION static
    void sourcew(double sw[],
                 const double /*x*/[],  const double /*uq*/[],
                 const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                 const double /*uinf*/[], double /*t*/) {
        if constexpr (Self::ncw > 0) zero_fill_<Self::ncw>(sw);
    }

    KOKKOS_INLINE_FUNCTION static
    void tdfunc(double m[],
                const double /*x*/[],  const double /*uq*/[],
                const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        // Identity mass weighting by default.
        for (int k = 0; k < Self::ncu; ++k) m[k] = 1.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void avfield(double av[],
                 const double /*x*/[],  const double /*uq*/[],
                 const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                 const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(av);
    }

    // ---- Boundary terms (default: zero) ----

    // Boundary kernels: pointwise functions seen by every kernel call from
    // <exasim/kernels/boundary.hpp>. Args (after ub/fb output and `ib` tag):
    //   x[nd], uq[Nq], v[nco], w[ncw], uh[ncu], n[nd], tau[ncu], mu[nparam], uinf, t
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(ub);
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
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
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[],  const double /*uq*/[],
                  const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                  const double /*n*/[],  const double /*tau*/[],
                  const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(fb);
    }

    // Interface kernels see both sides (uq1, uq2, v1, v2, w1, w2) plus the
    // trace uh, normal n, stabilization tau, and physics params mu. uinf is
    // pointer-passed and may be nullptr — model methods that need
    // free-stream values dereference at their own risk.
    KOKKOS_INLINE_FUNCTION static
    void fhat(double fh[],     const double /*x*/[],
              const double /*uq1*/[],  const double /*uq2*/[],
              const double /*v1*/[],   const double /*v2*/[],
              const double /*w1*/[],   const double /*w2*/[],
              const double /*uh*/[],   const double /*n*/[],   const double /*tau*/[],
              const double /*mu*/[],   const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(fh);
    }

    KOKKOS_INLINE_FUNCTION static
    void uhat(double uh[],     const double /*x*/[],
              const double /*uq1*/[],  const double /*uq2*/[],
              const double /*v1*/[],   const double /*v2*/[],
              const double /*w1*/[],   const double /*w2*/[],
              const double /*trace*/[],const double /*n*/[],   const double /*tau*/[],
              const double /*mu*/[],   const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(uh);
    }

    KOKKOS_INLINE_FUNCTION static
    void stab(double tau[],    const double /*x*/[],
              const double /*uq1*/[],  const double /*uq2*/[],
              const double /*v1*/[],   const double /*v2*/[],
              const double /*w1*/[],   const double /*w2*/[],
              const double /*uh*/[],   const double /*n*/[],   const double /*tau_in*/[],
              const double /*mu*/[],   const double /*uinf*/[], double /*t*/) {
        for (int k = 0; k < Self::ncu; ++k) tau[k] = 1.0;
    }

    // ---- Equation of state (default: identity / zero) ----

    KOKKOS_INLINE_FUNCTION static
    void eos(double e[],
             const double /*x*/[],  const double /*uq*/[],
             const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
             const double /*uinf*/[], double /*t*/) {
        zero_fill_<Self::ncu>(e);
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_du(double ed[],
                const double /*x*/[],  const double /*uq*/[],
                const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) ed[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_dw(double ew[],
                const double /*x*/[],  const double /*uq*/[],
                const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) ew[k] = 0.0;
        }
    }

    // ---- Initial conditions (default: zero, except `initu` which is required) ----

    KOKKOS_INLINE_FUNCTION static
    void initq(double q[], const double[], const double[], const double[]) {
        zero_fill_<Self::ncu * Self::nd>(q);
    }

    KOKKOS_INLINE_FUNCTION static
    void initwdg(double w[], const double[], const double[], const double[]) {
        if constexpr (Self::ncw > 0) zero_fill_<Self::ncw>(w);
    }

    KOKKOS_INLINE_FUNCTION static
    void initudg(double udg[], const double[], const double[], const double[]) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        zero_fill_<Nq>(udg);
    }

    KOKKOS_INLINE_FUNCTION static
    void initodg(double odg[], const double[], const double[], const double[]) {
        // No default; `nco` is determined by the discretization, not by Self.
        // Concrete kernels pass the correct count in.
        (void)odg;
    }

    // ---- Visualization, QoI, monitor/output (default: empty / zero) ----
    //
    // Volume-shape methods (`vis_*`, `qoi_volume`, `monitor`, `output`)
    // take the same args as `flux`/`source`:
    //   (out, x, uq, w, mu, uinf, t)
    // `qoi_boundary` adds an `ib` tag and the boundary-trace args
    // (uh, n, tau).

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(double[], const double /*x*/[], const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(double[], const double /*x*/[], const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_tensors(double[], const double /*x*/[], const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(double[], const double /*x*/[], const double /*uq*/[],
                    const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                    const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(double[], int /*ib*/,
                      const double /*x*/[],  const double /*uq*/[],
                      const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                      const double /*n*/[],  const double /*tau*/[],
                      const double /*mu*/[], const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void monitor(double[], const double /*x*/[], const double /*uq*/[],
                 const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                 const double /*uinf*/[], double /*t*/) { }

    KOKKOS_INLINE_FUNCTION static
    void output(double[], const double /*x*/[], const double /*uq*/[],
                const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) { }

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
    void flux_jac_uq(double f_uq[],
                     const double /*x*/[],  const double /*uq*/[],
                     const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Self::nd * Nq; ++k) f_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_w(double f_w[],
                    const double /*x*/[],  const double /*uq*/[],
                    const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                    const double /*uinf*/[], double /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::nd * Self::ncw; ++k) f_w[k] = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_uq(double s_uq[],
                       const double /*x*/[],  const double /*uq*/[],
                       const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                       const double /*uinf*/[], double /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) s_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_w(double s_w[],
                      const double /*x*/[],  const double /*uq*/[],
                      const double /*v*/[],  const double /*w*/[],  const double /*mu*/[],
                      const double /*uinf*/[], double /*t*/) {
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
    void fbou_jac_uq(double fb_uq[], int /*ib*/,
                     const double /*x*/[],  const double /*uq*/[],
                     const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                     const double /*n*/[],  const double /*tau*/[],
                     const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) fb_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_uh(double fb_uh[], int /*ib*/,
                     const double /*x*/[],  const double /*uq*/[],
                     const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                     const double /*n*/[],  const double /*tau*/[],
                     const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        for (int k = 0; k < Self::ncu * Self::ncu; ++k) fb_uh[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_w(double fb_w[], int /*ib*/,
                    const double /*x*/[],  const double /*uq*/[],
                    const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                    const double /*n*/[],  const double /*tau*/[],
                    const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
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
    void fbou_hdg_jac_uq(double fb_uq[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double /*tau*/[],
                         const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) fb_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double /*tau*/[],
                         const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        for (int k = 0; k < Self::ncu * Self::ncu; ++k) fb_uh[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_w(double fb_w[], int /*ib*/,
                        const double /*x*/[],  const double /*uq*/[],
                        const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                        const double /*n*/[],  const double /*tau*/[],
                        const double /*mu*/[], const double /*uinf*/[], double /*t*/) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) fb_w[k] = 0.0;
        }
    }
    // (fhat_jac_*, stab_jac_*, eos_jac_* follow the same shape; added
    // when the corresponding kernel template is wired up.)
};

// ===========================================================================
// Deferred surface — recorded for v2 / multi-domain coupling
// ===========================================================================
//
// `apps/poisson/poisson2d/pdemodel2.txt` and the `HdgFint*` / `HdgFext*`
// libpdemodel.hpp ABI symbols document a multi-domain HDG path that the
// header-only port does not yet cover:
//
//   M::fint        (out, ib, x, uq, w, uh, n, tau, mu, uinf, t)
//   M::fint_jac_uq (...)
//   M::fint_jac_w  (...)
//   M::fint_jac_uh (...)
//   M::fext        (out, ib, x, uq, w, uh, n, uext, tau, mu, uinf, t)
//   M::fext_jac_uq (...)
//   M::fext_jac_w  (...)
//   M::fext_jac_uh (...)
//
// `Fint` is the HDG trace residual on internal model interfaces; `Fext`
// adds an extra input array `uext` carrying the trace from a
// neighbouring PDE model so the user can write coupling residuals like
// `fb[0] = uext[0] - uh[0]`.
//
// Out of scope for v1 because:
//   - single-domain HDG (Poisson, single Navier-Stokes, etc.) does not
//     reach these code paths;
//   - `uext` requires plumbing a foreign-model trace through the kernel
//     dispatch, which interacts with how multiple `CSolution<Model>`
//     instances are orchestrated by the runtime;
//   - none of the canonical apps in apps/poisson/* or apps/navierstokes/*
//     except poisson2d/pdemodel2.txt define `Fint`/`Fext`.
//
// When implementing: add `static constexpr bool has_external_coupling`
// to the contract, gate the `fext`-related kernel calls behind it, and
// extend `<exasim/kernels/boundary.hpp>` with `fint_kernel<M>` /
// `fext_kernel<M, OtherModel>` (the second template arg supplies the
// foreign model whose trace lands in `uext`).
//
// Also deferred to v2:
//   - `hessian` declarations from text2code's DSL: second-derivative
//     methods for Galerkin schemes that need them. None of the
//     reference apps populate `hessian …`.
//   - Multi-model dispatch via runtime `modelnumber` (text2code's
//     KokkosFlux1/KokkosFlux2 split). In the templated world a user
//     orchestrates by instantiating `CSolution<Model1>`,
//     `CSolution<Model2>` separately, so kernels see only one Model
//     at a time and `modelnumber` becomes a no-op.

// ---- Compile-time concept-style checks ----
//
// Static asserts that fire if a user struct is missing the required
// surface. Used inside kernel templates: `static_assert(is_model_v<M>, …);`

namespace detail {

template <class, class = void> struct has_required_constants : std::false_type {};
template <class M>
struct has_required_constants<M,
    std::void_t<
        decltype(M::nd), decltype(M::ncu),
        decltype(M::ncw), decltype(M::nparam)
    >> : std::true_type {};

} // namespace detail

template <class M>
inline constexpr bool is_model_v = detail::has_required_constants<M>::value;

} // namespace exasim
