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
//   // Derived (provided by ModelDefaults<Self>; users may override):
//   static constexpr int Nq = ncu * (1 + nd);
//
// Required pointwise functions
// ----------------------------
// The minimal set every Model must define (no defaults):
//
//   KOKKOS_INLINE_FUNCTION static
//   void flux  (double f [/*ncu*nd*/], const double x[/*nd*/],
//               const double uq[/*Nq*/], const double w[/*ncw*/],
//               const double mu[/*nparam*/], double t);
//
//   KOKKOS_INLINE_FUNCTION static
//   void initu (double ui[/*ncu*/], const double x[/*nd*/],
//               const double uinf[/*ncu*/], const double mu[/*nparam*/]);
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

    // Helper: zero-fill an output buffer of size N at compile time.
    template <int N>
    KOKKOS_INLINE_FUNCTION static void zero_fill_(double f[]) {
        for (int k = 0; k < N; ++k) f[k] = 0.0;
    }

    // ---- Volume terms (default: zero) ----

    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double[], const double[], const double[],
                const double[], double) {
        zero_fill_<Self::ncu>(s);
    }

    KOKKOS_INLINE_FUNCTION static
    void sourcew(double sw[], const double[], const double[], const double[],
                 const double[], double) {
        if constexpr (Self::ncw > 0) zero_fill_<Self::ncw>(sw);
    }

    KOKKOS_INLINE_FUNCTION static
    void tdfunc(double m[], const double[], const double[], const double[],
                const double[], double) {
        // Identity mass weighting by default.
        for (int k = 0; k < Self::ncu; ++k) m[k] = 1.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void avfield(double av[], const double[], const double[], const double[],
                 const double[], double) {
        zero_fill_<Self::ncu>(av);
    }

    // ---- Boundary terms (default: zero) ----

    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/, const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], double) {
        zero_fill_<Self::ncu>(ub);
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/, const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], double) {
        zero_fill_<Self::ncu>(fb);
    }

    KOKKOS_INLINE_FUNCTION static
    void fhat(double fh[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], double) {
        zero_fill_<Self::ncu>(fh);
    }

    KOKKOS_INLINE_FUNCTION static
    void uhat(double uh[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], double) {
        zero_fill_<Self::ncu>(uh);
    }

    KOKKOS_INLINE_FUNCTION static
    void stab(double tau[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], const double[], const double[],
              const double[], const double[], double) {
        for (int k = 0; k < Self::ncu; ++k) tau[k] = 1.0;
    }

    // ---- Equation of state (default: identity / zero) ----

    KOKKOS_INLINE_FUNCTION static
    void eos(double e[], const double[], const double[], const double[],
             const double[], const double[], double) {
        zero_fill_<Self::ncu>(e);
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_du(double ed[], const double[], const double[], const double[],
                const double[], const double[], double) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) ed[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void eos_dw(double ew[], const double[], const double[], const double[],
                const double[], const double[], double) {
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

    // ---- Visualization & QoI (default: empty / zero) ----

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(double[], const double[], const double[], const double[],
                     const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(double[], const double[], const double[], const double[],
                     const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void vis_tensors(double[], const double[], const double[], const double[],
                     const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(double[], const double[], const double[], const double[],
                    const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(double[], int /*ib*/, const double[], const double[],
                      const double[], const double[], const double[], const double[],
                      const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void monitor(double[], const double[], const double[], const double[],
                 const double[], const double[], double) { }

    KOKKOS_INLINE_FUNCTION static
    void output(double[], const double[], const double[], const double[],
                const double[], const double[], double) { }

    // ---- HDG Jacobians (default: zero) ----
    //
    // For HDG models these MUST be overridden by the user — the
    // defaults give zero Jacobians, which makes the linearized system
    // zero and Newton fail to converge. The defaults exist only so
    // that LDG models (which don't call into the *_jac_* paths) don't
    // need to write empty stubs.

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double[], const double[], const double[],
                     const double[], double) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Self::nd * Nq; ++k) f_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_w(double f_w[], const double[], const double[], const double[],
                    const double[], double) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::nd * Self::ncw; ++k) f_w[k] = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_uq(double s_uq[], const double[], const double[], const double[],
                       const double[], double) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) s_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void source_jac_w(double s_w[], const double[], const double[], const double[],
                      const double[], double) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) s_w[k] = 0.0;
        }
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_uq(double fb_uq[], int /*ib*/, const double[], const double[],
                     const double[], const double[], const double[], const double[],
                     const double[], const double[], double) {
        constexpr int Nq = Self::ncu * (1 + Self::nd);
        for (int k = 0; k < Self::ncu * Nq; ++k) fb_uq[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_uh(double fb_uh[], int /*ib*/, const double[], const double[],
                     const double[], const double[], const double[], const double[],
                     const double[], const double[], double) {
        for (int k = 0; k < Self::ncu * Self::ncu; ++k) fb_uh[k] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_jac_w(double fb_w[], int /*ib*/, const double[], const double[],
                    const double[], const double[], const double[], const double[],
                    const double[], const double[], double) {
        if constexpr (Self::ncw > 0) {
            for (int k = 0; k < Self::ncu * Self::ncw; ++k) fb_w[k] = 0.0;
        }
    }
    // (fhat_jac_*, stab_jac_*, eos_jac_* follow the same shape; added
    // when the corresponding kernel template is wired up.)
};

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
