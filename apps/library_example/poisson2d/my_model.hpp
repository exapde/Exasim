// SPDX-License-Identifier: see Exasim LICENSE
//
// my_model.hpp — Poisson 2D PDE model, written by hand against the
// `<exasim/model.hpp>` contract. No DSL, no autodiff, no codegen.
//
// Solves   -∇·(μ∇u) = 2π² sin(πx) sin(πy)   on the unit square
//                u = 0                       on the boundary
//
// Exact solution:   u(x,y) = sin(πx) sin(πy)
// Physics param:    μ = mu[0]  (set to 1 in pdeapp.txt -> physicsparam = [1])
//
// The math here matches what text2code's SymEngine pipeline produces
// for apps/poisson/poisson2d/pdemodel.txt, byte-for-byte in evaluation
// order, so the templated path can reproduce the existing baseline
// numerics bit-for-bit when the runtime context is the same.

#pragma once

#include <exasim/model.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    // ---- Compile-time configuration -------------------------------
    static constexpr int nd     = 2;   // spatial dimension
    static constexpr int ncu    = 1;   // # primary unknowns
    static constexpr int ncw    = 0;   // # auxiliary scalar fields
    static constexpr int nco    = 0;   // # other DG fields (`v` / `odg`)
    static constexpr int nparam = 1;   // # physics parameters

    static constexpr auto disc = exasim::Discretization::HDG;

    // Derived: Nq = ncu*(1 + nd) = 3
    static constexpr int Nq = ncu * (1 + nd);

    // ---- Volume terms ---------------------------------------------

    // Flux f = μ ∇u
    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double /*x*/[], const double uq[],
              const double /*w*/[], const double mu[],
              const double /*uinf*/[], double /*t*/) {
        const double mu0  = mu[0];
        const double udg2 = uq[1];
        const double udg3 = uq[2];
        f[0] = mu0 * udg2;
        f[1] = mu0 * udg3;
    }

    // Source s = 2π² sin(πx) sin(πy)
    // The constant 1.973920880217872E+1 = 2*π² and matches text2code's
    // SymEngine output exactly.
    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double /*uq*/[],
                const double /*w*/[], const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        const double xdg1 = x[0];
        const double xdg2 = x[1];
        s[0] = Kokkos::sin(xdg1 * 3.141592653589793)
             * Kokkos::sin(xdg2 * 3.141592653589793)
             * 1.973920880217872E+1;
    }

    // tdfunc: empty body in text2code's KokkosTdfunc1 (steady problem,
    // torder=1, never called). The ModelDefaults<Poisson2D>::tdfunc
    // identity-fill default is fine — it's just not exercised at
    // runtime for this app.

    // ---- Required initial condition -------------------------------

    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double /*x*/[],
               const double /*uinf*/[], const double /*mu*/[]) {
        ui[0] = 0.0;
    }

    // ---- HDG path -------------------------------------------------

    // Volume-term Jacobian ∂f/∂uq:
    //
    // Layout (column-major; uq index is outer):
    //   f_uq[ j * (ncu*nd) + i ] = ∂f[i] / ∂uq[j]
    //
    // For Poisson 2D (ncu=1, nd=2, Nq=3 -> shape [ncu*nd=2 × Nq=3] = 6
    // entries):
    //   ∂f[0]/∂uq[1] = mu  ->  index = 1*2 + 0 = 2
    //   ∂f[1]/∂uq[2] = mu  ->  index = 2*2 + 1 = 5
    // All others zero.  This matches text2code's HdgFlux output exactly.
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double /*x*/[],
                     const double /*uq*/[], const double /*w*/[],
                     const double mu[], const double /*uinf*/[],
                     double /*t*/) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];   // ∂f[0]/∂uq[1]
        f_uq[2 * (ncu * nd) + 1] = mu[0];   // ∂f[1]/∂uq[2]
    }

    // Source has no uq dependence — defaults zero-fill from ModelDefaults.

    // ---- HDG boundary residual (FbouHdg in pdemodel.txt) ----------
    //
    // FbouHdg = τ*(0 - uhat) = -τ * uhat   (Dirichlet via stabilization)
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[],  const double /*uq*/[],
                  const double /*w*/[],  const double uh[],
                  const double /*n*/[],  const double tau[],
                  const double /*mu*/[], const double /*uinf*/[],
                  double /*t*/) {
        const double tau0  = tau[0];
        const double uhat0 = uh[0];
        fb[0] = -tau0 * uhat0;
    }

    // ∂fbou_hdg/∂uq = 0 (Poisson's HDG residual doesn't depend on uq).
    // ∂fbou_hdg/∂uh = -τ (1×1 block).
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double tau[],
                         const double /*mu*/[], const double /*uinf*/[],
                         double /*t*/) {
        fb_uh[0] = -tau[0];
    }

    // fbou_hdg_jac_uq, fbou_hdg_jac_w: zero — defaults from ModelDefaults.

    // ---- LDG boundary flux Fbou (used by some assembly paths) -----
    //
    // Fbou = f·n + τ*(uq[0] - uhat[0])  — exact form from pdemodel.txt
    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/,
              const double x[],  const double uq[],
              const double w[],  const double uh[],
              const double n[],  const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }

    // ---- Ubou: zero Dirichlet -------------------------------------
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[],
              double /*t*/) {
        ub[0] = 0.0;
    }

    // ---- Visualization & QoI (cheap; rarely on the hot path) ------

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(double s[], const double /*x*/[], const double uq[],
                     const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[0];
        s[1] = uq[1] + uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(double s[], const double /*x*/[], const double uq[],
                     const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[1];
        s[1] = uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(double s[], const double x[], const double uq[],
                    const double /*w*/[], const double /*mu*/[],
                    const double /*uinf*/[], double /*t*/) {
        // s[0] = (uq[0] - u_exact)^2  with u_exact = sin(πx) sin(πy)
        // s[1] = uq[0]
        const double t1 = 3.141592653589793;
        const double t2 = Kokkos::sin(t1 * x[0]);
        const double t3 = Kokkos::sin(t1 * x[1]);
        const double uexact = t2 * t3;
        s[0] = (uq[0] - uexact) * (uq[0] - uexact);
        s[1] = uq[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(double fb[], int /*ib*/,
                      const double x[],  const double uq[],
                      const double w[],  const double uh[],
                      const double n[],  const double tau[],
                      const double mu[], const double uinf[],
                      double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }
};
