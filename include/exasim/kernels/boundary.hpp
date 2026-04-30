// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/boundary.hpp> — boundary kernels: fbou, ubou + their HDG variants.
//
// Boundary kernels operate on faces tagged with `ib` and additionally
// see the trace `uhg`, the outward normal `nlg`, and the stabilization
// `tau`. The HDG variants fill four Jacobian outputs:
//   f_udg = ∂fb/∂uq    (ncu × Nq    × ng)
//   f_wdg = ∂fb/∂w     (ncu × ncw   × ng)  [only if ncw > 0]
//   f_uhg = ∂fb/∂uh    (ncu × ncu   × ng)

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void fbou_kernel(dstype*       fb,
                 const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg,
                 const dstype* tau, const dstype* /*uinf*/, const dstype* param,
                 dstype t, int /*modelnumber*/, int ib, int ng,
                 int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    Kokkos::parallel_for("exasim::fbou_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw_buf], uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];   // tau is per-component, not per-i

        double fb_local[ncu];
        M::fbou(fb_local, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) fb[k * ng + i] = fb_local[k];
    });
}

template <class M>
void hdg_fbou_kernel(dstype* fb, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                     const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                     const dstype* wdg, const dstype* uhg, const dstype* nlg,
                     const dstype* tau, const dstype* /*uinf*/, const dstype* param,
                     dstype t, int /*modelnumber*/, int ib, int ng,
                     int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    Kokkos::parallel_for("exasim::hdg_fbou_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw_buf], uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        // Value
        double fb_local[ncu];
        M::fbou(fb_local, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) fb[k * ng + i] = fb_local[k];

        // ∂fb/∂uq
        double fb_uq[ncu * Nq];
        M::fbou_jac_uq(fb_uq, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu * Nq; ++k) f_udg[k * ng + i] = fb_uq[k];

        // ∂fb/∂w
        if constexpr (ncw > 0) {
            double fb_w[ncu * ncw];
            M::fbou_jac_w(fb_w, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
            for (int k = 0; k < ncu * ncw; ++k) f_wdg[k * ng + i] = fb_w[k];
        }

        // ∂fb/∂uh
        double fb_uh[ncu * ncu];
        M::fbou_jac_uh(fb_uh, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu * ncu; ++k) f_uhg[k * ng + i] = fb_uh[k];
    });
}

template <class M>
void ubou_kernel(dstype* ub,
                 const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg,
                 const dstype* tau, const dstype* /*uinf*/, const dstype* param,
                 dstype t, int /*modelnumber*/, int ib, int ng,
                 int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    Kokkos::parallel_for("exasim::ubou_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw_buf], uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        double ub_local[ncu];
        M::ubou(ub_local, ib, x, uq, w, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) ub[k * ng + i] = ub_local[k];
    });
}

} // namespace exasim
