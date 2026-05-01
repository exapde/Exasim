// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/qoi.hpp> — Quantity-of-Interest integrand kernels.
//
// Replaces KokkosQoIvolume / KokkosQoIboundary from libpdemodel.hpp.
// QoI integrands are pointwise expressions evaluated at quadrature
// points; the runtime integrates them in element / face mass loops.
// Output count is model-defined (number of QoIs the user wants).

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void qoi_volume_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* /*uinf*/, const dstype* param, dstype t,
                       int /*modelnumber*/, int ng,
                       int nc_runtime, int /*ncu*/, int /*nd*/,
                       int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;
    constexpr int kMax = 16;

    Kokkos::parallel_for("exasim::qoi_volume_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[kMax];
        M::qoi_volume(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);

        for (int k = 0; k < nc_runtime; ++k) f[k * ng + i] = out_local[k];
    });
}

template <class M>
void qoi_boundary_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                         const dstype* odg, const dstype* wdg,
                         const dstype* uhg, const dstype* nlg, const dstype* tau,
                         const dstype* /*uinf*/, const dstype* param, dstype t,
                         int /*modelnumber*/, int ib, int ng,
                         int nc_runtime, int /*ncu*/, int /*nd*/,
                         int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;
    constexpr int kMax = 16;

    Kokkos::parallel_for("exasim::qoi_boundary_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf], uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        double out_local[kMax];
        M::qoi_boundary(out_local, ib, x, uq, v, w, uh, n, t_, param, /*uinf=*/nullptr, t);

        for (int k = 0; k < nc_runtime; ++k) f[k * ng + i] = out_local[k];
    });
}

} // namespace exasim
