// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/output.hpp> — monitor / output diagnostic kernels.
//
// Replaces KokkosMonitor / KokkosOutput from libpdemodel.hpp. Same
// pointwise shape as `source` (volume integrand); output count is
// supplied by the runtime.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

namespace detail {

template <class M, class Method>
void output_dispatch(dstype* f, const dstype* xdg, const dstype* udg,
                     const dstype* odg, const dstype* wdg,
                     const dstype* param, dstype t,
                     int ng, int nc_out, const char* label, Method&& m)
{
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;
    constexpr int kMax = 16;

    Kokkos::parallel_for(label, ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if constexpr (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[kMax];
        m(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);

        for (int k = 0; k < nc_out; ++k) f[k * ng + i] = out_local[k];
    });
}

} // namespace detail

template <class M>
void monitor_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                    const dstype* odg, const dstype* wdg,
                    const dstype* /*uinf*/, const dstype* param, dstype t,
                    int /*modelnumber*/, int ng,
                    int nc_runtime, int /*ncu*/, int /*nd*/,
                    int /*ncx*/, int /*nco*/, int /*ncw*/,
                    int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    detail::output_dispatch<M>(f, xdg, udg, odg, wdg, param, t, ng, nc_runtime,
        "exasim::monitor_kernel",
        [](double out[], const double x[], const double uq[],
           const double v[], const double w[],
           const double mu[], const double uinf[], double tt) {
            M::monitor(out, x, uq, v, w, mu, uinf, tt);
        });
}

template <class M>
void output_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                   const dstype* odg, const dstype* wdg,
                   const dstype* /*uinf*/, const dstype* param, dstype t,
                   int /*modelnumber*/, int ng,
                   int nc_runtime, int /*ncu*/, int /*nd*/,
                   int /*ncx*/, int /*nco*/, int /*ncw*/,
                   int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    detail::output_dispatch<M>(f, xdg, udg, odg, wdg, param, t, ng, nc_runtime,
        "exasim::output_kernel",
        [](double out[], const double x[], const double uq[],
           const double v[], const double w[],
           const double mu[], const double uinf[], double tt) {
            M::output(out, x, uq, v, w, mu, uinf, tt);
        });
}

} // namespace exasim
