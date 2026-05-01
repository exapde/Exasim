// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/visualization.hpp> — visualization output kernels.
//
// Replaces KokkosVis{Scalars,Vectors,Tensors} from libpdemodel.hpp.
// The user-supplied vis methods produce outputs of model-specific
// length (e.g. nscalars, nvectors*nd, ntensors*nd*nd). Output sizes
// are inferred from the runtime `nc` parameter the kernel receives.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

namespace detail {

// Iterate the user method `m(out, x, uq, v, w, mu, uinf, t)` over all
// `ng` quadrature points; output buffer sized `nc_out` per point.
template <class M, class Method>
void vis_dispatch(dstype* f, const dstype* xdg, const dstype* udg,
                  const dstype* odg, const dstype* wdg,
                  const dstype* param, dstype t,
                  int ng, int nc_out, const char* label, Method&& m)
{
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for(label, ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if constexpr (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        // Cap the output stack buffer at a generous max (16). Models
        // emitting more than 16 vis components per category should
        // route through a non-stack path; rare in practice.
        constexpr int kMax = 16;
        double out_local[kMax];
        m(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);

        for (int k = 0; k < nc_out; ++k) f[k * ng + i] = out_local[k];
    });
}

} // namespace detail

template <class M>
void vis_scalars_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        const dstype* /*uinf*/, const dstype* param, dstype t,
                        int /*modelnumber*/, int ng,
                        int nc_runtime, int /*ncu*/, int /*nd*/,
                        int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    detail::vis_dispatch<M>(f, xdg, udg, odg, wdg, param, t, ng, nc_runtime,
        "exasim::vis_scalars_kernel",
        [](double out[], const double x[], const double uq[],
           const double v[], const double w[],
           const double mu[], const double uinf[], double tt) {
            M::vis_scalars(out, x, uq, v, w, mu, uinf, tt);
        });
}

template <class M>
void vis_vectors_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        const dstype* /*uinf*/, const dstype* param, dstype t,
                        int /*modelnumber*/, int ng,
                        int nc_runtime, int /*ncu*/, int /*nd*/,
                        int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    detail::vis_dispatch<M>(f, xdg, udg, odg, wdg, param, t, ng, nc_runtime,
        "exasim::vis_vectors_kernel",
        [](double out[], const double x[], const double uq[],
           const double v[], const double w[],
           const double mu[], const double uinf[], double tt) {
            M::vis_vectors(out, x, uq, v, w, mu, uinf, tt);
        });
}

template <class M>
void vis_tensors_kernel(dstype* f, const dstype* xdg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        const dstype* /*uinf*/, const dstype* param, dstype t,
                        int /*modelnumber*/, int ng,
                        int nc_runtime, int /*ncu*/, int /*nd*/,
                        int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    detail::vis_dispatch<M>(f, xdg, udg, odg, wdg, param, t, ng, nc_runtime,
        "exasim::vis_tensors_kernel",
        [](double out[], const double x[], const double uq[],
           const double v[], const double w[],
           const double mu[], const double uinf[], double tt) {
            M::vis_tensors(out, x, uq, v, w, mu, uinf, tt);
        });
}

} // namespace exasim
