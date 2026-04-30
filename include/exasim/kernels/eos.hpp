// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/eos.hpp> — equation-of-state and AV field kernels.
//
// Replaces KokkosEoS, KokkosEoSdu, KokkosEoSdw, KokkosAvfield, and the
// HDG variant HdgEoS from libpdemodel.hpp. EoS pointwise functions
// produce the same shape as source (ncu output); the *_du / *_dw
// flavors emit Jacobian blocks.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

namespace detail {

template <class M, int OutputSize, class Method>
KOKKOS_INLINE_FUNCTION
void eos_dispatch(dstype* f, size_t i, int ng,
                  const dstype* xdg, const dstype* udg, const dstype* wdg,
                  const dstype* param, dstype t, Method&& m)
{
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    double x[nd], uq[Nq], w[ncw_buf];
    for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
    for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
    if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

    double out_local[OutputSize];
    m(out_local, x, uq, w, param, t);
    for (int k = 0; k < OutputSize; ++k) f[k * ng + i] = out_local[k];
}

} // namespace detail

template <class M>
void eos_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    Kokkos::parallel_for("exasim::eos_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        detail::eos_dispatch<M, M::ncu>(f, i, ng, xdg, udg, wdg, param, t,
            [](double out[], const double x[], const double uq[], const double w[],
               const double mu[], double tt) {
                M::eos(out, x, uq, w, mu, /*uinf=*/nullptr, tt);
            });
    });
}

template <class M>
void eos_du_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                   const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                   int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                   int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int Nq = M::ncu * (1 + M::nd);
    Kokkos::parallel_for("exasim::eos_du_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        detail::eos_dispatch<M, M::ncu * Nq>(f, i, ng, xdg, udg, wdg, param, t,
            [](double out[], const double x[], const double uq[], const double w[],
               const double mu[], double tt) {
                M::eos_du(out, x, uq, w, mu, /*uinf=*/nullptr, tt);
            });
    });
}

template <class M>
void eos_dw_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                   const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                   int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                   int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    if constexpr (M::ncw > 0) {
        Kokkos::parallel_for("exasim::eos_dw_kernel", ng, KOKKOS_LAMBDA(size_t i) {
            detail::eos_dispatch<M, M::ncu * M::ncw>(f, i, ng, xdg, udg, wdg, param, t,
                [](double out[], const double x[], const double uq[], const double w[],
                   const double mu[], double tt) {
                    M::eos_dw(out, x, uq, w, mu, /*uinf=*/nullptr, tt);
                });
        });
    } else {
        (void)f; (void)xdg; (void)udg; (void)wdg; (void)param; (void)t; (void)ng;
    }
}

template <class M>
void avfield_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                    const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                    int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                    int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    Kokkos::parallel_for("exasim::avfield_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        detail::eos_dispatch<M, M::ncu>(f, i, ng, xdg, udg, wdg, param, t,
            [](double out[], const double x[], const double uq[], const double w[],
               const double mu[], double tt) {
                M::avfield(out, x, uq, w, mu, /*uinf=*/nullptr, tt);
            });
    });
}

} // namespace exasim
