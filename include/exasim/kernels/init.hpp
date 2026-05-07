// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/init.hpp> — initial-condition kernels.
//
// Replaces KokkosInit{u,q,udg,wdg,odg} from libpdemodel.hpp. The init
// kernels iterate over total nodes (`ng = npe * ne`) rather than
// quadrature points; the SoA layout is element-major:
//   f[node + npe*comp + npe*nc*elem]

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

namespace detail {

// Common scaffolding for every init kernel. The user method `Method`
// takes (out, x, uinf, mu) and fills `OutputSize` components.
template <class M, int OutputSize, class Method>
KOKKOS_INLINE_FUNCTION
void init_dispatch(dstype* f, const dstype* xdg, const dstype* uinf,
                   const dstype* param, int ng, int npe, Method&& m)
{
    constexpr int nd = M::nd;
    Kokkos::parallel_for("exasim::init_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        const int j    = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);

        double x[nd];
        for (int k = 0; k < nd; ++k) {
            // xdg layout: [npe x ncx x ne], column-major; x at this node:
            x[k] = xdg[j + npe * k + npe * nd * elem];
        }

        double out_local[OutputSize];
        m(out_local, x, uinf, param);

        for (int k = 0; k < OutputSize; ++k) {
            f[j + npe * k + npe * OutputSize * elem] = out_local[k];
        }
    });
}

} // namespace detail

template <class M>
void initu_kernel(dstype* f, const dstype* xdg, const dstype* uinf,
                  const dstype* param, int /*modelnumber*/, int ng,
                  int /*ncx*/, int /*nce*/, int npe, int /*ne*/)
{
    static_assert(is_model_v<M>);
    detail::init_dispatch<M, M::ncu>(
        f, xdg, uinf, param, ng, npe,
        [](double out[], const double x[], const double uinf[], const double mu[])
 { M::initu(out, x, uinf, mu); });
}

template <class M>
void initq_kernel(dstype* f, const dstype* xdg, const dstype* uinf,
                  const dstype* param, int /*modelnumber*/, int ng,
                  int /*ncx*/, int /*nce*/, int npe, int /*ne*/)
{
    static_assert(is_model_v<M>);
    detail::init_dispatch<M, M::ncu * M::nd>(
        f, xdg, uinf, param, ng, npe,
        [](double out[], const double x[], const double uinf[], const double mu[])
 { M::initq(out, x, uinf, mu); });
}

template <class M>
void initudg_kernel(dstype* f, const dstype* xdg, const dstype* uinf,
                    const dstype* param, int /*modelnumber*/, int ng,
                    int /*ncx*/, int /*nce*/, int npe, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int Nq = M::ncu * (1 + M::nd);
    detail::init_dispatch<M, Nq>(
        f, xdg, uinf, param, ng, npe,
        [](double out[], const double x[], const double uinf[], const double mu[])
 { M::initudg(out, x, uinf, mu); });
}

template <class M>
void initwdg_kernel(dstype* f, const dstype* xdg, const dstype* uinf,
                    const dstype* param, int /*modelnumber*/, int ng,
                    int /*ncx*/, int /*nce*/, int npe, int /*ne*/)
{
    static_assert(is_model_v<M>);
    if constexpr (M::ncw > 0) {
        detail::init_dispatch<M, M::ncw>(
            f, xdg, uinf, param, ng, npe,
            [](double out[], const double x[], const double uinf[], const double mu[])
 { M::initwdg(out, x, uinf, mu); });
    } else {
        (void)f; (void)xdg; (void)uinf; (void)param; (void)ng; (void)npe;
    }
}

// initodg's output count is determined by the discretization (`nco`),
// not by Self. The runtime nce arg carries the count.
template <class M>
void initodg_kernel(dstype* f, const dstype* xdg, const dstype* uinf,
                    const dstype* param, int /*modelnumber*/, int ng,
                    int /*ncx*/, int nce, int npe, int /*ne*/)
{
    static_assert(is_model_v<M>);
    Kokkos::parallel_for("exasim::initodg_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        const int j    = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);
        for (int k = 0; k < nce; ++k) {
            f[j + npe * k + npe * nce * elem] = 0.0;
        }
        // Models that need a non-zero odg field can override M::initodg;
        // we don't pre-pull the user override here because the size isn't
        // known at compile time.
        (void)xdg; (void)uinf; (void)param;
    });
}

} // namespace exasim
