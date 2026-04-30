// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/source.hpp> — templated source-term kernels.
//
// Replaces `KokkosSource` and `HdgSource` from libpdemodel.hpp.
//
// The forcing/source pointwise function emits an `ncu`-vector. The
// HDG variant additionally fills `f_udg = ∂s/∂uq` and (when ncw > 0)
// `f_wdg = ∂s/∂w` from the user's hand-written `M::source_jac_uq`
// and `M::source_jac_w` methods.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

// Forward path — value only.
template <class M>
void source_kernel(dstype*       s,
                   const dstype* xdg,
                   const dstype* udg,
                   const dstype* /*odg*/,
                   const dstype* wdg,
                   const dstype* /*uinf*/,
                   const dstype* param,
                   dstype        t,
                   int           /*modelnumber*/,
                   int           ng,
                   int           /*nc*/,
                   int           ncu_runtime,
                   int           nd_runtime,
                   int           /*ncx*/,
                   int           /*nco*/,
                   int           ncw_runtime)
{
    static_assert(is_model_v<M>, "source_kernel<M>: M must satisfy the Model contract.");

    constexpr int nd  = M::nd;
    constexpr int ncu = M::ncu;
    constexpr int ncw = M::ncw;
    constexpr int Nq  = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    assert(ncu_runtime == ncu && nd_runtime == nd && ncw_runtime == ncw);
    (void)ncu_runtime; (void)nd_runtime; (void)ncw_runtime;

    Kokkos::parallel_for("exasim::source_kernel", ng,
        KOKKOS_LAMBDA(const size_t i) {
            double x [nd];
            double uq[Nq];
            double w [ncw_buf];

            for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
            for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
            if constexpr (ncw > 0) {
                for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
            }

            double s_local[ncu];
            M::source(s_local, x, uq, w, param, t);

            for (int k = 0; k < ncu; ++k) s[k * ng + i] = s_local[k];
        });
}

// HDG path — value + ∂s/∂uq + ∂s/∂w.
template <class M>
void hdg_source_kernel(dstype*       s,
                       dstype*       s_udg,
                       dstype*       s_wdg,
                       const dstype* xdg,
                       const dstype* udg,
                       const dstype* /*odg*/,
                       const dstype* wdg,
                       const dstype* /*uinf*/,
                       const dstype* param,
                       dstype        t,
                       int           /*modelnumber*/,
                       int           ng,
                       int           /*nc*/,
                       int           ncu_runtime,
                       int           nd_runtime,
                       int           /*ncx*/,
                       int           /*nco*/,
                       int           ncw_runtime)
{
    static_assert(is_model_v<M>, "hdg_source_kernel<M>: M must satisfy the Model contract.");

    constexpr int nd  = M::nd;
    constexpr int ncu = M::ncu;
    constexpr int ncw = M::ncw;
    constexpr int Nq  = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    assert(ncu_runtime == ncu && nd_runtime == nd && ncw_runtime == ncw);
    (void)ncu_runtime; (void)nd_runtime; (void)ncw_runtime;

    Kokkos::parallel_for("exasim::hdg_source_kernel", ng,
        KOKKOS_LAMBDA(const size_t i) {
            double x [nd];
            double uq[Nq];
            double w [ncw_buf];

            for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
            for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
            if constexpr (ncw > 0) {
                for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
            }

            // Value
            double s_local[ncu];
            M::source(s_local, x, uq, w, param, t);
            for (int k = 0; k < ncu; ++k) s[k * ng + i] = s_local[k];

            // ∂s/∂uq
            double s_uq[ncu * Nq];
            M::source_jac_uq(s_uq, x, uq, w, param, t);
            for (int k = 0; k < ncu * Nq; ++k) s_udg[k * ng + i] = s_uq[k];

            // ∂s/∂w (only when present)
            if constexpr (ncw > 0) {
                double s_w[ncu * ncw];
                M::source_jac_w(s_w, x, uq, w, param, t);
                for (int k = 0; k < ncu * ncw; ++k) s_wdg[k * ng + i] = s_w[k];
            }
        });
}

} // namespace exasim
