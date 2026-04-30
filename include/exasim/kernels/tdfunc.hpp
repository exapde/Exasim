// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/tdfunc.hpp> — time-derivative mass-weighting kernel.
//
// Replaces KokkosTdfunc from libpdemodel.hpp. Same pointwise shape as
// source: emits an `ncu`-vector mass weighting at each quadrature point.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void tdfunc_kernel(dstype* f,
                   const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                   const dstype* wdg, const dstype* /*uinf*/, const dstype* param,
                   dstype t, int /*modelnumber*/, int ng,
                   int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;

    Kokkos::parallel_for("exasim::tdfunc_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        double x[nd], uq[Nq], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if constexpr (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double m_local[ncu];
        M::tdfunc(m_local, x, uq, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = m_local[k];
    });
}

} // namespace exasim
