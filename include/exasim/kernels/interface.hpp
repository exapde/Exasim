// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/interface.hpp> — interface kernels: fhat, uhat, stab.
//
// Internal-face kernels see both sides (`udg1`, `udg2`, `wdg1`, `wdg2`,
// `odg1`, `odg2`) plus the trace `uhg`, normal `nlg`, and stabilization
// `tau`. libpdemodel.hpp does not provide HDG variants for these
// (interface fluxes are evaluated only at the value level).

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void fhat_kernel(dstype* f, const dstype* xdg,
                 const dstype* udg1, const dstype* udg2,
                 const dstype* odg1, const dstype* odg2,
                 const dstype* wdg1,  const dstype* wdg2,
                 const dstype* uhg,   const dstype* nlg, const dstype* tau,
                 const dstype* /*uinf*/, const dstype* param,
                 dstype t, int /*modelnumber*/, int ng,
                 int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::fhat_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg1; (void)odg2; (void)wdg1; (void)wdg2;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq1[Nq], uq2[Nq];
        double v1[nco_buf], v2[nco_buf], w1[ncw_buf], w2[ncw_buf];
        double uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd; ++k)  x  [k] = xdg [k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq1[k] = udg1[k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq2[k] = udg2[k * ng + i];
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) v1[k] = odg1[k * ng + i];
            for (int k = 0; k < nco; ++k) v2[k] = odg2[k * ng + i];
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) w1[k] = wdg1[k * ng + i];
            for (int k = 0; k < ncw; ++k) w2[k] = wdg2[k * ng + i];
        }
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        double f_local[ncu];
        M::fhat(f_local, x, uq1, uq2, v1, v2, w1, w2, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = f_local[k];
    });
}

template <class M>
void uhat_kernel(dstype* f, const dstype* xdg,
                 const dstype* udg1, const dstype* udg2,
                 const dstype* odg1, const dstype* odg2,
                 const dstype* wdg1,  const dstype* wdg2,
                 const dstype* uhg,   const dstype* nlg, const dstype* tau,
                 const dstype* /*uinf*/, const dstype* param,
                 dstype t, int /*modelnumber*/, int ng,
                 int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::uhat_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg1; (void)odg2; (void)wdg1; (void)wdg2;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq1[Nq], uq2[Nq];
        double v1[nco_buf], v2[nco_buf], w1[ncw_buf], w2[ncw_buf];
        double uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd; ++k)  x  [k] = xdg [k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq1[k] = udg1[k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq2[k] = udg2[k * ng + i];
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) v1[k] = odg1[k * ng + i];
            for (int k = 0; k < nco; ++k) v2[k] = odg2[k * ng + i];
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) w1[k] = wdg1[k * ng + i];
            for (int k = 0; k < ncw; ++k) w2[k] = wdg2[k * ng + i];
        }
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        double f_local[ncu];
        M::uhat(f_local, x, uq1, uq2, v1, v2, w1, w2, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = f_local[k];
    });
}

template <class M>
void stab_kernel(dstype* f, const dstype* xdg,
                 const dstype* udg1, const dstype* udg2,
                 const dstype* odg1, const dstype* odg2,
                 const dstype* wdg1,  const dstype* wdg2,
                 const dstype* uhg,   const dstype* nlg, const dstype* tau,
                 const dstype* /*uinf*/, const dstype* param,
                 dstype t, int /*modelnumber*/, int ng,
                 int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/, int /*ncw*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::stab_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg1; (void)odg2; (void)wdg1; (void)wdg2;  // HOT.6.2 nvcc force-capture: see /tmp/patch_constexpr_capture.py
        double x[nd], uq1[Nq], uq2[Nq];
        double v1[nco_buf], v2[nco_buf], w1[ncw_buf], w2[ncw_buf];
        double uh[ncu], n[nd], t_[ncu];
        for (int k = 0; k < nd; ++k)  x  [k] = xdg [k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq1[k] = udg1[k * ng + i];
        for (int k = 0; k < Nq; ++k)  uq2[k] = udg2[k * ng + i];
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) v1[k] = odg1[k * ng + i];
            for (int k = 0; k < nco; ++k) v2[k] = odg2[k * ng + i];
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) w1[k] = wdg1[k * ng + i];
            for (int k = 0; k < ncw; ++k) w2[k] = wdg2[k * ng + i];
        }
        for (int k = 0; k < ncu; ++k) uh[k] = uhg[k * ng + i];
        for (int k = 0; k < nd;  ++k) n [k] = nlg[k * ng + i];
        for (int k = 0; k < ncu; ++k) t_[k] = tau[k];

        double f_local[ncu];
        M::stab(f_local, x, uq1, uq2, v1, v2, w1, w2, uh, n, t_, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = f_local[k];
    });
}

} // namespace exasim
