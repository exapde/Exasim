// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/eos.hpp> — equation-of-state and AV field kernels.
//
// Replaces KokkosEoS, KokkosEoSdu, KokkosEoSdw, KokkosAvfield, and the
// HDG variant HdgEoS from libpdemodel.hpp. EoS pointwise functions
// produce the same shape as source (ncu output); the *_du / *_dw
// flavors emit Jacobian blocks.
//
// HOT.6.3: dispatch helper inlined into each kernel — see
// `output.hpp` for the same nvcc-extended-lambda restriction.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void eos_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::eos_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[ncu];
        M::eos(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = out_local[k];
    });
}

template <class M>
void eos_du_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                   int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                   int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::eos_du_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[ncu * Nq];
        M::eos_du(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu * Nq; ++k) f[k * ng + i] = out_local[k];
    });
}

template <class M>
void eos_dw_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                   const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                   int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                   int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    if constexpr (M::ncw > 0) {
        constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
        constexpr int Nq = ncu * (1 + nd);
        constexpr int nco_buf = (nco > 0) ? nco : 1;

        Kokkos::parallel_for("exasim::eos_dw_kernel", ng, KOKKOS_LAMBDA(size_t i) {
            (void)odg; (void)wdg;
            double x[nd], uq[Nq], v[nco_buf], w[ncw];
            for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
            for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
            if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
            for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

            double out_local[ncu * ncw];
            M::eos_dw(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
            for (int k = 0; k < ncu * ncw; ++k) f[k * ng + i] = out_local[k];
        });
    } else {
        (void)f; (void)xdg; (void)udg; (void)odg; (void)wdg; (void)param; (void)t; (void)ng;
    }
}

template <class M>
void avfield_kernel(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg,
                    const dstype* wdg, const dstype* /*uinf*/, const dstype* param, dstype t,
                    int /*modelnumber*/, int ng, int /*nc*/, int /*ncu*/, int /*nd*/,
                    int /*ncx*/, int /*nco*/, int /*ncw*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;

    Kokkos::parallel_for("exasim::avfield_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[ncu];
        M::avfield(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = out_local[k];
    });
}

} // namespace exasim
