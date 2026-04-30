// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/sourcew.hpp> — source kernel for the auxiliary `w` field.
//
// Replaces KokkosSourcew / HdgSourcew. Output size is `ncw` (auxiliary
// field count) rather than `ncu`. Only meaningful when `ncw > 0`.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void sourcew_kernel(dstype* sw,
                    const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                    const dstype* wdg, const dstype* /*uinf*/, const dstype* param,
                    dstype t, int /*modelnumber*/, int ng,
                    int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/,
                    int /*ncw_runtime*/, int /*nce*/, int /*npe*/, int /*ne*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);

    if constexpr (ncw > 0) {
        Kokkos::parallel_for("exasim::sourcew_kernel", ng, KOKKOS_LAMBDA(size_t i) {
            double x[nd], uq[Nq], w[ncw];
            for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
            for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
            for (int k = 0; k < ncw; ++k) w [k] = wdg[k * ng + i];

            double sw_local[ncw];
            M::sourcew(sw_local, x, uq, w, param, /*uinf=*/nullptr, t);
            for (int k = 0; k < ncw; ++k) sw[k * ng + i] = sw_local[k];
        });
    } else {
        (void)sw; (void)xdg; (void)udg; (void)wdg; (void)param; (void)t; (void)ng;
    }
}

template <class M>
void hdg_sourcew_kernel(dstype* sw, dstype* sw_udg, dstype* sw_wdg,
                        const dstype* xdg, const dstype* udg, const dstype* /*odg*/,
                        const dstype* wdg, const dstype* /*uinf*/, const dstype* param,
                        dstype t, int /*modelnumber*/, int ng,
                        int /*nc*/, int /*ncu*/, int /*nd*/, int /*ncx*/, int /*nco*/,
                        int /*ncw_runtime*/)
{
    static_assert(is_model_v<M>);
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw;
    constexpr int Nq = ncu * (1 + nd);

    if constexpr (ncw > 0) {
        Kokkos::parallel_for("exasim::hdg_sourcew_kernel", ng, KOKKOS_LAMBDA(size_t i) {
            double x[nd], uq[Nq], w[ncw];
            for (int k = 0; k < nd;  ++k) x [k] = xdg[k * ng + i];
            for (int k = 0; k < Nq;  ++k) uq[k] = udg[k * ng + i];
            for (int k = 0; k < ncw; ++k) w [k] = wdg[k * ng + i];

            double sw_local[ncw];
            M::sourcew(sw_local, x, uq, w, param, /*uinf=*/nullptr, t);
            for (int k = 0; k < ncw; ++k) sw[k * ng + i] = sw_local[k];

            // ∂sw/∂uq, ∂sw/∂w — when text2code generates these the
            // user model's *_jac_* methods must be present.
            // Defaults in ModelDefaults zero-fill.
            // (Kept minimal here — model.hpp defines the contract.)
            (void)sw_udg; (void)sw_wdg;
        });
    } else {
        (void)sw; (void)sw_udg; (void)sw_wdg;
        (void)xdg; (void)udg; (void)wdg; (void)param; (void)t; (void)ng;
    }
}

} // namespace exasim
