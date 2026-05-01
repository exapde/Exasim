// SPDX-License-Identifier: see LICENSE
//
// <exasim/kernels/visualization.hpp> — visualization output kernels.
//
// Replaces KokkosVis{Scalars,Vectors,Tensors} from libpdemodel.hpp.
// The user-supplied vis methods produce outputs of model-specific
// length (e.g. nscalars, nvectors*nd, ntensors*nd*nd). Output sizes
// are inferred from the runtime `nc` parameter the kernel receives.
//
// HOT.6.3: dispatch helper inlined into each kernel — see
// `output.hpp` for the same nvcc-extended-lambda restriction.

#pragma once

#include <Kokkos_Core.hpp>

#include "../common.h"
#include "../model.hpp"

namespace exasim {

template <class M>
void vis_scalars_kernel(dstype* f, const dstype* xdg, const dstype* udg,
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

    Kokkos::parallel_for("exasim::vis_scalars_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[kMax];
        M::vis_scalars(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < nc_runtime; ++k) f[k * ng + i] = out_local[k];
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
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;
    constexpr int kMax = 16;

    Kokkos::parallel_for("exasim::vis_vectors_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[kMax];
        M::vis_vectors(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < nc_runtime; ++k) f[k * ng + i] = out_local[k];
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
    constexpr int nd = M::nd, ncu = M::ncu, ncw = M::ncw, nco = M::nco;
    constexpr int Nq = ncu * (1 + nd);
    constexpr int ncw_buf = (ncw > 0) ? ncw : 1;
    constexpr int nco_buf = (nco > 0) ? nco : 1;
    constexpr int kMax = 16;

    Kokkos::parallel_for("exasim::vis_tensors_kernel", ng, KOKKOS_LAMBDA(size_t i) {
        (void)odg; (void)wdg;
        double x[nd], uq[Nq], v[nco_buf], w[ncw_buf];
        for (int k = 0; k < nd; ++k) x [k] = xdg[k * ng + i];
        for (int k = 0; k < Nq; ++k) uq[k] = udg[k * ng + i];
        if (nco > 0) for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        if (ncw > 0) for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        double out_local[kMax];
        M::vis_tensors(out_local, x, uq, v, w, param, /*uinf=*/nullptr, t);
        for (int k = 0; k < nc_runtime; ++k) f[k * ng + i] = out_local[k];
    });
}

} // namespace exasim
