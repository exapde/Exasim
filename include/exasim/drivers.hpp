// SPDX-License-Identifier: see LICENSE
//
// <exasim/drivers.hpp> — templated *Driver<M> wrappers.
//
// This header is the bridge between the FEM internals (which use
// `meshstruct`, `appstruct`, `commonstruct`, … to carry runtime
// shape) and the templated kernels in <exasim/kernels/*.hpp> (which
// take `M::ncu`, `M::nd`, `M::ncw` from the user's Model type at
// compile time).
//
// Each function here mirrors the corresponding non-templated driver
// in `backend/Model/{KokkosDrivers,ModelDrivers}.cpp` — same name,
// same arg list, plus a leading `<class M>` template parameter. The
// body extracts sizes from `commonstruct` (just like today's drivers)
// and calls `exasim::*_kernel<M>(...)` instead of the libpdemodel.hpp
// ABI symbol.
//
// The existing non-templated drivers in `backend/Model/` are
// untouched. This header lives alongside them; HOT.2 step 3.2 will
// switch the FEM internals (CDiscretization & friends) to the
// templated path, which then makes the libpdemodel.hpp ABI
// boundary obsolete.

#pragma once

#include <utility>   // std::forward (for Fint/Fext stubs)

#include <Kokkos_Core.hpp>

#include "common.h"
#include "model.hpp"

#include "kernels/flux.hpp"
#include "kernels/source.hpp"
#include "kernels/sourcew.hpp"
#include "kernels/boundary.hpp"
#include "kernels/interface.hpp"
#include "kernels/init.hpp"
#include "kernels/tdfunc.hpp"
#include "kernels/eos.hpp"
#include "kernels/visualization.hpp"
#include "kernels/qoi.hpp"
#include "kernels/output.hpp"

namespace exasim {

namespace detail {

// Convenience accessor for the per-call shape. Mirrors the inlined
// extraction at the top of every existing *Driver function.
struct DriverShape {
    Int nc;
    Int ncu;
    Int ncw;
    Int nco;
    Int ncx;
    Int nd;
    Int numPoints;
    dstype time;
    Int modelnumber;
};

inline DriverShape compute_shape(commonstruct& common, Int nge, Int e1, Int e2) {
    return DriverShape{
        common.nc, common.ncu, common.ncw, common.nco,
        common.ncx, common.nd,
        nge * (e2 - e1),
        common.time,
        common.modelnumber
    };
}

} // namespace detail

// ===== Volume drivers =====

// LDG flux: f only.
template <class M>
inline void FluxDriver(dstype* f, const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nge, Int e1, Int e2,
                       Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    flux_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                   s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                   s.nd, s.ncx, s.nco, s.ncw);
}

// HDG flux: f + ∂f/∂uq + ∂f/∂w.
template <class M>
inline void FluxDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                       const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nge, Int e1, Int e2,
                       Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    hdg_flux_kernel<M>(f, f_udg, f_wdg, xg, udg, odg, wdg,
                       app.uinf, app.physicsparam,
                       s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                       s.nd, s.ncx, s.nco, s.ncw);
}

// LDG source.
template <class M>
inline void SourceDriver(dstype* f, const dstype* xg, const dstype* udg,
                         const dstype* odg, const dstype* wdg,
                         meshstruct& /*mesh*/, masterstruct& /*master*/,
                         appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                         commonstruct& common, Int nge, Int e1, Int e2,
                         Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    source_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                     s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                     s.nd, s.ncx, s.nco, s.ncw);
}

// HDG source: f + ∂s/∂uq + ∂s/∂w.
template <class M>
inline void SourceDriver(dstype* f, dstype* f_udg, dstype* f_wdg,
                         const dstype* xg, const dstype* udg,
                         const dstype* odg, const dstype* wdg,
                         meshstruct& /*mesh*/, masterstruct& /*master*/,
                         appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                         commonstruct& common, Int nge, Int e1, Int e2,
                         Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    hdg_source_kernel<M>(f, f_udg, f_wdg, xg, udg, odg, wdg,
                         app.uinf, app.physicsparam,
                         s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                         s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void TdfuncDriver(dstype* f, const dstype* xg, const dstype* udg,
                         const dstype* odg, const dstype* wdg,
                         meshstruct& /*mesh*/, masterstruct& /*master*/,
                         appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                         commonstruct& common, Int nge, Int e1, Int e2,
                         Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    tdfunc_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                     s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                     s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void AvfieldDriver(dstype* f, const dstype* xg, const dstype* udg,
                          const dstype* odg, const dstype* wdg,
                          meshstruct& /*mesh*/, masterstruct& /*master*/,
                          appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                          commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    avfield_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                      common.time, common.modelnumber, numPoints,
                      common.nc, common.ncu, common.nd, common.ncx,
                      common.nco, common.ncw, common.nce, common.npe, common.ne);
}

template <class M>
inline void OutputDriver(dstype* f, const dstype* xg, const dstype* udg,
                         const dstype* odg, const dstype* wdg,
                         meshstruct& /*mesh*/, masterstruct& /*master*/,
                         appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                         commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    output_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                     common.time, common.modelnumber, numPoints,
                     common.nc, common.ncu, common.nd, common.ncx,
                     common.nco, common.ncw, common.nce, common.npe, common.ne);
}

template <class M>
inline void MonitorDriver(dstype* f, Int nc_sol, const dstype* xg, const dstype* udg,
                          const dstype* odg, const dstype* wdg,
                          meshstruct& /*mesh*/, masterstruct& /*master*/,
                          appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                          commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    monitor_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                      common.time, common.modelnumber, numPoints,
                      nc_sol, common.ncu, common.nd, common.ncx,
                      common.nco, common.ncw, common.nce, common.npe, common.ne);
}

// ===== EoS drivers =====

template <class M>
inline void EosDriver(dstype* f, const dstype* xg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      meshstruct& /*mesh*/, masterstruct& /*master*/,
                      appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                      commonstruct& common, Int npe, Int e1, Int e2, Int /*backend*/)
{
    Int numPoints = npe * (e2 - e1);
    eos_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                  common.time, common.modelnumber, numPoints,
                  common.nc, common.ncu, common.nd, common.ncx,
                  common.nco, common.ncw, common.nce, npe, e2 - e1);
}

template <class M>
inline void EosduDriver(dstype* f, const dstype* xg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        meshstruct& /*mesh*/, masterstruct& /*master*/,
                        appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                        commonstruct& common, Int npe, Int e1, Int e2, Int /*backend*/)
{
    Int numPoints = npe * (e2 - e1);
    eos_du_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                     common.time, common.modelnumber, numPoints,
                     common.nc, common.ncu, common.nd, common.ncx,
                     common.nco, common.ncw, common.nce, npe, e2 - e1);
}

template <class M>
inline void EosdwDriver(dstype* f, const dstype* xg, const dstype* udg,
                        const dstype* odg, const dstype* wdg,
                        meshstruct& /*mesh*/, masterstruct& /*master*/,
                        appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                        commonstruct& common, Int npe, Int e1, Int e2, Int /*backend*/)
{
    Int numPoints = npe * (e2 - e1);
    eos_dw_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                     common.time, common.modelnumber, numPoints,
                     common.nc, common.ncu, common.nd, common.ncx,
                     common.nco, common.ncw, common.nce, npe, e2 - e1);
}

// ===== Sourcew drivers (auxiliary `w` field) =====

template <class M>
inline void SourcewDriver(dstype* f, const dstype* xg, const dstype* udg,
                          const dstype* odg, const dstype* wdg,
                          meshstruct& /*mesh*/, masterstruct& /*master*/,
                          appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                          commonstruct& common, Int npe, Int e1, Int e2, Int /*backend*/)
{
    Int numPoints = npe * (e2 - e1);
    sourcew_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                      common.time, common.modelnumber, numPoints,
                      common.nc, common.ncu, common.nd, common.ncx,
                      common.nco, common.ncw, common.nce, npe, e2 - e1);
}

// ===== Visualization & QoI drivers =====

template <class M>
inline void VisScalarsDriver(dstype* f, const dstype* xg, const dstype* udg,
                             const dstype* odg, const dstype* wdg,
                             meshstruct& /*mesh*/, masterstruct& /*master*/,
                             appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                             commonstruct& common, Int nge, Int e1, Int e2,
                             Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    vis_scalars_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                          s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                          s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void VisVectorsDriver(dstype* f, const dstype* xg, const dstype* udg,
                             const dstype* odg, const dstype* wdg,
                             meshstruct& /*mesh*/, masterstruct& /*master*/,
                             appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                             commonstruct& common, Int nge, Int e1, Int e2,
                             Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    vis_vectors_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                          s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                          s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void VisTensorsDriver(dstype* f, const dstype* xg, const dstype* udg,
                             const dstype* odg, const dstype* wdg,
                             meshstruct& /*mesh*/, masterstruct& /*master*/,
                             appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                             commonstruct& common, Int nge, Int e1, Int e2,
                             Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    vis_tensors_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                          s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                          s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void QoIvolumeDriver(dstype* f, const dstype* xg, const dstype* udg,
                            const dstype* odg, const dstype* wdg,
                            meshstruct& /*mesh*/, masterstruct& /*master*/,
                            appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                            commonstruct& common, Int nge, Int e1, Int e2,
                            Int /*backend*/)
{
    auto s = detail::compute_shape(common, nge, e1, e2);
    qoi_volume_kernel<M>(f, xg, udg, odg, wdg, app.uinf, app.physicsparam,
                         s.time, s.modelnumber, s.numPoints, s.nc, s.ncu,
                         s.nd, s.ncx, s.nco, s.ncw);
}

template <class M>
inline void QoIboundaryDriver(dstype* fb, const dstype* xg, const dstype* udg,
                              const dstype* odg, const dstype* wdg,
                              const dstype* uhg, const dstype* nl,
                              meshstruct& /*mesh*/, masterstruct& /*master*/,
                              appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                              commonstruct& common, Int ngf, Int f1, Int f2, Int ib, Int /*backend*/)
{
    Int ng = ngf * (f2 - f1);
    qoi_boundary_kernel<M>(fb, xg, udg, odg, wdg, uhg, nl, app.tau,
                           app.uinf, app.physicsparam,
                           common.time, common.modelnumber, ib, ng,
                           common.nc, common.ncu, common.nd, common.ncx,
                           common.nco, common.ncw);
}

// ===== Surface drivers =====

// LDG boundary flux (no Jacobians — `KokkosFbou`).
template <class M>
inline void FbouDriver(dstype* fb, const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int ngf, Int f1, Int f2, Int ib, Int /*backend*/)
{
    Int ng = ngf * (f2 - f1);
    fbou_kernel<M>(fb, xg, udg, odg, wdg, uhg, nl, app.tau,
                   app.uinf, app.physicsparam,
                   common.time, common.modelnumber, ib, ng,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

// HDG boundary residual + 3 Jacobians — calls M::fbou_hdg, NOT M::fbou.
// Single-block form takes `nga, ib, backend`.
template <class M>
inline void FbouDriver(dstype* fb, dstype* fb_udg, dstype* fb_wdg, dstype* fb_uhg,
                       const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nga, Int ib, Int /*backend*/)
{
    hdg_fbou_kernel<M>(fb, fb_udg, fb_wdg, fb_uhg, xg, udg, odg, wdg,
                       uhg, nl, app.tau, app.uinf, app.physicsparam,
                       common.time, common.modelnumber, ib, nga,
                       common.nc, common.ncu, common.nd, common.ncx,
                       common.nco, common.ncw);
}

template <class M>
inline void UbouDriver(dstype* ub, const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int ngf, Int f1, Int f2, Int ib, Int /*backend*/)
{
    Int ng = ngf * (f2 - f1);
    ubou_kernel<M>(ub, xg, udg, odg, wdg, uhg, nl, app.tau,
                   app.uinf, app.physicsparam,
                   common.time, common.modelnumber, ib, ng,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

// FhatDriver — LDG path (interior faces): `ngf, f1, f2, backend`.
template <class M>
inline void FhatDriver(dstype* fg, const dstype* xg,
                       const dstype* ug1, const dstype* ug2,
                       const dstype* og1, const dstype* og2,
                       const dstype* wg1, const dstype* wg2,
                       const dstype* uh, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int ngf, Int f1, Int f2, Int /*backend*/)
{
    Int ng = ngf * (f2 - f1);
    fhat_kernel<M>(fg, xg, ug1, ug2, og1, og2, wg1, wg2,
                   uh, nl, app.tau, app.uinf, app.physicsparam,
                   common.time, common.modelnumber, ng,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

// FhatDriver — single-block form (HDG path): `nga, backend`.
template <class M>
inline void FhatDriver(dstype* fg, const dstype* xg,
                       const dstype* ug1, const dstype* ug2,
                       const dstype* og1, const dstype* og2,
                       const dstype* wg1, const dstype* wg2,
                       const dstype* uh, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nga, Int /*backend*/)
{
    fhat_kernel<M>(fg, xg, ug1, ug2, og1, og2, wg1, wg2,
                   uh, nl, app.tau, app.uinf, app.physicsparam,
                   common.time, common.modelnumber, nga,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

// FhatDriver — HDG single-block, value + 3 Jacobians. Mirrors the
// legacy `FhatDriver(f, f_udg, f_wdg, f_uhg, ...)` in
// `backend/Model/ModelDrivers.cpp:558`. The HDG numerical flux is
// auto-derived from the user's Flux:
//
//   Fhat = Flux(uhat, q) · n + tau * (u - uhat)
//
// The Jacobians are derived consistently. This is the standard HDG
// numerical flux — the user does NOT define a separate `M::fhat`
// for it; instead it's built from `M::flux` + `M::flux_jac_uq`
// (which are already required) plus a fixed stabilization term.
template <class M>
inline void FhatDriver(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                       const dstype* xg, dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nga, Int /*backend*/)
{
    Int nc = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd  = common.nd;
    Int numPoints = nga;
    Int Mn = numPoints * ncu;
    Int N  = numPoints * ncu * nd;
    dstype time = common.time;

    // Stash u-component of udg into f_uhg; replace u with uhg
    // (so udg = (uh, q) for the flux call).
    ArrayCopy(f_uhg, udg, numPoints * ncu);
    ArrayCopy(udg,   uhg, numPoints * ncu);

    // f, f_udg, f_wdg via templated HDG flux kernel at (uh, q).
    hdg_flux_kernel<M>(f, f_udg, f_wdg, xg, udg, odg, wdg,
                       app.uinf, app.physicsparam, time,
                       common.modelnumber, numPoints, nc, ncu, nd,
                       ncx, nco, ncw);

    // Restore u-component of udg.
    ArrayCopy(udg, f_uhg, numPoints * ncu);

    // Dot-normal: f.n, f_udg.n (per nc-block), f_wdg.n (per ncw-block).
    FluxDotNormal(f, f, nl, Mn, numPoints, nd);
    for (int n = 0; n < nc; n++) {
        FluxDotNormal(&f_udg[Mn * n], &f_udg[N * n], nl, Mn, numPoints, nd);
    }
    if ((ncw > 0) & (common.wave == 0)) {
        for (int n = 0; n < ncw; n++) {
            FluxDotNormal(&f_wdg[Mn * n], &f_wdg[N * n], nl, Mn, numPoints, nd);
        }
    }

    // Move u-block of f_udg (∂(f.n)/∂u) into f_uhg as starting point;
    // zero f_udg's u-block; then add stabilization: f += tau*(u - uh),
    // f_udg += tau (on u-block), f_uhg -= tau (on uh-block).
    ArrayCopy(f_uhg, f_udg, numPoints * ncu * ncu);
    ArraySetValue(f_udg, zero, numPoints * ncu * ncu);
    AddStabilization1(f, f_udg, f_uhg, udg, uhg, app.tau, Mn, numPoints);
}

// FhatDriver — HDG single-block, value-only, with `u` scratch.
// Pattern: FhatDriver(fh, u, xg, udg, odg, wdg, uhg, nlg, ..., nga, backend).
// Mirrors `backend/Model/ModelDrivers.cpp:614`. Computes the standard
// HDG numerical flux Fhat = Flux(uhat, q) · n + tau * (u - uhat).
// The `u` argument is workspace for stashing udg's u-component while
// the flux is evaluated at (uhat, q).
template <class M>
inline void FhatDriver(dstype* fh, dstype* u, const dstype* xg,
                       dstype* udg, const dstype* odg, const dstype* wdg,
                       const dstype* uhg, dstype* nlg,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nga, Int /*backend*/)
{
    Int nc  = common.nc;
    Int ncu = common.ncu;
    Int ncw = common.ncw;
    Int nco = common.nco;
    Int ncx = common.ncx;
    Int nd  = common.nd;
    Int numPoints = nga;
    Int Mn = numPoints * ncu;
    dstype time = common.time;

    // Stash u-component, swap in uhg.
    ArrayCopy(u,   udg, numPoints * ncu);
    ArrayCopy(udg, uhg, numPoints * ncu);

    // Flux at (uh, q), value-only.
    flux_kernel<M>(fh, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                   common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);

    // Restore u-component.
    ArrayCopy(udg, u, numPoints * ncu);

    // Dot-normal then add stabilization.
    FluxDotNormal(fh, fh, nlg, Mn, numPoints, nd);
    AddStabilization1(fh, udg, uhg, app.tau, Mn);
}

// FbouDriver — single-block form (HDG path), value-only: `nga, ib, backend`.
// Calls `M::fbou_hdg` (boundary condition like `-tau*uhat`), matching
// legacy `HdgFbouonly` from `backend/Model/HdgFbouonly.cpp`. NOT
// `M::fbou`, which is the LDG numerical flux used by the LDG path.
template <class M>
inline void FbouDriver(dstype* fb, const dstype* xg, const dstype* udg,
                       const dstype* odg, const dstype* wdg,
                       const dstype* uhg, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int nga, Int ib, Int /*backend*/)
{
    hdg_fbou_only_kernel<M>(fb, xg, udg, odg, wdg, uhg, nl, app.tau,
                            app.uinf, app.physicsparam,
                            common.time, common.modelnumber, ib, nga,
                            common.nc, common.ncu, common.nd, common.ncx,
                            common.nco, common.ncw);
}

template <class M>
inline void UhatDriver(dstype* fg, const dstype* xg,
                       const dstype* ug1, const dstype* ug2,
                       const dstype* og1, const dstype* og2,
                       const dstype* wg1, const dstype* wg2,
                       const dstype* uh, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int ng, Int /*backend*/)
{
    uhat_kernel<M>(fg, xg, ug1, ug2, og1, og2, wg1, wg2,
                   uh, nl, app.tau, app.uinf, app.physicsparam,
                   common.time, common.modelnumber, ng,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

template <class M>
inline void StabDriver(dstype* fg, const dstype* xg,
                       const dstype* ug1, const dstype* ug2,
                       const dstype* og1, const dstype* og2,
                       const dstype* wg1, const dstype* wg2,
                       const dstype* uh, const dstype* nl,
                       meshstruct& /*mesh*/, masterstruct& /*master*/,
                       appstruct& app, solstruct& /*sol*/, tempstruct& /*temp*/,
                       commonstruct& common, Int ng, Int /*backend*/)
{
    stab_kernel<M>(fg, xg, ug1, ug2, og1, og2, wg1, wg2,
                   uh, nl, app.tau, app.uinf, app.physicsparam,
                   common.time, common.modelnumber, ng,
                   common.nc, common.ncu, common.nd, common.ncx,
                   common.nco, common.ncw);
}

// ===== Initialization drivers =====

template <class M>
inline void InitodgDriver(dstype* f, const dstype* xg,
                          appstruct& app, commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    initodg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      common.modelnumber, numPoints,
                      common.ncx, common.nce, common.npe, common.ne);
}

template <class M>
inline void InitqDriver(dstype* f, const dstype* xg,
                        appstruct& app, commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    initq_kernel<M>(f, xg, app.uinf, app.physicsparam,
                    common.modelnumber, numPoints,
                    common.ncx, common.nce, common.npe, common.ne);
}

template <class M>
inline void InitudgDriver(dstype* f, const dstype* xg,
                          appstruct& app, commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    initudg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      common.modelnumber, numPoints,
                      common.ncx, common.nce, common.npe, common.ne);
}

template <class M>
inline void InituDriver(dstype* f, const dstype* xg,
                        appstruct& app, commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    initu_kernel<M>(f, xg, app.uinf, app.physicsparam,
                    common.modelnumber, numPoints,
                    common.ncx, common.nce, common.npe, common.ne);
}

template <class M>
inline void InitwdgDriver(dstype* f, const dstype* xg,
                          appstruct& app, commonstruct& common, Int /*backend*/)
{
    Int numPoints = common.npe * common.ne;
    initwdg_kernel<M>(f, xg, app.uinf, app.physicsparam,
                      common.modelnumber, numPoints,
                      common.ncx, common.nce, common.npe, common.ne);
}

} // namespace exasim — end of templated *Driver<M> wrappers

// ===== Fint / Fext (multi-domain coupling, deferred to v2) =====
//
// Variadic stubs that forward to whatever ::FintDriver / ::FextDriver
// overload exists at the call site (resolved at instantiation time
// once Model/ModelDrivers.cpp has provided definitions).
//
// The real templated kernels for multi-domain coupling are recorded
// in <exasim/model.hpp>'s "Deferred surface" notes; until then, the
// stubs keep the templating cascade compilable so single-domain
// models (which never reach the multi-domain code paths at runtime)
// still build and run unchanged.

// Forward-declare the global Fint/Fext drivers so the variadic stubs
// below can reference them. Declaration order: callers in
// backend/Discretization/*.hpp will see them via the discretization.hpp
// chain (Model/ModelDrivers.cpp included before residual.hpp etc).
// We skip declaring exact signatures here — the actual overloads
// (LDG-only, HDG with Jacobians, etc.) live in ModelDrivers.cpp and
// resolve via overload resolution at the call site.

namespace exasim::detail {

// Trick: the body of these forward stubs uses unqualified
// FintDriver / FextDriver and forwards via a helper that is
// instantiation-dependent on M. This delays name lookup until
// instantiation time, by which point Model/ModelDrivers.cpp has
// provided the global definitions.
template <class M>
struct multidomain_forward {
    template <class... Args>
    static void Fint(Args&&... args) {
        FintDriver(std::forward<Args>(args)...);   // ADL + unqualified -> picks ::FintDriver at instantiation
    }
    template <class... Args>
    static void Fext(Args&&... args) {
        FextDriver(std::forward<Args>(args)...);
    }
};

} // namespace exasim::detail

namespace exasim {

template <class M, class... Args>
inline void FintDriver(Args&&... args) {
    ::exasim::detail::multidomain_forward<M>::Fint(std::forward<Args>(args)...);
}

template <class M, class... Args>
inline void FextDriver(Args&&... args) {
    ::exasim::detail::multidomain_forward<M>::Fext(std::forward<Args>(args)...);
}

} // namespace exasim — end of Fint/Fext stubs
