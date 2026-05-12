/**
 * @file model_ops.hpp
 * @brief Runtime dispatch table for model-dependent operations in the hybrid Exasim architecture.
 *
 * This header defines ModelOps, a lightweight collection of host-callable
 * function pointers used by the non-templated backend to reach model-specific
 * kernels. The slot layout mirrors the 42 legacy *Driver definitions shared by
 * KokkosDrivers.cpp, ModelDrivers.cpp, and BuiltinModelDrivers.cpp so the
 * provider-specific dispatch adapters can bind those implementations with full
 * type checking.
 */
#ifndef __EXASIM_MODEL_DISPATCH_MODEL_OPS_HPP__
#define __EXASIM_MODEL_DISPATCH_MODEL_OPS_HPP__

#include "../../Common/common.h"

namespace exasim {

struct ModelOps {
    using ElementDriverFn =
        void (*)(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, meshstruct& mesh,
                 masterstruct& master, appstruct& app, solstruct& sol,
                 tempstruct& temp, commonstruct& common, Int nge, Int e1,
                 Int e2, Int backend);

    using GlobalElementDriverFn =
        void (*)(dstype* f, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, meshstruct& mesh,
                 masterstruct& master, appstruct& app, solstruct& sol,
                 tempstruct& temp, commonstruct& common, Int backend);

    using MonitorDriverFn =
        void (*)(dstype* f, Int nc_sol, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, meshstruct& mesh,
                 masterstruct& master, appstruct& app, solstruct& sol,
                 tempstruct& temp, commonstruct& common, Int backend);

    using FaceCoupledDriverFn =
        void (*)(dstype* fg, const dstype* xg, const dstype* ug1,
                 const dstype* ug2, const dstype* og1, const dstype* og2,
                 const dstype* wg1, const dstype* wg2, const dstype* uh,
                 const dstype* nl, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& tmp,
                 commonstruct& common, Int ngf, Int f1, Int f2, Int backend);

    using BoundaryDriverFn =
        void (*)(dstype* fb, const dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nl, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int ngf, Int f1, Int f2, Int ib,
                 Int backend);

    using UhatDriverFn =
        void (*)(dstype* fg, dstype* xg, dstype* ug1, dstype* ug2,
                 const dstype* og1, const dstype* og2, const dstype* wg1,
                 const dstype* wg2, const dstype* uh, const dstype* nl,
                 meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& tmp, commonstruct& common,
                 Int ngf, Int f1, Int f2, Int backend);

    using InitDriverFn =
        void (*)(dstype* f, const dstype* xg, appstruct& app, Int ncx,
                 Int nfield, Int npe, Int ne, Int backend);

    using FluxJacDriverFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg,
                 const dstype* xg, dstype* udg, const dstype* odg,
                 const dstype* wdg, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int nge, Int e1, Int e2,
                 Int backend);

    using ElementJacDriverFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg,
                 const dstype* xg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int nge, Int e1, Int e2,
                 Int backend);

    using SourcewWdgDriverFn =
        void (*)(dstype* f, dstype* f_wdg, const dstype* xg,
                 const dstype* udg, const dstype* odg, const dstype* wdg,
                 meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& temp, commonstruct& common,
                 Int nge, Int e1, Int e2, Int backend);

    using BoundaryJacDriverFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                 dstype* xg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, dstype* uhg, const dstype* nl,
                 meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& temp, commonstruct& common,
                 Int nga, Int ib, Int backend);

    using BoundaryStateDriverFn =
        void (*)(dstype* f, dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, dstype* uhg,
                 const dstype* nl, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int nga, Int ib, Int backend);

    using BoundaryExternalJacDriverFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                 dstype* xg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, dstype* uhg, const dstype* nl,
                 const dstype* uext, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int nga, Int ib, Int backend);

    using BoundaryExternalStateDriverFn =
        void (*)(dstype* f, dstype* xg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, dstype* uhg,
                 const dstype* nl, const dstype* uext, meshstruct& mesh,
                 masterstruct& master, appstruct& app, solstruct& sol,
                 tempstruct& temp, commonstruct& common, Int nga, Int ib,
                 Int backend);

    using FaceJacDriverFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                 const dstype* xg, dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, dstype* nl,
                 meshstruct& mesh, masterstruct& master, appstruct& app,
                 solstruct& sol, tempstruct& temp, commonstruct& common,
                 Int nga, Int backend);

    using FaceStateDriverFn =
        void (*)(dstype* f, dstype* u, const dstype* xg, dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 dstype* nl, meshstruct& mesh, masterstruct& master,
                 appstruct& app, solstruct& sol, tempstruct& temp,
                 commonstruct& common, Int nga, Int backend);

    ElementDriverFn flux = nullptr;
    ElementDriverFn source = nullptr;
    ElementDriverFn sourcew = nullptr;
    GlobalElementDriverFn output = nullptr;
    MonitorDriverFn monitor = nullptr;
    GlobalElementDriverFn avfield = nullptr;
    ElementDriverFn eos = nullptr;
    ElementDriverFn eosdu = nullptr;
    ElementDriverFn eosdw = nullptr;
    ElementDriverFn tdfunc = nullptr;
    FaceCoupledDriverFn fhatLDG = nullptr;
    BoundaryDriverFn fbouLDG = nullptr;
    UhatDriverFn uhat = nullptr;
    BoundaryDriverFn ubou = nullptr;

    InitDriverFn initodg = nullptr;
    InitDriverFn initq = nullptr;
    InitDriverFn initudg = nullptr;
    InitDriverFn initu = nullptr;
    InitDriverFn initwdg = nullptr;
    InitDriverFn cpuInitodg = nullptr;
    InitDriverFn cpuInitq = nullptr;
    InitDriverFn cpuInitudg = nullptr;
    InitDriverFn cpuInitu = nullptr;
    InitDriverFn cpuInitwdg = nullptr;

    FluxJacDriverFn fluxJac = nullptr;
    ElementJacDriverFn sourceJac = nullptr;
    ElementJacDriverFn sourcewJac = nullptr;
    SourcewWdgDriverFn sourcewWdg = nullptr;
    ElementJacDriverFn eosJac = nullptr;
    BoundaryJacDriverFn fbouHDGJac = nullptr;
    BoundaryStateDriverFn fbouHDG = nullptr;
    BoundaryJacDriverFn fintJac = nullptr;
    BoundaryStateDriverFn fintState = nullptr;
    BoundaryExternalJacDriverFn fextJac = nullptr;
    BoundaryExternalStateDriverFn fextState = nullptr;
    FaceJacDriverFn fhatHDGJac = nullptr;
    FaceStateDriverFn fhatHDG = nullptr;

    ElementDriverFn visScalars = nullptr;
    ElementDriverFn visVectors = nullptr;
    ElementDriverFn visTensors = nullptr;
    ElementDriverFn qoiVolume = nullptr;
    BoundaryDriverFn qoiBoundary = nullptr;

    [[nodiscard]] bool hasCpuInitialization() const {
        return cpuInitu != nullptr || cpuInitq != nullptr ||
               cpuInitudg != nullptr || cpuInitodg != nullptr ||
               cpuInitwdg != nullptr;
    }
};

} // namespace exasim

#endif
