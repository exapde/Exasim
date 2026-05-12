/**
 * @file driver_abi.hpp
 * @brief Stable low-level provider ABI for Exasim model modules.
 *
 * This header adapts the small exported-ABI pattern from StaticLibraryABI to
 * Exasim. Provider modules expose one versioned DriverABI object containing
 * typed kernel entry points. The reusable core backend can then build
 * higher-level driver wrappers on top of this low-level table without
 * compiling provider code into the core static library.
 */
#ifndef __EXASIM_DRIVER_ABI_H__
#define __EXASIM_DRIVER_ABI_H__

#include <cstdint>
#include <Kokkos_Core.hpp>

#ifdef USE_FLOAT
using dstype = float;
#else
using dstype = double;
#endif

inline constexpr std::uint32_t kExasimDriverABIVersion = 1;

struct ExasimDriverABI {
    using KokkosElementFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ng,
                 int nc, int ncu, int nd, int ncx, int nco, int ncw);

    using KokkosGlobalElementFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ng,
                 int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce,
                 int npe, int ne);

    using KokkosBoundaryFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ib,
                 int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

    using KokkosFaceCoupledFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg1,
                 const dstype* udg2, const dstype* odg1, const dstype* odg2,
                 const dstype* wdg1, const dstype* wdg2, const dstype* uhg,
                 const dstype* nlg, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ng,
                 int nc, int ncu, int nd, int ncx, int nco, int ncw);

    using KokkosInitFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* uinf,
                 const dstype* param, int modelnumber, int ng, int ncx,
                 int nfield, int npe, int ne);

    using CpuInitFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* uinf,
                 const dstype* param, int modelnumber, int ng, int ncx,
                 int nfield_runtime, int npe, int ne);

    using HdgElementJacFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg,
                 const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uinf, const dstype* param,
                 dstype time, int modelnumber, int ng, int nc, int ncu,
                 int nd, int ncx, int nco, int ncw);

    using HdgSourcewOnlyFn =
        void (*)(dstype* f, dstype* f_wdg, const dstype* xdg,
                 const dstype* udg, const dstype* odg, const dstype* wdg,
                 const dstype* uinf, const dstype* param, dstype time,
                 int modelnumber, int ng, int nc, int ncu, int nd, int ncx,
                 int nco, int ncw);

    using HdgBoundaryJacFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                 const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg,
                 const dstype* tau, const dstype* uinf, const dstype* param,
                 dstype time, int modelnumber, int ib, int ng, int nc,
                 int ncu, int nd, int ncx, int nco, int ncw);

    using HdgBoundaryStateFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ib,
                 int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

    using HdgBoundaryExternalJacFn =
        void (*)(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                 const dstype* xdg, const dstype* udg, const dstype* odg,
                 const dstype* wdg, const dstype* uhg, const dstype* nlg,
                 const dstype* uext, const dstype* tau, const dstype* uinf,
                 const dstype* param, dstype time, int modelnumber, int ib,
                 int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);

    using HdgBoundaryExternalStateFn =
        void (*)(dstype* f, const dstype* xdg, const dstype* udg,
                 const dstype* odg, const dstype* wdg, const dstype* uhg,
                 const dstype* nlg, const dstype* uext, const dstype* tau,
                 const dstype* uinf, const dstype* param, dstype time,
                 int modelnumber, int ib, int ng, int nc, int ncu, int nd,
                 int ncx, int nco, int ncw);

    std::uint32_t abi_version = 0;
    std::uint32_t struct_size = 0;

    KokkosElementFn KokkosFlux = nullptr;
    KokkosElementFn KokkosSource = nullptr;
    KokkosGlobalElementFn KokkosSourcew = nullptr;
    KokkosElementFn KokkosTdfunc = nullptr;
    KokkosGlobalElementFn KokkosAvfield = nullptr;
    KokkosGlobalElementFn KokkosEoS = nullptr;
    KokkosGlobalElementFn KokkosEoSdu = nullptr;
    KokkosGlobalElementFn KokkosEoSdw = nullptr;
    KokkosBoundaryFn KokkosFbou = nullptr;
    KokkosBoundaryFn KokkosUbou = nullptr;
    KokkosFaceCoupledFn KokkosFhat = nullptr;
    KokkosFaceCoupledFn KokkosUhat = nullptr;
    KokkosFaceCoupledFn KokkosStab = nullptr;
    KokkosGlobalElementFn KokkosOutput = nullptr;
    KokkosGlobalElementFn KokkosMonitor = nullptr;
    KokkosElementFn KokkosVisScalars = nullptr;
    KokkosElementFn KokkosVisVectors = nullptr;
    KokkosElementFn KokkosVisTensors = nullptr;
    KokkosElementFn KokkosQoIvolume = nullptr;
    KokkosBoundaryFn KokkosQoIboundary = nullptr;

    KokkosInitFn KokkosInitu = nullptr;
    KokkosInitFn KokkosInitq = nullptr;
    KokkosInitFn KokkosInitudg = nullptr;
    KokkosInitFn KokkosInitwdg = nullptr;
    KokkosInitFn KokkosInitodg = nullptr;
    CpuInitFn cpuInitu = nullptr;
    CpuInitFn cpuInitq = nullptr;
    CpuInitFn cpuInitudg = nullptr;
    CpuInitFn cpuInitwdg = nullptr;
    CpuInitFn cpuInitodg = nullptr;

    HdgElementJacFn HdgFlux = nullptr;
    HdgElementJacFn HdgSource = nullptr;
    HdgElementJacFn HdgSourcew = nullptr;
    HdgSourcewOnlyFn HdgSourcewonly = nullptr;
    HdgElementJacFn HdgEoS = nullptr;
    HdgBoundaryJacFn HdgFbou = nullptr;
    HdgBoundaryStateFn HdgFbouonly = nullptr;
    HdgBoundaryJacFn HdgFint = nullptr;
    HdgBoundaryStateFn HdgFintonly = nullptr;
    HdgBoundaryExternalJacFn HdgFext = nullptr;
    HdgBoundaryExternalStateFn HdgFextonly = nullptr;
};

#endif
