/**
 * @file provider.hpp
 * @brief Low-level kernel provider ABI export for Exasim.
 *
 * This translation unit intentionally exposes only model metadata and
 * low-level kernel entry points. It does not include backend-facing driver
 * wrappers such as FluxDriver/SourceDriver, so it stays independent of
 * mesh/master/sol/temp/common runtime state. Higher-level adapters in the
 * reusable core library can dispatch through this ABI later.
 */

#pragma once

#include "driver_abi.hpp"
#include "modeldefaults.hpp"

#ifndef EXASIM_KOKKOSKERNELS_HEADER
#define EXASIM_KOKKOSKERNELS_HEADER "kokkoskernels.hpp"
#endif

#include EXASIM_KOKKOSKERNELS_HEADER

namespace kokkos_kernel_source {

using ::PdeModel;

#include "kernels/KokkosFlux.hpp"
#include "kernels/KokkosFhat.hpp"
#include "kernels/KokkosFbou.hpp"
#include "kernels/KokkosUbou.hpp"
#include "kernels/KokkosUhat.hpp"
#include "kernels/KokkosStab.hpp"
#include "kernels/KokkosSource.hpp"
#include "kernels/KokkosVisScalars.hpp"
#include "kernels/KokkosVisVectors.hpp"
#include "kernels/KokkosVisTensors.hpp"
#include "kernels/KokkosQoIvolume.hpp"
#include "kernels/KokkosQoIboundary.hpp"
#include "kernels/KokkosSourcew.hpp"
#include "kernels/KokkosOutput.hpp"
#include "kernels/KokkosMonitor.hpp"
#include "kernels/KokkosInitu.hpp"
#include "kernels/KokkosInitq.hpp"
#include "kernels/KokkosInitwdg.hpp"
#include "kernels/KokkosInitudg.hpp"
#include "kernels/KokkosInitodg.hpp"
#include "kernels/KokkosEoS.hpp"
#include "kernels/KokkosEoSdu.hpp"
#include "kernels/KokkosEoSdw.hpp"
#include "kernels/KokkosAvfield.hpp"
#include "kernels/KokkosTdfunc.hpp"

#include "kernels/cpuInitu.hpp"
#include "kernels/cpuInitq.hpp"
#include "kernels/cpuInitwdg.hpp"
#include "kernels/cpuInitudg.hpp"
#include "kernels/cpuInitodg.hpp"

#include "kernels/HdgFlux.hpp"
#include "kernels/HdgSource.hpp"
#include "kernels/HdgSourcew.hpp"
#include "kernels/HdgSourcewonly.hpp"
#include "kernels/HdgFbou.hpp"
#include "kernels/HdgFbouonly.hpp"
#include "kernels/HdgFint.hpp"
#include "kernels/HdgFintonly.hpp"
#include "kernels/HdgFext.hpp"
#include "kernels/HdgFextonly.hpp"
#include "kernels/HdgEoS.hpp"

} 

inline const ExasimDriverABI& getKokkosKernelExasimDriverABI()
{
    static const ExasimDriverABI abi = [] {
        ExasimDriverABI value;
        value.abi_version = kExasimDriverABIVersion;
        value.struct_size = sizeof(ExasimDriverABI);

        value.KokkosFlux = &kokkos_kernel_source::KokkosFlux;
        value.KokkosSource = &kokkos_kernel_source::KokkosSource;
        value.KokkosSourcew = &kokkos_kernel_source::KokkosSourcew;
        value.KokkosTdfunc = &kokkos_kernel_source::KokkosTdfunc;
        value.KokkosAvfield = &kokkos_kernel_source::KokkosAvfield;
        value.KokkosEoS = &kokkos_kernel_source::KokkosEoS;
        value.KokkosEoSdu = &kokkos_kernel_source::KokkosEoSdu;
        value.KokkosEoSdw = &kokkos_kernel_source::KokkosEoSdw;
        value.KokkosFbou = &kokkos_kernel_source::KokkosFbou;
        value.KokkosUbou = &kokkos_kernel_source::KokkosUbou;
        value.KokkosFhat = &kokkos_kernel_source::KokkosFhat;
        value.KokkosUhat = &kokkos_kernel_source::KokkosUhat;
        value.KokkosStab = &kokkos_kernel_source::KokkosStab;
        value.KokkosOutput = &kokkos_kernel_source::KokkosOutput;
        value.KokkosMonitor = &kokkos_kernel_source::KokkosMonitor;
        value.KokkosVisScalars = &kokkos_kernel_source::KokkosVisScalars;
        value.KokkosVisVectors = &kokkos_kernel_source::KokkosVisVectors;
        value.KokkosVisTensors = &kokkos_kernel_source::KokkosVisTensors;
        value.KokkosQoIvolume = &kokkos_kernel_source::KokkosQoIvolume;
        value.KokkosQoIboundary = &kokkos_kernel_source::KokkosQoIboundary;

        value.KokkosInitu = &kokkos_kernel_source::KokkosInitu;
        value.KokkosInitq = &kokkos_kernel_source::KokkosInitq;
        value.KokkosInitudg = &kokkos_kernel_source::KokkosInitudg;
        value.KokkosInitwdg = &kokkos_kernel_source::KokkosInitwdg;
        value.KokkosInitodg = &kokkos_kernel_source::KokkosInitodg;
        value.cpuInitu = &kokkos_kernel_source::cpuInitu;
        value.cpuInitq = &kokkos_kernel_source::cpuInitq;
        value.cpuInitudg = &kokkos_kernel_source::cpuInitudg;
        value.cpuInitwdg = &kokkos_kernel_source::cpuInitwdg;
        value.cpuInitodg = &kokkos_kernel_source::cpuInitodg;

        value.HdgFlux = &kokkos_kernel_source::HdgFlux;
        value.HdgSource = &kokkos_kernel_source::HdgSource;
        value.HdgSourcew = &kokkos_kernel_source::HdgSourcew;
        value.HdgSourcewonly = &kokkos_kernel_source::HdgSourcewonly;
        value.HdgEoS = &kokkos_kernel_source::HdgEoS;
        value.HdgFbou = &kokkos_kernel_source::HdgFbou;
        value.HdgFbouonly = &kokkos_kernel_source::HdgFbouonly;
        value.HdgFint = &kokkos_kernel_source::HdgFint;
        value.HdgFintonly = &kokkos_kernel_source::HdgFintonly;
        value.HdgFext = &kokkos_kernel_source::HdgFext;
        value.HdgFextonly = &kokkos_kernel_source::HdgFextonly;

        return value;
    }();

    return abi;
}
