/**
 * @file provider.cpp
 * @brief Low-level FrontendGenerated provider ABI export for Exasim.
 *
 * This translation unit exposes only generated low-level model kernels through
 * GetExasimDriverABI(). It intentionally avoids backend-facing driver wrappers and
 * therefore stays independent of mesh/master/sol/temp/common runtime state.
 */

#include "../ModelDispatch/driver_abi.h"

namespace frontend_generated_source {

#include "KokkosFlux.cpp"
#include "KokkosFhat.cpp"
#include "KokkosFbou.cpp"
#include "KokkosUbou.cpp"
#include "KokkosUhat.cpp"
#include "KokkosStab.cpp"
#include "KokkosSource.cpp"
#include "KokkosVisScalars.cpp"
#include "KokkosVisVectors.cpp"
#include "KokkosVisTensors.cpp"
#include "KokkosQoIvolume.cpp"
#include "KokkosQoIboundary.cpp"
#include "KokkosSourcew.cpp"
#include "KokkosOutput.cpp"
#include "KokkosMonitor.cpp"
#include "KokkosInitu.cpp"
#include "KokkosInitq.cpp"
#include "KokkosInitwdg.cpp"
#include "KokkosInitudg.cpp"
#include "KokkosInitodg.cpp"
#include "KokkosEoS.cpp"
#include "KokkosEoSdu.cpp"
#include "KokkosEoSdw.cpp"
#include "KokkosAvfield.cpp"
#include "KokkosTdfunc.cpp"

#include "cpuInitu.cpp"
#include "cpuInitq.cpp"
#include "cpuInitwdg.cpp"
#include "cpuInitudg.cpp"
#include "cpuInitodg.cpp"

#include "HdgFlux.cpp"
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgSourcewonly.cpp"
#include "HdgFbou.cpp"
#include "HdgFbouonly.cpp"
#include "HdgFint.cpp"
#include "HdgFintonly.cpp"
#include "HdgFext.cpp"
#include "HdgFextonly.cpp"
#include "HdgEoS.cpp"

} // namespace frontend_generated_source

const ExasimDriverABI& getFrontendGeneratedExasimDriverABI()
{
    static const ExasimDriverABI abi = [] {
        ExasimDriverABI value;
        value.abi_version = kExasimDriverABIVersion;
        value.struct_size = sizeof(ExasimDriverABI);

        value.KokkosFlux = &frontend_generated_source::KokkosFlux;
        value.KokkosSource = &frontend_generated_source::KokkosSource;
        value.KokkosSourcew = &frontend_generated_source::KokkosSourcew;
        value.KokkosTdfunc = &frontend_generated_source::KokkosTdfunc;
        value.KokkosAvfield = &frontend_generated_source::KokkosAvfield;
        value.KokkosEoS = &frontend_generated_source::KokkosEoS;
        value.KokkosEoSdu = &frontend_generated_source::KokkosEoSdu;
        value.KokkosEoSdw = &frontend_generated_source::KokkosEoSdw;
        value.KokkosFbou = &frontend_generated_source::KokkosFbou;
        value.KokkosUbou = &frontend_generated_source::KokkosUbou;
        value.KokkosFhat = &frontend_generated_source::KokkosFhat;
        value.KokkosUhat = &frontend_generated_source::KokkosUhat;
        value.KokkosStab = &frontend_generated_source::KokkosStab;
        value.KokkosOutput = &frontend_generated_source::KokkosOutput;
        value.KokkosMonitor = &frontend_generated_source::KokkosMonitor;
        value.KokkosVisScalars = &frontend_generated_source::KokkosVisScalars;
        value.KokkosVisVectors = &frontend_generated_source::KokkosVisVectors;
        value.KokkosVisTensors = &frontend_generated_source::KokkosVisTensors;
        value.KokkosQoIvolume = &frontend_generated_source::KokkosQoIvolume;
        value.KokkosQoIboundary = &frontend_generated_source::KokkosQoIboundary;

        value.KokkosInitu = &frontend_generated_source::KokkosInitu;
        value.KokkosInitq = &frontend_generated_source::KokkosInitq;
        value.KokkosInitudg = &frontend_generated_source::KokkosInitudg;
        value.KokkosInitwdg = &frontend_generated_source::KokkosInitwdg;
        value.KokkosInitodg = &frontend_generated_source::KokkosInitodg;
        value.cpuInitu = &frontend_generated_source::cpuInitu;
        value.cpuInitq = &frontend_generated_source::cpuInitq;
        value.cpuInitudg = &frontend_generated_source::cpuInitudg;
        value.cpuInitwdg = &frontend_generated_source::cpuInitwdg;
        value.cpuInitodg = &frontend_generated_source::cpuInitodg;

        value.HdgFlux = &frontend_generated_source::HdgFlux;
        value.HdgSource = &frontend_generated_source::HdgSource;
        value.HdgSourcew = &frontend_generated_source::HdgSourcew;
        value.HdgSourcewonly = &frontend_generated_source::HdgSourcewonly;
        value.HdgEoS = &frontend_generated_source::HdgEoS;
        value.HdgFbou = &frontend_generated_source::HdgFbou;
        value.HdgFbouonly = &frontend_generated_source::HdgFbouonly;
        value.HdgFint = &frontend_generated_source::HdgFint;
        value.HdgFintonly = &frontend_generated_source::HdgFintonly;
        value.HdgFext = &frontend_generated_source::HdgFext;
        value.HdgFextonly = &frontend_generated_source::HdgFextonly;

        return value;
    }();

    return abi;
}
