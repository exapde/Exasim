/**
 * @file provider.cpp
 * @brief Low-level UserDefined provider ABI export for Exasim.
 *
 * This translation unit intentionally exposes only model metadata and
 * low-level kernel entry points. It does not include backend-facing driver
 * wrappers such as FluxDriver/SourceDriver, so it stays independent of
 * mesh/master/sol/temp/common runtime state. Higher-level adapters in the
 * reusable core library can dispatch through this ABI later.
 */

#include "my_model.hpp"

namespace user_defined_source {

using ::PdeModel;

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

} 

const ExasimDriverABI& getUserDefinedExasimDriverABI()
{
    static const ExasimDriverABI abi = [] {
        ExasimDriverABI value;
        value.abi_version = kExasimDriverABIVersion;
        value.struct_size = sizeof(ExasimDriverABI);

        value.KokkosFlux = &user_defined_source::KokkosFlux;
        value.KokkosSource = &user_defined_source::KokkosSource;
        value.KokkosSourcew = &user_defined_source::KokkosSourcew;
        value.KokkosTdfunc = &user_defined_source::KokkosTdfunc;
        value.KokkosAvfield = &user_defined_source::KokkosAvfield;
        value.KokkosEoS = &user_defined_source::KokkosEoS;
        value.KokkosEoSdu = &user_defined_source::KokkosEoSdu;
        value.KokkosEoSdw = &user_defined_source::KokkosEoSdw;
        value.KokkosFbou = &user_defined_source::KokkosFbou;
        value.KokkosUbou = &user_defined_source::KokkosUbou;
        value.KokkosFhat = &user_defined_source::KokkosFhat;
        value.KokkosUhat = &user_defined_source::KokkosUhat;
        value.KokkosStab = &user_defined_source::KokkosStab;
        value.KokkosOutput = &user_defined_source::KokkosOutput;
        value.KokkosMonitor = &user_defined_source::KokkosMonitor;
        value.KokkosVisScalars = &user_defined_source::KokkosVisScalars;
        value.KokkosVisVectors = &user_defined_source::KokkosVisVectors;
        value.KokkosVisTensors = &user_defined_source::KokkosVisTensors;
        value.KokkosQoIvolume = &user_defined_source::KokkosQoIvolume;
        value.KokkosQoIboundary = &user_defined_source::KokkosQoIboundary;

        value.KokkosInitu = &user_defined_source::KokkosInitu;
        value.KokkosInitq = &user_defined_source::KokkosInitq;
        value.KokkosInitudg = &user_defined_source::KokkosInitudg;
        value.KokkosInitwdg = &user_defined_source::KokkosInitwdg;
        value.KokkosInitodg = &user_defined_source::KokkosInitodg;
        value.cpuInitu = &user_defined_source::cpuInitu;
        value.cpuInitq = &user_defined_source::cpuInitq;
        value.cpuInitudg = &user_defined_source::cpuInitudg;
        value.cpuInitwdg = &user_defined_source::cpuInitwdg;
        value.cpuInitodg = &user_defined_source::cpuInitodg;

        value.HdgFlux = &user_defined_source::HdgFlux;
        value.HdgSource = &user_defined_source::HdgSource;
        value.HdgSourcew = &user_defined_source::HdgSourcew;
        value.HdgSourcewonly = &user_defined_source::HdgSourcewonly;
        value.HdgEoS = &user_defined_source::HdgEoS;
        value.HdgFbou = &user_defined_source::HdgFbou;
        value.HdgFbouonly = &user_defined_source::HdgFbouonly;
        value.HdgFint = &user_defined_source::HdgFint;
        value.HdgFintonly = &user_defined_source::HdgFintonly;
        value.HdgFext = &user_defined_source::HdgFext;
        value.HdgFextonly = &user_defined_source::HdgFextonly;

        return value;
    }();

    return abi;
}
