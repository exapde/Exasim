/**
 * @file provider.cpp
 * @brief Low-level Text2codeGenerated provider ABI export for Exasim.
 *
 * This translation unit exposes only generated low-level model kernels through
 * the provider ABI table. It intentionally avoids backend-facing driver
 * wrappers and therefore stays independent of mesh/master/sol/temp/common
 * runtime state.
 */

#include "../ModelDispatch/driver_abi.hpp"

namespace exasim::text2code_generated_source {

#include <cmath>
#include <Kokkos_Core.hpp>

#ifdef USE_FLOAT
typedef float dstype;
#else
typedef double dstype;
#endif

using namespace std;

#include "KokkosFlux.cpp"
#include "KokkosFhat.cpp"
#include "KokkosFbou.cpp"
#include "KokkosUbou.cpp"
#include "KokkosUhat.cpp"
#include "KokkosStab.cpp"
#include "KokkosSource.cpp"
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
#include "HdgFbou.cpp"
#include "HdgFbouonly.cpp"
#include "HdgFint.cpp"
#include "HdgFintonly.cpp"
#include "HdgFext.cpp"
#include "HdgFextonly.cpp"
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgSourcewonly.cpp"
#include "HdgEoS.cpp"
#include "KokkosVisScalars.cpp"
#include "KokkosVisVectors.cpp"
#include "KokkosVisTensors.cpp"
#include "KokkosQoIvolume.cpp"
#include "KokkosQoIboundary.cpp"

} // namespace exasim::text2code_generated_source

namespace exasim {

namespace {

const DriverABI& getText2codeGeneratedDriverABI()
{
    static const DriverABI abi = [] {
        DriverABI value;
        value.abi_version = kDriverAbiVersion;
        value.struct_size = sizeof(DriverABI);

        value.KokkosFlux = &text2code_generated_source::KokkosFlux;
        value.KokkosSource = &text2code_generated_source::KokkosSource;
        value.KokkosSourcew = &text2code_generated_source::KokkosSourcew;
        value.KokkosTdfunc = &text2code_generated_source::KokkosTdfunc;
        value.KokkosAvfield = &text2code_generated_source::KokkosAvfield;
        value.KokkosEoS = &text2code_generated_source::KokkosEoS;
        value.KokkosEoSdu = &text2code_generated_source::KokkosEoSdu;
        value.KokkosEoSdw = &text2code_generated_source::KokkosEoSdw;
        value.KokkosFbou = &text2code_generated_source::KokkosFbou;
        value.KokkosUbou = &text2code_generated_source::KokkosUbou;
        value.KokkosFhat = &text2code_generated_source::KokkosFhat;
        value.KokkosUhat = &text2code_generated_source::KokkosUhat;
        value.KokkosStab = &text2code_generated_source::KokkosStab;
        value.KokkosOutput = &text2code_generated_source::KokkosOutput;
        value.KokkosMonitor = &text2code_generated_source::KokkosMonitor;
        value.KokkosVisScalars = &text2code_generated_source::KokkosVisScalars;
        value.KokkosVisVectors = &text2code_generated_source::KokkosVisVectors;
        value.KokkosVisTensors = &text2code_generated_source::KokkosVisTensors;
        value.KokkosQoIvolume = &text2code_generated_source::KokkosQoIvolume;
        value.KokkosQoIboundary = &text2code_generated_source::KokkosQoIboundary;

        value.KokkosInitu = &text2code_generated_source::KokkosInitu;
        value.KokkosInitq = &text2code_generated_source::KokkosInitq;
        value.KokkosInitudg = &text2code_generated_source::KokkosInitudg;
        value.KokkosInitwdg = &text2code_generated_source::KokkosInitwdg;
        value.KokkosInitodg = &text2code_generated_source::KokkosInitodg;
        value.cpuInitu = &text2code_generated_source::cpuInitu;
        value.cpuInitq = &text2code_generated_source::cpuInitq;
        value.cpuInitudg = &text2code_generated_source::cpuInitudg;
        value.cpuInitwdg = &text2code_generated_source::cpuInitwdg;
        value.cpuInitodg = &text2code_generated_source::cpuInitodg;

        value.HdgFlux = &text2code_generated_source::HdgFlux;
        value.HdgSource = &text2code_generated_source::HdgSource;
        value.HdgSourcew = &text2code_generated_source::HdgSourcew;
        value.HdgSourcewonly = &text2code_generated_source::HdgSourcewonly;
        value.HdgEoS = &text2code_generated_source::HdgEoS;
        value.HdgFbou = &text2code_generated_source::HdgFbou;
        value.HdgFbouonly = &text2code_generated_source::HdgFbouonly;
        value.HdgFint = &text2code_generated_source::HdgFint;
        value.HdgFintonly = &text2code_generated_source::HdgFintonly;
        value.HdgFext = &text2code_generated_source::HdgFext;
        value.HdgFextonly = &text2code_generated_source::HdgFextonly;

        return value;
    }();

    return abi;
}

} // namespace

} // namespace exasim
