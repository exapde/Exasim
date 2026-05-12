#include <cmath>
#include <Kokkos_Core.hpp>
#include "../ModelDispatch/driver_abi.h"

#ifdef USE_FLOAT
typedef float dstype;
#else
typedef double dstype; //  double is default precision 
#endif

using namespace std;

namespace text2code_shared_lib {

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

} // namespace text2code_shared_lib

extern "C" const ExasimDriverABI* GetText2CodeExasimDriverABI()
{
    static const ExasimDriverABI abi = [] {
        ExasimDriverABI value;
        value.abi_version = kExasimDriverABIVersion;
        value.struct_size = sizeof(ExasimDriverABI);

        value.KokkosFlux = &text2code_shared_lib::KokkosFlux;
        value.KokkosSource = &text2code_shared_lib::KokkosSource;
        value.KokkosSourcew = &text2code_shared_lib::KokkosSourcew;
        value.KokkosTdfunc = &text2code_shared_lib::KokkosTdfunc;
        value.KokkosAvfield = &text2code_shared_lib::KokkosAvfield;
        value.KokkosEoS = &text2code_shared_lib::KokkosEoS;
        value.KokkosEoSdu = &text2code_shared_lib::KokkosEoSdu;
        value.KokkosEoSdw = &text2code_shared_lib::KokkosEoSdw;
        value.KokkosFbou = &text2code_shared_lib::KokkosFbou;
        value.KokkosUbou = &text2code_shared_lib::KokkosUbou;
        value.KokkosFhat = &text2code_shared_lib::KokkosFhat;
        value.KokkosUhat = &text2code_shared_lib::KokkosUhat;
        value.KokkosStab = &text2code_shared_lib::KokkosStab;
        value.KokkosOutput = &text2code_shared_lib::KokkosOutput;
        value.KokkosMonitor = &text2code_shared_lib::KokkosMonitor;
        value.KokkosVisScalars = &text2code_shared_lib::KokkosVisScalars;
        value.KokkosVisVectors = &text2code_shared_lib::KokkosVisVectors;
        value.KokkosVisTensors = &text2code_shared_lib::KokkosVisTensors;
        value.KokkosQoIvolume = &text2code_shared_lib::KokkosQoIvolume;
        value.KokkosQoIboundary = &text2code_shared_lib::KokkosQoIboundary;

        value.KokkosInitu = &text2code_shared_lib::KokkosInitu;
        value.KokkosInitq = &text2code_shared_lib::KokkosInitq;
        value.KokkosInitudg = &text2code_shared_lib::KokkosInitudg;
        value.KokkosInitwdg = &text2code_shared_lib::KokkosInitwdg;
        value.KokkosInitodg = &text2code_shared_lib::KokkosInitodg;
        value.cpuInitu = &text2code_shared_lib::cpuInitu;
        value.cpuInitq = &text2code_shared_lib::cpuInitq;
        value.cpuInitudg = &text2code_shared_lib::cpuInitudg;
        value.cpuInitwdg = &text2code_shared_lib::cpuInitwdg;
        value.cpuInitodg = &text2code_shared_lib::cpuInitodg;

        value.HdgFlux = &text2code_shared_lib::HdgFlux;
        value.HdgSource = &text2code_shared_lib::HdgSource;
        value.HdgSourcew = &text2code_shared_lib::HdgSourcew;
        value.HdgSourcewonly = &text2code_shared_lib::HdgSourcewonly;
        value.HdgEoS = &text2code_shared_lib::HdgEoS;
        value.HdgFbou = &text2code_shared_lib::HdgFbou;
        value.HdgFbouonly = &text2code_shared_lib::HdgFbouonly;
        value.HdgFint = &text2code_shared_lib::HdgFint;
        value.HdgFintonly = &text2code_shared_lib::HdgFintonly;
        value.HdgFext = &text2code_shared_lib::HdgFext;
        value.HdgFextonly = &text2code_shared_lib::HdgFextonly;

        return value;
    }();

    return &abi;
}
