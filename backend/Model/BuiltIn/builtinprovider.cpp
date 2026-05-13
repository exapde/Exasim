/**
 * @file provider.cpp
 * @brief Low-level BuiltIn provider ABI export for Exasim.
 *
 * This translation unit binds the reusable ExasimDriverABI table directly to the
 * built-in low-level kernel dispatch functions implemented in
 * libbuiltinmodel.cpp. It stays independent of backend-facing driver wrappers
 * such as BuiltinModelDrivers.cpp and therefore does not depend on
 * mesh/master/sol/temp/common runtime state.
 */

#include "../ModelDispatch/driver_abi.h"
#include "builtinmodelsource.cpp"

const ExasimDriverABI& getBuiltInExasimDriverABI()
{
    static const ExasimDriverABI abi = [] {
        ExasimDriverABI value;
        value.abi_version = kExasimDriverABIVersion;
        value.struct_size = sizeof(ExasimDriverABI);

        value.KokkosFlux = &::builtinKokkosFlux;
        value.KokkosSource = &::builtinKokkosSource;
        value.KokkosSourcew = &::builtinKokkosSourcew;
        value.KokkosTdfunc = &::builtinKokkosTdfunc;
        value.KokkosAvfield = &::builtinKokkosAvfield;
        value.KokkosEoS = &::builtinKokkosEoS;
        value.KokkosEoSdu = &::builtinKokkosEoSdu;
        value.KokkosEoSdw = &::builtinKokkosEoSdw;
        value.KokkosFbou = &::builtinKokkosFbou;
        value.KokkosUbou = &::builtinKokkosUbou;
        value.KokkosFhat = &::builtinKokkosFhat;
        value.KokkosUhat = &::builtinKokkosUhat;
        value.KokkosStab = &::builtinKokkosStab;
        value.KokkosOutput = &::builtinKokkosOutput;
        value.KokkosMonitor = &::builtinKokkosMonitor;
        value.KokkosVisScalars = &::builtinKokkosVisScalars;
        value.KokkosVisVectors = &::builtinKokkosVisVectors;
        value.KokkosVisTensors = &::builtinKokkosVisTensors;
        value.KokkosQoIvolume = &::builtinKokkosQoIvolume;
        value.KokkosQoIboundary = &::builtinKokkosQoIboundary;

        value.KokkosInitu = &::builtinKokkosInitu;
        value.KokkosInitq = &::builtinKokkosInitq;
        value.KokkosInitudg = &::builtinKokkosInitudg;
        value.KokkosInitwdg = &::builtinKokkosInitwdg;
        value.KokkosInitodg = &::builtinKokkosInitodg;
        value.cpuInitu = &::builtincpuInitu;
        value.cpuInitq = &::builtincpuInitq;
        value.cpuInitudg = &::builtincpuInitudg;
        value.cpuInitwdg = &::builtincpuInitwdg;
        value.cpuInitodg = &::builtincpuInitodg;

        value.HdgFlux = &::builtinHdgFlux;
        value.HdgSource = &::builtinHdgSource;
        value.HdgSourcew = &::builtinHdgSourcew;
        value.HdgSourcewonly = &::builtinHdgSourcewonly;
        value.HdgEoS = &::builtinHdgEoS;
        value.HdgFbou = &::builtinHdgFbou;
        value.HdgFbouonly = &::builtinHdgFbouonly;
        value.HdgFint = &::builtinHdgFint;
        value.HdgFintonly = &::builtinHdgFintonly;
        value.HdgFext = &::builtinHdgFext;
        value.HdgFextonly = &::builtinHdgFextonly;

        return value;
    }();

    return abi;
}
