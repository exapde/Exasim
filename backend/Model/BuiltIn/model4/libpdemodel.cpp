#include <cmath>
#include <Kokkos_Core.hpp>

#ifdef USE_FLOAT
typedef float dstype;
#else
typedef double dstype; //  double is default precision 
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
