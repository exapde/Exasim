#include <cmath>
#include <Kokkos_Core.hpp>



namespace Poisson3D {
  using namespace std;
#ifdef USE_FLOAT
  typedef float dstype;
#else
  typedef double dstype; //  double is default precision 
#endif

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
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgSourcewonly.cpp"
#include "HdgEoS.cpp"
#include "KokkosVisScalars.cpp"
#include "KokkosVisVectors.cpp"
#include "KokkosVisTensors.cpp"
#include "KokkosQoIvolume.cpp"
#include "KokkosQoIboundary.cpp"

  struct Poisson3D {

    // ---------------------------
    // Volume flux / source (no HDG)
    // ---------------------------

    struct FluxFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw) const
      {
	KokkosFlux(f, xdg, udg, odg, wdg,
                   uinf, param, time, modelnumber,
                   ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct SourceFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw) const
      {
	KokkosSource(f, xdg, udg, odg, wdg,
                     uinf, param, time,
                     modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct SourcewFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosSourcew(f, xdg, udg, odg, wdg,
                      uinf, param, time, modelnumber,
                      ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct AvfieldFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosAvfield(f, xdg, udg, odg, wdg,
                      uinf, param, time, modelnumber,
                      ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct EoSFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosEoS(f, xdg, udg, odg, wdg,
                  uinf, param, time, modelnumber,
                  ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct EoSduFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosEoSdu(f, xdg, udg, odg, wdg,
                    uinf, param, time, modelnumber,
                    ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct EoSdwFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosEoSdw(f, xdg, udg, odg, wdg,
                    uinf, param, time, modelnumber,
                    ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct TdfuncFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw) const
      {
	KokkosTdfunc(f, xdg, udg, odg, wdg,
                     uinf, param, time,
                     modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    // ---------------------------
    // Output / monitor
    // ---------------------------

    struct OutputFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int nce, Int npe, Int ne) const
      {
	KokkosOutput(f, xdg, udg, odg, wdg,
                     uinf, param, time, modelnumber,
                     ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
      }
    };

    struct MonitorFn {
      void operator()(dstype* f, const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int ncSol, Int ncu, Int nd, Int ncx, Int nco, Int ncw,
                      Int ncm, Int npe, Int ne) const
      {
	// note: ncSol is nc_sol argument in MonitorDriver
	KokkosMonitor(f, xdg, udg, odg, wdg,
                      uinf, param, time,
                      modelnumber, ng, ncSol, ncu, nd,
                      ncx, nco, ncw, ncm, npe, ne);
      }
    };

    // ---------------------------
    // Boundary / interface (non-HDG)
    // ---------------------------

    struct FhatFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg1, const dstype* udg2,
                      const dstype* odg1, const dstype* odg2,
                      const dstype* wdg1, const dstype* wdg2,
                      const dstype* uhg, const dstype* nlg,
                      const dstype* tau, const dstype* uinf,
                      const dstype* param, dstype time,
                      Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosFhat(f, xdg, udg1, udg2,
                   odg1, odg2, wdg1, wdg2,
                   uhg, nlg, tau, uinf, param,
                   time, modelnumber,
                   ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct FbouFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uhg,
                      const dstype* nlg, const dstype* tau,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosFbou(f, xdg, udg, odg, wdg, uhg, nlg,
                   tau, uinf, param, time,
                   modelnumber, ib, ng,
                   nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct UbouFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uhg,
                      const dstype* nlg, const dstype* tau,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosUbou(f, xdg, udg, odg, wdg, uhg, nlg,
                   tau, uinf, param, time,
                   modelnumber, ib, ng,
                   nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct UhatFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg1, const dstype* udg2,
                      const dstype* odg1, const dstype* odg2,
                      const dstype* wdg1, const dstype* wdg2,
                      const dstype* uhg, const dstype* nlg,
                      const dstype* tau, const dstype* uinf,
                      const dstype* param, dstype time,
                      Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosUhat(f, xdg, udg1, udg2,
                   odg1, odg2, wdg1, wdg2,
                   uhg, nlg, tau, uinf, param,
                   time, modelnumber,
                   ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct StabFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg1, const dstype* udg2,
                      const dstype* odg1, const dstype* odg2,
                      const dstype* wdg1, const dstype* wdg2,
                      const dstype* uhg, const dstype* nlg,
                      const dstype* tau, const dstype* uinf,
                      const dstype* param, dstype time,
                      Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosStab(f, xdg, udg1, udg2,
                   odg1, odg2, wdg1, wdg2,
                   uhg, nlg, tau, uinf, param,
                   time, modelnumber,
                   ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    // ---------------------------
    // Init (Kokkos)
    // ---------------------------

    struct InitodgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	KokkosInitodg(f, xdg, uinf, param,
                      modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct InitqFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	KokkosInitq(f, xdg, uinf, param,
                    modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct InituFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	KokkosInitu(f, xdg, uinf, param,
                    modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct InitudgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	KokkosInitudg(f, xdg, uinf, param,
                      modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct InitwdgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	KokkosInitwdg(f, xdg, uinf, param,
                      modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    // ---------------------------
    // Init (CPU)
    // ---------------------------

    struct CpuInitodgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	cpuInitodg(f, xdg, uinf, param,
                   modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct CpuInitqFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	cpuInitq(f, xdg, uinf, param,
		 modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct CpuInituFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	cpuInitu(f, xdg, uinf, param,
		 modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct CpuInitudgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	cpuInitudg(f, xdg, uinf, param,
                   modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    struct CpuInitwdgFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* uinf, const dstype* param,
                      Int modelnumber, Int ng,
                      Int ncx, Int nce, Int npe, Int ne) const
      {
	cpuInitwdg(f, xdg, uinf, param,
                   modelnumber, ng, ncx, nce, npe, ne);
      }
    };

    // ---------------------------
    // HDG volume/source/etc.
    // ---------------------------

    struct HdgFluxFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgFlux(f, f_udg, f_wdg,
		xdg, udg, odg, wdg,
		uinf, param, time,
		modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgSourceFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgSource(f, f_udg, f_wdg,
                  xdg, udg, odg, wdg,
                  uinf, param, time,
                  modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgSourcewFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgSourcew(f, f_udg, f_wdg,
                   xdg, udg, odg, wdg,
                   uinf, param, time,
                   modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgSourcewonlyFn {
      void operator()(dstype* f, dstype* f_wdg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgSourcewonly(f, f_wdg,
                       xdg, udg, odg, wdg,
                       uinf, param, time,
                       modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgEoSFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgEoS(f, f_udg, f_wdg,
               xdg, udg, odg, wdg,
               uinf, param, time,
               modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgFbouFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uhg, const dstype* nlg,
                      const dstype* tau, const dstype* uinf,
                      const dstype* param, dstype time,
                      Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgFbou(f, f_udg, f_wdg, f_uhg,
		xdg, udg, odg, wdg,
		uhg, nlg, tau, uinf, param,
		time, modelnumber, ib,
		ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgFbouonlyFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uhg,
                      const dstype* nlg, const dstype* tau,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgFbouonly(f, xdg, udg, odg, wdg,
                    uhg, nlg, tau, uinf, param,
                    time, modelnumber, ib,
                    ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgFintFn {
      void operator()(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg,
                      const dstype* xdg, const dstype* udg,
                      const dstype* odg, const dstype* wdg,
                      const dstype* uhg, const dstype* nlg,
                      const dstype* tau, const dstype* uinf,
                      const dstype* param, dstype time,
                      Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgFint(f, f_udg, f_wdg, f_uhg,
		xdg, udg, odg, wdg,
		uhg, nlg, tau, uinf, param,
		time, modelnumber, ib,
		ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct HdgFintonlyFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uhg,
                      const dstype* nlg, const dstype* tau,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	HdgFintonly(f, xdg, udg, odg, wdg,
                    uhg, nlg, tau, uinf, param,
                    time, modelnumber, ib,
                    ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };

    // ---------------------------
    // Visualization & QoI
    // ---------------------------

    struct VisScalarsFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosVisScalars(f, xdg, udg, odg, wdg,
			 uinf, param, time,
			 modelnumber, ng,
			 nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct VisVectorsFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosVisVectors(f, xdg, udg, odg, wdg,
			 uinf, param, time,
			 modelnumber, ng,
			 nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct VisTensorsFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosVisTensors(f, xdg, udg, odg, wdg,
			 uinf, param, time,
			 modelnumber, ng,
			 nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct QoIvolumeFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosQoIvolume(f, xdg, udg, odg, wdg,
			uinf, param, time,
			modelnumber, ng,
			nc, ncu, nd, ncx, nco, ncw);
      }
    };

    struct QoIboundaryFn {
      void operator()(dstype* f, const dstype* xdg,
                      const dstype* udg, const dstype* odg,
                      const dstype* wdg, const dstype* uhg,
                      const dstype* nlg, const dstype* tau,
                      const dstype* uinf, const dstype* param,
                      dstype time, Int modelnumber, Int ib, Int ng,
                      Int nc, Int ncu, Int nd,
                      Int ncx, Int nco, Int ncw) const
      {
	KokkosQoIboundary(f, xdg, udg, odg, wdg,
                          uhg, nlg, tau, uinf, param,
                          time, modelnumber, ib,
                          ng, nc, ncu, nd, ncx, nco, ncw);
      }
    };
  };
}
