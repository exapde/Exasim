#include <type_traits>
#include "libpdemodel.hpp"
#pragma once
struct DefaultModel {

  struct Flux {
    void operator()(dstype* f, const dstype* xg, const dstype* udg,
                    const dstype* odg, const dstype* wdg,
                    const dstype* uinf, const dstype* physicsparam,
                    dstype time, Int modelnumber, Int numPoints,
                    Int nc, Int ncu, Int nd, Int ncx, Int nco, Int ncw) const
    {
      KokkosFlux(
		 f, xg, udg, odg, wdg,
		 uinf, physicsparam, time,
		 modelnumber, numPoints,
		 nc, ncu, nd, ncx, nco, ncw
		 );
    }
  };

  
  struct HdgFlux {
    void operator ()(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw) const
    {
      HdgFlux(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf, param, time,  modelnumber,  ng,  nc,  ncu,  nd,  ncx,  nco,  ncw);
    }
  };
  
  struct Source {
    void operator ()(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* physicsparam, const dstype time, const int modelnumber, const int numPoints, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw) const {
       KokkosSource(f, xdg, udg, odg, wdg, uinf, physicsparam, time, 
              modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
  };
};
