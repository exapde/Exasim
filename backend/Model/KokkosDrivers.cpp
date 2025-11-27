
/*
  This file contains driver functions for various numerical operations in the Exasim backend,
  including flux, source, output, monitor, boundary, initialization, and equation-of-state computations.
  The drivers are designed to interface with both DG-based and HDG-based implementations.

  Included Files:
  - DG-based implementations: KokkosFlux.cpp, KokkosFhat.cpp, KokkosFbou.cpp, KokkosUbou.cpp, KokkosUhat.cpp,
    KokkosStab.cpp, KokkosSource.cpp, KokkosSourcew.cpp, KokkosOutput.cpp, KokkosMonitor.cpp, KokkosInitu.cpp,
    KokkosInitq.cpp, KokkosInitwdg.cpp, KokkosInitudg.cpp, KokkosInitodg.cpp, KokkosEoS.cpp, KokkosEoSdu.cpp,
    KokkosEoSdw.cpp, KokkosAvfield.cpp, KokkosTdfunc.cpp
  - CPU-based initialization: cpuInitu.cpp, cpuInitq.cpp, cpuInitwdg.cpp, cpuInitudg.cpp, cpuInitodg.cpp
  - HDG-based implementations: HdgFlux.cpp, HdgFbou.cpp, HdgFbouonly.cpp, HdgFint.cpp, HdgFintonly.cpp,
    HdgSource.cpp, HdgSourcew.cpp, HdgSourcewonly.cpp, HdgEoS.cpp

  Driver Functions:
  - FluxDriver: Computes fluxes using Kokkos or HDG implementations.
  - SourceDriver: Computes source terms.
  - SourcewDriver: Computes source terms related to w variables.
  - OutputDriver: Computes output quantities.
  - MonitorDriver: Computes monitoring quantities.
  - AvfieldDriver: Computes artificial viscosity field.
  - EosDriver, EosduDriver, EosdwDriver: Computes equation-of-state and its derivatives.
  - TdfuncDriver: Computes time-dependent functions.
  - FhatDriver: Computes numerical fluxes at element interfaces.
  - FbouDriver: Computes boundary fluxes.
  - UhatDriver: Computes numerical traces at interfaces.
  - UbouDriver: Computes boundary traces.
  - Initialization Drivers (InitodgDriver, InitqDriver, InitudgDriver, InituDriver, InitwdgDriver): 
    Initialize solution variables.
  - CPU Initialization Drivers (cpuInitodgDriver, cpuInitqDriver, cpuInitudgDriver, cpuInituDriver, cpuInitwdgDriver):
    CPU-based initialization routines.

  HDG-specific Drivers:
  - HDG versions of FluxDriver, SourceDriver, SourcewDriver, EosDriver, FbouDriver, FintDriver, FhatDriver
    provide additional outputs such as derivatives with respect to udg, wdg, and uhg.

  Parameters:
  - All drivers accept pointers to solution arrays, mesh, master, app, sol, temp, and common structures,
    as well as integer parameters specifying grid sizes, element indices, and backend selection.

  Notes:
  - The drivers abstract away the details of the underlying numerical implementation, allowing for
    flexible switching between Kokkos and HDG approaches.
  - Some drivers are overloaded to provide additional outputs (e.g., derivatives).
  - The file is intended for internal use within the Exasim backend and is not a standalone module.
*/
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
  
void FluxDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;

  KokkosFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
}

void SourceDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  KokkosSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void QoIvolumeDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  KokkosQoIvolume(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void VisScalarsDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  KokkosVisScalars(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void VisVectorsDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  KokkosVisVectors(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}


void VisTensorsDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  KokkosVisTensors(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void SourcewDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int ne = e2-e1;
  Int numPoints = npe*ne;              
  dstype time = common.time;            
  
  KokkosSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
              numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne);                
}

void OutputDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)     
  Int nce = common.nce;// number of compoments of (e)
  Int nd = common.nd;     // spatial dimension    
  Int npe = common.npe;
  Int ne = common.ne;
  Int numPoints = npe*ne;              

  dstype time = common.time;

  KokkosOutput(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);                
}

void MonitorDriver(dstype* f, Int nc_sol, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
  masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
  commonstruct &common, Int backend)
{ 
  //Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)     
  Int ncm = common.ncm;// number of compoments of (e)
  Int nd = common.nd;     // spatial dimension    
  Int npe = common.npe;
  Int ne = common.ne;
  Int numPoints = npe*ne;              

  dstype time = common.time;

  KokkosMonitor(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
          common.modelnumber, numPoints, nc_sol, ncu, nd, ncx, nco, ncw, ncm, npe, ne);                
}

void AvfieldDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  //Int nce = common.nce;// number of compoments of (e)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int npe = common.npe;
  Int ne = common.ne;
  Int numPoints = npe*ne;              

  dstype time = common.time;

  KokkosAvfield(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber, 
              numPoints, nc, ncu, nd, ncx, nco, ncw, common.ncAV, npe, ne);                
}

void EosDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int ne = e2-e1;
  Int numPoints = npe*ne;              
  dstype time = common.time;            
  
  KokkosEoS(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
              numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne);                                          
}

void EosduDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int ne = e2-e1;
  Int numPoints = npe*ne;              
  dstype time = common.time;            
  
  KokkosEoSdu(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
              numPoints, nc, ncu, nd, ncx, nco, ncw, ncw*ncu, npe, ne);                  
}

void EosdwDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int npe, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int ne = e2-e1;
  Int numPoints = npe*ne;              
  dstype time = common.time;            
  
  KokkosEoSdw(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
              numPoints, nc, ncu, nd, ncx, nco, ncw, ncw*ncw, npe, ne);                        
}

void TdfuncDriver(dstype* f, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;

  KokkosTdfunc(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
}

void FhatDriver(dstype* fg, const dstype* xg, const dstype* ug1, const dstype* ug2, const dstype*  og1, const dstype*  og2, 
    const dstype* wg1, const dstype* wg2, const dstype* uh, const dstype* nl, meshstruct &mesh, masterstruct &master, appstruct &app, 
    solstruct &sol, tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{
  Int nc = common.nc; // number of compoments of (u, q, p)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension            
  Int numPoints = ngf*(f2-f1);
  Int M = numPoints * ncu;
  Int N = numPoints*ncu*nd;
  Int ntau = common.ntau;
  dstype time = common.time; 
  
  if (common.extFhat==1) { 
    KokkosFhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
            time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
  }
  else {     
    // left flux
    FluxDriver(fg, xg, ug1, og1, wg1, mesh, master, app, sol, tmp, common,
                    ngf, f1, f2, backend);        
    // right flux
    //auto fg2 = Kokkos::subview(fg, Kokkos::make_pair(N, 2*N));
    dstype *fg2 = &fg[N];
    FluxDriver(fg2, xg, ug2, og2, wg2, mesh, master, app, sol, tmp, common,
                    ngf, f1, f2, backend);        
    
    // Part 1: fh = fg dot nl
    AverageFlux(fg, N);    
    AverageFluxDotNormal(fg, nl, N, M, numPoints, nd);            

    // Part 2: Contribution due to tau*(U-UH)
    if (common.extStab>=1) { 
      KokkosStab(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
              time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
    else if (ntau==0) { 
      AddStabilization1(fg, ug1, ug2, app.tau, M);
    }
    else if (ntau==1) { // constant scalar  
      AddStabilization1(fg, ug1, ug2, app.tau, M);
    }
    else if (ntau==ncu) { // constant diagonal tensor
      AddStabilization2(fg, ug1, ug2, app.tau, M, numPoints);
    }
    else if (ntau==ncu*ncu) { // constant full tensor      
      AddStabilization3(fg, ug1, ug2, app.tau, M, numPoints, ncu);
    }
    else {
      printf("Stabilization option is not implemented");
      exit(-1);
    }            
  }
}

void FbouDriver(dstype* fb, const dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, const dstype* uhg, const dstype* nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    KokkosFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void QoIboundaryDriver(dstype* fb, const dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, const dstype* uhg, const dstype* nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    KokkosQoIboundary(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void UhatDriver(dstype* fg, dstype* xg, dstype* ug1, dstype* ug2, const dstype*  og1, 
     const dstype*  og2, const dstype* wg1, const dstype* wg2, const dstype* uh, const dstype* nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{        
    
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension            
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time; 
    
    if (common.extUhat==1) { 
      KokkosUhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
              time, common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
    else {         
      ArrayAXPBY(fg, ug1, ug2, (dstype) 0.5, (dstype) 0.5, ngf*common.ncu*(f2-f1));
    }
}

void UbouDriver(dstype* ub, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    KokkosUbou(ub, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void InitodgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nco, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];
  
  KokkosInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne);                
}

void InitqDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;     
  Int modelnumber = app.flag[12];

  KokkosInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
    
}

void InitudgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  KokkosInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
}

void InituDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  KokkosInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
}

void InitwdgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int ncw, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  KokkosInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne);                
}

void cpuInitodgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nco, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];
  
  cpuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne);                
}

void cpuInitqDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;     
  Int modelnumber = app.flag[12];

  cpuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
    
}

void cpuInitudgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  cpuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
}

void cpuInituDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  cpuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
}

void cpuInitwdgDriver(dstype* f, const dstype* xg, appstruct &app, Int ncx, Int ncw, Int npe, Int ne, Int backend)
{     
  Int numPoints = npe*ne;              
  Int modelnumber = app.flag[12];

  cpuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne);                
}


#include "HdgFlux.cpp"
#include "HdgFbou.cpp"
#include "HdgFbouonly.cpp"
#include "HdgFint.cpp"
#include "HdgFintonly.cpp"
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgSourcewonly.cpp"
#include "HdgEoS.cpp"

void FluxDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg,  dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;

  HdgFlux(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);         
}

void SourceDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);              
  dstype time = common.time;
            
  HdgSource(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void SourcewDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);    
  dstype time = common.time;            
  
  HdgSourcew(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void SourcewDriver(dstype* f, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);    
  dstype time = common.time;            
  
  HdgSourcewonly(f, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                       
}

void EosDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
  Int nc = common.nc; // number of compoments of (u, q)
  Int ncu = common.ncu;// number of compoments of (u)
  Int ncw = common.ncw;// number of compoments of (w)
  Int nco = common.nco;// number of compoments of (o)
  Int ncx = common.ncx;// number of compoments of (xdg)        
  Int nd = common.nd;     // spatial dimension    
  Int numPoints = nge*(e2-e1);    
  dstype time = common.time;            
  
  HdgEoS(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
              common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                           
}

void FbouDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg, dstype* xg, const dstype* udg, 
        const dstype*  odg, const dstype*  wdg, dstype* uhg, const dstype* nl, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, commonstruct &common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    dstype time = common.time;    


    HdgFbou(f, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);    
}


void FbouDriver(dstype* f, dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, dstype* uhg, 
        const dstype* nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    dstype time = common.time;    

    HdgFbouonly(f, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void FintDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg, dstype* xg, const dstype* udg, 
        const dstype*  odg, const dstype*  wdg, dstype* uhg, const dstype* nl, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, commonstruct &common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    dstype time = common.time;    


    HdgFint(f, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void FintDriver(dstype* f, dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, dstype* uhg, 
        const dstype* nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nga, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    dstype time = common.time;    

    HdgFintonly(f, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}

void FhatDriver(dstype* f,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xg, dstype* udg, 
        const dstype*  odg, const dstype*  wdg, const dstype* uhg,  dstype* nl, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, commonstruct &common, Int nga, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    Int M = numPoints * ncu;
    Int N = numPoints*ncu*nd;

    dstype time = common.time;    
    
    // copy u-compoment of udg to f_uhg
    ArrayCopy(f_uhg, udg, numPoints*ncu);

    // copy uhg to udg, so udg = (uh, q)  
    ArrayCopy(udg, uhg, numPoints*ncu);

    // f = ng * ncu * nd 
    // f_udg = ng * ncu * nd * nc
    // f_wdg = ng * ncu * nd * ncw
    HdgFlux(f, f_udg, f_wdg, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                      

    // copy f_uhg to u-compoment of udg 
    ArrayCopy(udg, f_uhg, numPoints*ncu);

    // f = ng * ncu  
    FluxDotNormal(f, f, nl, M, numPoints, nd);         

    // f_udg = ng * ncu * nc
    for (int n=0; n< nc; n++) {
      FluxDotNormal(&f_udg[M*n], &f_udg[N*n], nl, M, numPoints, nd);        
    }

    // f_wdg = ng * ncu * ncw
    if ((ncw>0) & (common.wave==0)) {
      for (int n=0; n< ncw; n++) {
        FluxDotNormal(&f_wdg[M*n], &f_wdg[N*n], nl, M, numPoints, nd);    
      }      
    }  

    // f_uhg = ng * ncu * ncu
    ArrayCopy(f_uhg, f_udg, numPoints*ncu*ncu); // copy f_u to f_uhg

    // reset f_u to zero
    ArraySetValue(f_udg, zero, numPoints*ncu*ncu);        

    // add tau*(u - uh) to f
    AddStabilization1(f, f_udg, f_uhg, udg, uhg, app.tau, M, numPoints);  
}

void FhatDriver(dstype* f, dstype* u, const dstype* xg, dstype* udg, const dstype*  odg, const dstype*  wdg, const dstype* uhg,  
        dstype* nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nga, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nga;
    Int M = numPoints * ncu;
    //Int N = numPoints*ncu*nd;

    dstype time = common.time;    

    // copy u-compoment of udg to u
    ArrayCopy(u, udg, numPoints*ncu);

    // copy uhg to udg, so udg = (uh, q)  
    ArrayCopy(udg, uhg, numPoints*ncu);

    // f = ng * ncu * nd 
    KokkosFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                      

    // copy u to u-compoment of udg 
    ArrayCopy(udg, u, numPoints*ncu);

    // f = ng * ncu  
    FluxDotNormal(f, f, nl, M, numPoints, nd);         

    // add tau*(u - uh) to f
    AddStabilization1(f, udg, uhg, app.tau, M);      
}

