/*
        HdgDrivers.cpp

        This file contains driver functions for HDG (Hybridizable Discontinuous Galerkin) discretization components.
        Each driver function serves as an interface to the corresponding HDG computational kernel, passing the necessary
        arguments and extracting relevant parameters from the provided data structures.

        Included HDG kernels:
                - HdgFlux: Computes flux terms for the HDG method.
                - HdgFbou: Computes boundary fluxes for the HDG method.
                - HdgSource: Computes source terms for the HDG method.
                - HdgSourcew: Computes source terms with additional variables for the HDG method.
                - HdgEoS: Computes equation of state related terms for the HDG method.

        Driver Functions:
                - FluxDriver: Calls HdgFlux to compute fluxes.
                - SourceDriver: Calls HdgSource to compute source terms.
                - SourcewDriver: Calls HdgSourcew to compute source terms with additional variables.
                - EosDriver: Calls HdgEoS to compute equation of state terms.
                - FbouDriver: Calls HdgFbou to compute boundary fluxes.

        Arguments:
                - f, fb: Output arrays for computed terms.
                - f_udg, f_wdg, f_uhg: Output arrays for derivatives with respect to udg, wdg, uhg.
                - xg: Coordinates of grid points.
                - udg, odg, wdg, uhg: Solution and auxiliary variables.
                - nl: Normal vectors at boundary faces.
                - mesh, master, app, sol, temp, common: Data structures containing mesh, master element, application, solution, temporary, and common parameters.
                - nge, ngf: Number of grid points per element/face.
                - e1, e2, f1, f2: Element/face indices.
                - ib: Boundary identifier.
                - backend: Backend identifier for computational routines.

        Notes:
                - The driver functions extract relevant parameters from the 'common' structure (such as number of components, spatial dimension, and time).
                - The number of points for computation is determined by the range of elements/faces and the number of grid points per element/face.
                - These drivers provide a modular interface for HDG computations, facilitating code organization and maintainability.
*/
#include "HdgFlux.cpp"
#include "HdgFbou.cpp"
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgEoS.cpp"

void FluxDriver<Model>(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, 
        const dstype* odg, const dstype* wdg, meshstruct &mesh, 
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

void SourceDriver<Model>(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
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

void SourcewDriver<Model>(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
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

void EosDriver<Model>(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
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

void FbouDriver<Model>(dstype* fb,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, const dstype* uhg, const dstype* nl, 
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

    HdgFbou(fb, f_udg, f_wdg, f_uhg, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                      common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
}


