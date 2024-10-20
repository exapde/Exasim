#include "HdgFlux.cpp"
#include "HdgFbou.cpp"
#include "HdgSource.cpp"
#include "HdgSourcew.cpp"
#include "HdgEoS.cpp"

void FluxDriver(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xg, const dstype* udg, const dstype* odg, const dstype* wdg, meshstruct &mesh, 
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

void FbouDriver(dstype* fb,  dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xg, const dstype* udg, const dstype*  odg, const dstype*  wdg, const dstype* uhg, const dstype* nl, 
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


