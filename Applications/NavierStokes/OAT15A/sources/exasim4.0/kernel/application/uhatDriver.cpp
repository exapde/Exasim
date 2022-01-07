#ifndef __UHATDRIVER
#define __UHATDRIVER

void UhatDriver(dstype *fg, dstype *xg, dstype *uh, dstype *ug1, dstype *ug2, dstype * og1, 
     dstype * og2, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)        
{        
    
    ArrayAXPBY(fg, ug1, ug2, (dstype) 0.5, (dstype) 0.5, ngf*common.ncu*(f2-f1), backend);
}

#endif
