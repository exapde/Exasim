#ifndef __FLUXDRIVER
#define __FLUXDRIVER

#include <stdio.h>
void FluxDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int numPoints = nge*(e2-e1);              
    dstype time = common.time;

    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time,
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                     common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }
#endif    
}

#ifdef HAVE_ENZYME      

void FluxDriver(dstype *f, dstype *df, dstype *xg, dstype *udg, dstype *dudg, dstype *odg, dstype *dodg, dstype *wdg, 
        dstype *dwdg, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int numPoints = nge*(e2-e1);              
    dstype time = common.time;
    

    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int N = numPoints*ncu;
    Int Nf = N * nd;
    // ArraySetValue(df, numPoints*ncu, 0.0, backend);
    
    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuFluxEnzyme(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw); 
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuFluxEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuFluxEnzyme(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);               
    }
#endif           
}
#endif        

#endif
