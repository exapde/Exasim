#ifndef __FLUXDRIVER
#define __FLUXDRIVER

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
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuFlux(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }
#endif    
}

#endif
