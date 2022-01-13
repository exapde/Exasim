#ifndef __SOURCEWDRIVER
#define __SOURCEWDRIVER

void SourcewDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
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
    
    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
                    numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne, app.mix);
        #else
        opuSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
                    numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne);
        #endif                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
                    numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuSourcew(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber,
                    numPoints, nc, ncu, nd, ncx, nco, ncw, ncw, npe, ne);                        
    }
#endif    
}

#endif
