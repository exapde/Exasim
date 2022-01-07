#ifndef __OUTPUTDRIVER
#define __OUTPUTDRIVER

void OutputDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int backend)
{ 
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)     
    Int nce = common.nce;// number of compoments of (e)
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe*ne;              

    dstype time = common.time;

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuOutput(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuOutput(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuOutput(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);                        
    }
#endif    
}

#endif
