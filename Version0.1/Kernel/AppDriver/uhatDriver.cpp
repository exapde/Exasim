#ifndef __UHATDRIVER
#define __UHATDRIVER

void UhatDriver(dstype *fg, dstype *xg, dstype *ug1, dstype *ug2, dstype * og1, 
     dstype * og2, dstype *wg1, dstype *wg2, dstype *uh, dstype *nl, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, 
     tempstruct &tmp, commonstruct &common, Int ngf, Int f1, Int f2, Int backend)
{        
    
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension            
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time; 
    
    if (common.extUhat==1) { 
#ifdef HAVE_ONETHREAD            
        if (backend==0) {
            opuUhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif              
#ifdef HAVE_OPENMP                        
        if (backend==1) {
            cpuUhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif            
#ifdef HAVE_CUDA                             
        if (backend==2) {
            gpuUhat(fg, xg, ug1, ug2, og1, og2, wg1, wg2, uh, nl, app.tau, app.uinf, app.physicsparam, 
                    time, numPoints, nc, ncu, nd, ncx, nco, ncw);
        }
#endif                            
    }
    else {         
        ArrayAXPBY(fg, ug1, ug2, (dstype) 0.5, (dstype) 0.5, ngf*common.ncu*(f2-f1), backend);
    }
    
#ifdef CHECK_NAN                
    dstype nrmf = PNORM(common.cublasHandle, numPoints*ncu, fg, common.backend);
    if (isnan(nrmf) || nrmf > 1.0e14) {
        cout<<"Processor: "<<common.mpiRank<<", uhat norm: "<<nrmf<<endl;
        error("here");
    }
#endif    
    
}

#endif
