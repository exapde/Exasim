#ifndef __INITUDRIVER
#define __INITUDRIVER

void InituDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuInitu(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitu(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitu(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }
#endif    
    
#ifdef CHECK_NAN                
    dstype nrmf = PNORM(common.cublasHandle, numPoints*nc, f, common.backend);
    if (isnan(nrmf) || nrmf > 1.0e14) {
        cout<<"Processor: "<<common.mpiRank<<", u0 norm: "<<nrmf<<endl;
        error("here");
    }
#endif    
    
}

#endif
