#ifndef __INITQDRIVER
#define __INITQDRIVER

void InitqDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuInitq(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitq(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitq(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nc, npe, ne);             
    }
#endif    
    
#ifdef CHECK_NAN                
    dstype nrmf = PNORM(common.cublasHandle, numPoints*nc, f, common.backend);
    if (isnan(nrmf) || nrmf > 1.0e14) {
        cout<<"Processor: "<<common.mpiRank<<", q0 norm: "<<nrmf<<endl;
        error("here");
    }
#endif    
    
}

#endif
