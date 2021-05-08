#ifndef __INITODGDRIVER
#define __INITODGDRIVER

void InitodgDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nco, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuInitodg(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nco, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitodg(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nco, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitodg(f, xg, app.uinf, app.physicsparam, numPoints, ncx, nco, npe, ne);             
    }
#endif   
    
#ifdef CHECK_NAN                
    dstype nrmf = PNORM(0, numPoints*nco, f, backend);
    if (isnan(nrmf) || nrmf > 1.0e14) {
        cout<<"odg0 norm: "<<nrmf<<endl;
        error("here");
    }
#endif    
    
}

#endif
