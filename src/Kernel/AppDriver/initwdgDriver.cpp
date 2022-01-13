#ifndef __INITWDGDRIVER
#define __INITWDGDRIVER

void InitwdgDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int ncw, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              
    Int modelnumber = app.flag[12];
    
    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne, app.mix);      
        #else
        opuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne);
        #endif          
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitwdg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, ncw, npe, ne);             
    }
#endif    
}

#endif
