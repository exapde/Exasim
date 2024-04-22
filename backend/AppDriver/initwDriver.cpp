#ifndef __INITWDRIVER
#define __INITWDRIVER

void InitwDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              
    Int modelnumber = app.flag[12];
    
    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitw(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne, app.mix);
        #else
        opuInitw(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);
        #endif                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitw(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitw(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
