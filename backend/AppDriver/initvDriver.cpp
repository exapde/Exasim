#ifndef __INITVDRIVER
#define __INITVDRIVER

void InitvDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;  
    Int modelnumber = app.flag[12];

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitv(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne, app.mix);  
        #else
        opuInitv(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);
        #endif
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitv(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitv(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
