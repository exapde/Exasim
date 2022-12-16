#ifndef __INITUDGDRIVER
#define __INITUDGDRIVER

void InitudgDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              
    Int modelnumber = app.flag[12];

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne, app.mix);
        #else
        opuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);
        #endif                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitudg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
