#ifndef __INITQDRIVER
#define __INITQDRIVER

void InitqDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;     
    Int modelnumber = app.flag[12];

    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne, app.mix);
        #else
        opuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);
        #endif                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitq(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
