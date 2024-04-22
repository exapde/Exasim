#ifndef __INITUDRIVER
#define __INITUDRIVER

void InituDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nc, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              
    Int modelnumber = app.flag[12];
    
    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne, app.mix); 
        #else
        opuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);                
        #endif
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitu(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nc, npe, ne);             
    }
#endif    
}

#endif
