#ifndef __INITODGDRIVER
#define __INITODGDRIVER

void InitodgDriver(dstype *f, dstype *xg, appstruct &app, Int ncx, Int nco, Int npe, Int ne, Int backend)
{     
    Int numPoints = npe*ne;              
    Int modelnumber = app.flag[12];
    
    /* 2. Compute output field */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        #ifdef HAVE_MPP
        opuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne, app.mix);        
        #else        
        opuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne);   
        #endif
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne);             
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuInitodg(f, xg, app.uinf, app.physicsparam, modelnumber, numPoints, ncx, nco, npe, ne);             
    }
#endif    
}

#endif
