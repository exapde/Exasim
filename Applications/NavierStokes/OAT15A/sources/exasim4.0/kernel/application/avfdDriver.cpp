#ifndef __AVFDDRIVER
#define __AVFDDRIVER

void AvfdDriver(dstype *avField, dstype *xg, dstype *udg, dstype *odg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int backend)
{ 
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int ncp = common.ncp;// number of compoments of (p)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe*ne;              
    dstype time = common.time;

    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuAVfield(avField, xg, udg, odg, &app.physicsparam[0], time, 
                    numPoints, nc, ncu, nd, ncx, nco, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuAVfield(avField, xg, udg, odg, &app.physicsparam[0], time, 
                    numPoints, nc, ncu, nd, ncx, nco, npe, ne);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuAVfield(avField, xg, udg, odg, &app.physicsparam[0], time, 
                    numPoints, nc, ncu, nd, ncx, nco, npe, ne);                        
    }
#endif    
}

#endif
