#ifndef __TDFUNCDRIVER
#define __TDFUNCDRIVER

void TdfuncDriver(dstype *s, dstype *xg, dstype *udg, dstype *odg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int ncp = common.ncp;// number of compoments of (p)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = nge*(e2-e1);              
    dstype time = common.time;

    /* 2. Compute physical tdfunc */
#ifdef HAVE_ONETHREAD                 
    if (backend==0) {
        opuTdfunc(s, xg, udg, odg, &app.physicsparam[0], time, 
                            numPoints, nc, ncu, nd, ncx, nco, common.appname);        
    }
#endif              
#ifdef HAVE_OPENMP                    
    if (backend==1) {
        cpuTdfunc(s, xg, udg, odg, &app.physicsparam[0], time, 
                            numPoints, nc, ncu, nd, ncx, nco, common.appname);
    }
#endif            
#ifdef HAVE_CUDA                         
    if (backend==2) {
        gpuTdfunc(s, xg, udg, odg, &app.physicsparam[0], time, 
                            numPoints, nc, ncu, nd, ncx, nco, common.appname);
    }
#endif                
}

#endif
