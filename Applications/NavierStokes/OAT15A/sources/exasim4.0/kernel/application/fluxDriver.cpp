#ifndef __FLUXDRIVER
#define __FLUXDRIVER

void FluxDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, meshstruct &mesh, 
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

    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
//         if (common.appname == 1)            
//             opuFluxNS(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);                
//         else if (common.appname == 2)            
//             opuFluxPoisson(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);
//         else
            opuFlux(f, xg, udg, odg, &app.physicsparam[0], time, 
                        numPoints, nc, ncu, nd, ncx, nco, common.appname);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
//         if (common.appname == 1)  
//             cpuFluxNS(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);
//         else if (common.appname == 2) 
//             cpuFluxPoisson(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);
//         else
            cpuFlux(f, xg, udg, odg, &app.physicsparam[0], time, 
                        numPoints, nc, ncu, nd, ncx, nco, common.appname);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
//         if (common.appname == 1)            
//             gpuFluxNS(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);
//         else if (common.appname == 2) 
//             gpuFluxPoisson(f, xg, udg, odg, &app.physicsparam[0], time, 
//                         numPoints, nc, ncu, nd, ncx, nco);
//         else
            gpuFlux(f, xg, udg, odg, &app.physicsparam[0], time, 
                        numPoints, nc, ncu, nd, ncx, nco, common.appname);                        
    }
#endif    
}

#endif
