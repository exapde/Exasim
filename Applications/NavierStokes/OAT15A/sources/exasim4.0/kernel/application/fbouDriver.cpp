#ifndef __FBOUDRIVER
#define __FBOUDRIVER

void FbouDriver(dstype *fb, dstype *xg, dstype *uhg, dstype *udg, dstype * odg, dstype *nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int ncp = common.ncp;// number of compoments of (p)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    /* 3. Compute physical boundary conditions */
#ifdef HAVE_ONETHREAD             
    if (backend==0) {
//         if (common.appname == 2)       
//             opuFbouPoisson(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
//                              ib, numPoints, nc, ncu, nd, ncx, nco);
//         else
            opuFbou(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
                             ib, numPoints, nc, ncu, nd, ncx, nco, common.appname);
    }
#endif              
#ifdef HAVE_OPENMP                
    if (backend==1) {        
//         if (common.appname == 2)       
//             cpuFbouPoisson(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
//                              ib, numPoints, nc, ncu, nd, ncx, nco);
//         else
            cpuFbou(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
                             ib, numPoints, nc, ncu, nd, ncx, nco, common.appname);
    }
#endif            
#ifdef HAVE_CUDA                     
    if (backend==2) {
//         if (common.appname == 2)       
//             gpuFbouPoisson(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
//                              ib, numPoints, nc, ncu, nd, ncx, nco);
//         else
            gpuFbou(fb, xg, udg, uhg, odg, nl, app.tau, app.uinf, app.physicsparam, time, 
                             ib, numPoints, nc, ncu, nd, ncx, nco, common.appname);
    }
#endif            
}


#endif
