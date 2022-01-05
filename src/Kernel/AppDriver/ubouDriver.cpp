#ifndef __UBOUDRIVER
#define __UBOUDRIVER

void UbouDriver(dstype *fb, dstype *xg, dstype *udg, dstype * odg, dstype * wdg, dstype *uhg, dstype *nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    /* 3. Compute physical boundary conditions */
#ifdef HAVE_ONETHREAD             
    if (backend==0) {
        opuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif              
#ifdef HAVE_OPENMP                
    if (backend==1) {        
        cpuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
#ifdef HAVE_CUDA                     
    if (backend==2) {
        gpuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
}

#ifdef HAVE_ENZYME
void UbouDriver(dstype *fb, dstype *dfb, dstype *xg, dstype *udg, dstype *dudg, 
        dstype * odg, dstype * wdg, dstype *dwdg, dstype *uhg, dstype *nl, 
        meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int ngf, Int f1, Int f2, Int ib, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int numPoints = ngf*(f2-f1);
    dstype time = common.time;    

    /* 3. Compute physical boundary conditions */
#ifdef HAVE_ONETHREAD             
    if (backend==0) {
        // opuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                        //  common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
        opuUbouEnzyme(fb, dfb, xg, udg, dudg, odg, wdg, dwdg, uhg, nl, app.tau, app.uinf, 
                app.physicsparam, time, common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif              
#ifdef HAVE_OPENMP                
    if (backend==1) {        
        // cpuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                        //  common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
        cpuUbouEnzyme(fb, dfb, xg, udg, dudg, odg, wdg, dwdg, uhg, nl, app.tau, app.uinf, 
                app.physicsparam, time, common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
#ifdef HAVE_CUDA                     
    if (backend==2) {
        // gpuUbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                        //  common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
        gpuUbouEnzyme(fb, dfb, xg, udg, dudg, odg, wdg, dwdg, uhg, nl, app.tau, app.uinf, 
                app.physicsparam, time, common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
}
#endif

#endif
