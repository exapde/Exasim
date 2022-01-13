#ifndef __FBOUDRIVER
#define __FBOUDRIVER

void FbouDriver(dstype *fb, dstype *xg, dstype *udg, dstype * odg, dstype * wdg, dstype *uhg, dstype *nl, 
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
        #ifdef HAVE_MPP
        opuFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw, app.mix);
        #else
        opuFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
        #endif
    }
#endif              
#ifdef HAVE_OPENMP                
    if (backend==1) {        
        cpuFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
#ifdef HAVE_CUDA                     
    if (backend==2) {
        gpuFbou(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
}

#ifdef HAVE_ENZYME
void FbouDriver(dstype *fb, dstype *dfb, dstype *xg, dstype *udg, dstype *dudg, dstype * odg, dstype * wdg, dstype *dwdg, dstype *uhg, dstype *duhg, dstype *nl, 
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
        opuFbouEnzyme(fb, dfb, xg, udg, dudg, odg, wdg, dwdg, uhg, duhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif              
#ifdef HAVE_OPENMP                
    if (backend==1) {        
        cpuFbouEnzyme(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
#ifdef HAVE_CUDA                     
    if (backend==2) {
        gpuFbouEnzyme(fb, xg, udg, odg, wdg, uhg, nl, app.tau, app.uinf, app.physicsparam, time, 
                         common.modelnumber, ib, numPoints, nc, ncu, nd, ncx, nco, ncw);
    }
#endif            
}

#endif


#endif
