#ifndef __AVFIELDDRIVER
#define __AVFIELDDRIVER

void AvfieldDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int backend)
{ 
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int nce = common.nce;// number of compoments of (e)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe;
    Int ne = common.ne;
    Int numPoints = npe*ne;              

    dstype time = common.time;

    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuAvfield(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuAvfield(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuAvfield(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);         
        
    }
#endif    
}

#ifdef HAVE_ENZYME
    void AvfieldDriver(dstype *f, dstype *df, dstype *xg, dstype *udg, dstype *dudg, dstype *odg, dstype *wdg, dstype *dwdg, meshstruct &mesh, 
            masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
            commonstruct &common, Int backend)
    { 
        Int nc = common.nc; // number of compoments of (u, q)
        Int ncu = common.ncu;// number of compoments of (u)
        //Int ncq = common.ncq;// number of compoments of (q)
        Int ncw = common.ncw;// number of compoments of (w)
        Int nco = common.nco;// number of compoments of (o)
        Int nce = common.nce;// number of compoments of (e)
        Int ncx = common.ncx;// number of compoments of (xdg)        
        Int nd = common.nd;     // spatial dimension    
        Int npe = common.npe;
        Int ne = common.ne;
        Int numPoints = npe*ne;              

        dstype time = common.time;

        /* 2. Compute physical source */
    #ifdef HAVE_ONETHREAD        
        if (backend==0) {
            opuAvfieldEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                        numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);                
        }
    #endif              
    #ifdef HAVE_OPENMP        
        if (backend==1) {
            cpuAvfieldEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                        numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);                        
        }    
    #endif            
    #ifdef HAVE_CUDA             
        if (backend==2) {
            gpuAvfieldEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, common.modelnumber, 
                        numPoints, nc, ncu, nd, ncx, nco, ncw, nco, npe, ne);         
            
        }
    #endif    
}
#endif

#endif
