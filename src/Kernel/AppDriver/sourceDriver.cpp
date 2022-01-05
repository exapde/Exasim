#ifndef __SOURCEDRIVER
#define __SOURCEDRIVER

void SourceDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
        masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{ 
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int numPoints = nge*(e2-e1);              
    dstype time = common.time;

    //function app.uinf = (app.stgdata, app.stgparam);
//     if (common.gitm==1) {
//         // compute app.uinf
//     }
            
    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }
#endif    
}

#ifdef HAVE_ENZYME   

void SourceDriver(dstype *f, dstype *df, dstype *xg, dstype *udg, dstype *dudg, dstype *odg, dstype *wdg, 
        dstype *dwdg, meshstruct &mesh, masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
        commonstruct &common, Int nge, Int e1, Int e2, Int backend)
{
    Int nc = common.nc; // number of compoments of (u, q)
    Int ncu = common.ncu;// number of compoments of (u)
    //Int ncq = common.ncq;// number of compoments of (q)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int numPoints = nge*(e2-e1);              
    dstype time = common.time;
                
    ArraySetValue(df, numPoints*ncu, 0.0, backend);
    
    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuSourceEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuSourceEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuSourceEnzyme(f, df, xg, udg, dudg, odg, wdg, dwdg, app.uinf, app.physicsparam, time, 
                    common.modelnumber, numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif           
}
#endif        

#endif
