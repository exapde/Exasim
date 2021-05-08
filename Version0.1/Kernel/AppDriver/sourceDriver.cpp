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

    
    /* 2. Compute physical source */
#ifdef HAVE_ONETHREAD        
    if (backend==0) {
        opuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                
    }
#endif              
#ifdef HAVE_OPENMP        
    if (backend==1) {
        cpuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }    
#endif            
#ifdef HAVE_CUDA             
    if (backend==2) {
        gpuSource(f, xg, udg, odg, wdg, app.uinf, app.physicsparam, time, 
                    numPoints, nc, ncu, nd, ncx, nco, ncw);                        
    }
#endif
    
#ifdef CHECK_NAN            
    dstype nrmf = PNORM(common.cublasHandle, numPoints*ncu, f, common.backend);
    if (isnan(nrmf) || nrmf > 1.0e14) {
        cout<<"Processor: "<<common.mpiRank<<", source norm: "<<nrmf<<endl;
        error("here");
    }
#endif    
    
    // add EUV source term to f
    // if it is simple enough, we can do that in pdemodel file
    // s = sum_n a_n * f_n(x, u, q, mu)
    // s = 0
    // for n = 1:N
    //  f_n(x, u, q, mu)
    //  s = s + a_n * f_n
    // 
    // alpha*(nabla T dot n - q0(x)) + beta * (T - T0(x)) = 0    
    /// nabla T dot n = 0
    
    // v = 0
    // rho = rhohat 
}

#endif
