#ifndef __UPDATESOURCE
#define __UPDATESOURCE

void UpdateSourceDIRK(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncs = common.ncs;// number of compoments of (s)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne2; // number of elements in this subdomain 
    Int N = npe*ncs*ne;
    Int dirkStage = common.tstages;
    
    dstype *dirkd = &common.DIRKcoeff_d[0];    
    dstype dt = common.dt[common.currentstep];
    
    /* Update fc_u and fc_q */                
    dstype scalar = dirkd[common.currentstage*dirkStage+common.currentstage]/dt;
    common.dtfactor = scalar;
    // fc_u = scalar*dtcoef_u
    ArrayAXPB(app.fc_u, app.dtcoef_u, scalar, zero, common.ncu, backend);
    if (common.wave==1) 
        ArrayAXPB(app.fc_q, app.dtcoef_q, scalar, zero, common.ncq, backend);              
    
    // extract the current stage solution to res.Rq
    ArrayExtract(res.Rq, sol.udg, npe, nc, ne, 0, npe, 0, ncs, 0, ne, backend);          
             
    // update the source term
    switch (common.currentstage) {
        case 0:
            // source term for the first stage: sdg = (dirkd[0]/dt)*udg: udg
            ArrayAXPB(sol.sdg, res.Rq, dirkd[0]/dt, zero, N, backend);
            break;
        case 1:                            
            scalar = (dirkd[1]+dirkd[dirkStage+1])/dirkd[0];                
            // source term for the second stage: sdg = scalar*sdg - (dirkd[1]/dt)*udg
            ArrayAXPBY(sol.sdg, res.Rq, sol.sdg, -dirkd[1]/dt, scalar, N, backend);
            break;
        case 2:                
            dstype pcalar;
            scalar = (dirkd[2]/dirkd[1]);
            pcalar = dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2] - scalar*(dirkd[1]+dirkd[dirkStage+1]);                
            // source term for the third stage: sdg =  scalar*sdg - (dirkd[5]/dt)*udg +  (pcalar/dt)*udgprev
            ArrayAdd3Vectors(sol.sdg, res.Rq, sol.sdg, sys.udgprev, -dirkd[dirkStage+2]/dt, 
                    scalar, pcalar/dt, N, backend);                
            break;
        default:
            error("DIRK scheme not implemented yet.");
    }            
}

void UpdateSourceBDF(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncs = common.ncs;// number of compoments of (s)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne2; // number of elements in this subdomain 
    Int N = npe*ncs*ne;
    dstype dt = common.dt[common.currentstep];
    
    /* Update fc_u and fc_q */            
    dstype scalar = common.BDFcoeff_c[0]/dt;
    common.dtfactor = scalar;
    // fc_u = scalar*dtcoef_u
    ArrayAXPB(app.fc_u, app.dtcoef_u, scalar, zero, common.ncu, backend);
    if (common.wave==1) 
        ArrayAXPB(app.fc_q, app.dtcoef_q, scalar, zero, common.ncq, backend);              

    // update the source term
    switch (common.torder) {
        case 1:
            // source term for the BDF1 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev1
            ArrayAXPB(sol.sdg, sys.udgprev1, -common.BDFcoeff_c[1]/dt, zero, N, backend);
            break;
        case 2:                                    
            // source term for the BDF2 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev2 -(BDFcoeff_c[2]/dt)*udgprev1
            ArrayAXPBY(sol.sdg, sys.udgprev2, sys.udgprev1, -common.BDFcoeff_c[1]/dt, -common.BDFcoeff_c[2]/dt, N, backend);
            break;
        case 3:                            
            // source term for the BDF3 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev3 -(BDFcoeff_c[2]/dt)*udgprev2 -(BDFcoeff_c[3]/dt)*udgprev1
            ArrayAdd3Vectors(sol.sdg, sys.udgprev3, sys.udgprev2, sys.udgprev1, -common.BDFcoeff_c[1]/dt, 
                    -common.BDFcoeff_c[2]/dt, -common.BDFcoeff_c[3]/dt, N, backend);                
            break;
        default:
            error("BDF scheme not implemented yet.");
    }        
}

void UpdateSource(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
{           
    if (common.temporalScheme==0) // DIRK
        UpdateSourceDIRK(sol, sys, app, res, common, backend);
    else // BDF
        UpdateSourceBDF(sol, sys, app, res, common, backend);
}

#endif
