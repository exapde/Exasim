/*
    updatesource.cpp

    This file contains functions to update the source terms for time integration schemes
    in the Exasim backend. The source terms are updated for both DIRK (Diagonally Implicit Runge-Kutta)
    and BDF (Backward Differentiation Formula) schemes, including support for differential algebraic equations.

    Functions:

    - void UpdateSourceDIRK(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
        Updates the source term for the DIRK time integration scheme.
        - Computes the time derivative contributions to the source term for each DIRK stage.
        - Updates auxiliary vectors (fc_u, fc_q) and extracts the current stage solution.
        - Handles differential algebraic equations if present.

    - void UpdateSourceBDF(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
        Updates the source term for the BDF time integration scheme.
        - Computes the time derivative contributions to the source term for BDF1, BDF2, and BDF3 orders.
        - Updates auxiliary vectors (fc_u, fc_q).

    - void UpdateSource(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
        Dispatches the source term update to either DIRK or BDF scheme based on the temporalScheme parameter.

    Notes:
    - The functions rely on helper routines (ArrayAXPB, ArrayAXPBY, ArrayAdd3Vectors, ArrayExtract) for vector operations.
    - Error handling is provided for unsupported DIRK/BDF stages.
    - The source term update is essential for accurate time integration in PDE solvers.
*/
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
    ArrayAXPB(app.fc_u, app.dtcoef_u, scalar, zero, common.ncu);
    if (common.wave==1) 
        ArrayAXPB(app.fc_q, app.dtcoef_q, scalar, zero, common.ncq);              
    
    // extract the current stage solution to res.Rq
    ArrayExtract(res.Rq, sol.udg, npe, nc, ne, 0, npe, 0, ncs, 0, ne);          
             
    // update the source term due to the time derivative
    switch (common.currentstage) {
        case 0:
            // source term for the first stage: sdg = (dirkd[0]/dt)*udg
            ArrayAXPB(sol.sdg, res.Rq, dirkd[0]/dt, zero, N);
            break;
        case 1:                            
            scalar = (dirkd[1]+dirkd[dirkStage+1])/dirkd[0];                
            // source term for the second stage: sdg = scalar*sdg - (dirkd[1]/dt)*udg
            ArrayAXPBY(sol.sdg, res.Rq, sol.sdg, -dirkd[1]/dt, scalar, N);
            break;
        case 2:                
            dstype pcalar;
            scalar = (dirkd[2]/dirkd[1]);
            pcalar = dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2] - scalar*(dirkd[1]+dirkd[dirkStage+1]);                
            // source term for the third stage: sdg =  scalar*sdg - (dirkd[5]/dt)*udg +  (pcalar/dt)*udgprev
            ArrayAdd3Vectors(sol.sdg, res.Rq, sol.sdg, sys.udgprev, -dirkd[dirkStage+2]/dt, 
                    scalar, pcalar/dt, N);                
            break;
        default:
            error("DIRK scheme not implemented yet.");
    }       
    
    // update the source term due to the time derivative for differential algebraic equations
    if (common.ncw>0) { 
        Int ncw = common.ncw;
        N = npe*ncw*ne;
        // update the source term
        switch (common.currentstage) {
            case 0:
                ArrayAXPB(sol.wsrc, sol.wdg, dirkd[0]/dt, zero, N);
                break;
            case 1:                
                scalar = (dirkd[1]+dirkd[dirkStage+1])/dirkd[0];                
                ArrayAXPBY(sol.wsrc, sol.wdg, sol.wsrc, -dirkd[1]/dt, scalar, N);                                
                break;
            case 2:                
                dstype pcalar;
                scalar = (dirkd[2]/dirkd[1]);
                pcalar = dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2] - scalar*(dirkd[1]+dirkd[dirkStage+1]);                
                ArrayAdd3Vectors(sol.wsrc, sol.wdg, sol.wsrc, sys.wprev, -dirkd[dirkStage+2]/dt, 
                        scalar, pcalar/dt, N);                
                break;
            default:
                error("DIRK scheme not implemented yet.");
        }                 
    }       
}

void UpdateSourceBDF(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, commonstruct &common, Int backend)
{   
    //Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncs = common.ncs;// number of compoments of (s)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne2; // number of elements in this subdomain 
    Int N = npe*ncs*ne;
    dstype dt = common.dt[common.currentstep];
    
    /* Update fc_u and fc_q */            
    dstype scalar = common.BDFcoeff_c[0]/dt;
    common.dtfactor = scalar;
    // fc_u = scalar*dtcoef_u
    ArrayAXPB(app.fc_u, app.dtcoef_u, scalar, zero, common.ncu);
    if (common.wave==1) 
        ArrayAXPB(app.fc_q, app.dtcoef_q, scalar, zero, common.ncq);              

    // update the source term
    switch (common.torder) {
        case 1:
            // source term for the BDF1 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev1
            ArrayAXPB(sol.sdg, sys.udgprev1, -common.BDFcoeff_c[1]/dt, zero, N);
            break;
        case 2:                                    
            // source term for the BDF2 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev2 -(BDFcoeff_c[2]/dt)*udgprev1
            ArrayAXPBY(sol.sdg, sys.udgprev2, sys.udgprev1, -common.BDFcoeff_c[1]/dt, -common.BDFcoeff_c[2]/dt, N);
            break;
        case 3:                            
            // source term for the BDF3 scheme: sdg = -(BDFcoeff_c[1]/dt)*udgprev3 -(BDFcoeff_c[2]/dt)*udgprev2 -(BDFcoeff_c[3]/dt)*udgprev1
            ArrayAdd3Vectors(sol.sdg, sys.udgprev3, sys.udgprev2, sys.udgprev1, -common.BDFcoeff_c[1]/dt, 
                    -common.BDFcoeff_c[2]/dt, -common.BDFcoeff_c[3]/dt, N);                
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
