/*
    updatesolution.cpp

    This file contains functions for updating the solution vectors in time-dependent numerical simulations,
    supporting both DIRK (Diagonally Implicit Runge-Kutta) and BDF (Backward Differentiation Formula) temporal schemes.
    The update routines handle the propagation of solution variables (udg, wdg, etc.) through time integration stages,
    including special handling for wave problems and HDG (Hybridizable Discontinuous Galerkin) spatial schemes.

    Functions:

    - UpdateSolutionDIRK:
        Updates the solution for a DIRK time integration stage.
        Handles both standard and wave problems, updating solution vectors and temporary storage.
        After the last DIRK stage, copies the temporary solution to the main solution vector.

    - UpdateSolutionBDF:
        Updates the solution for a BDF time integration stage, specifically for wave problems.
        Computes the source term using previous solution states and BDF coefficients, then updates the solution.

    - UpdateSolution (overloaded):
        Selects the appropriate update routine (DIRK or BDF) based on the temporal scheme specified in 'common'.
        Handles both standard and wave problems.

    - UpdateSolution (extended, with appstruct, resstruct, tempstruct):
        Similar to the previous UpdateSolution, but includes additional arguments for application, residual, and temporary structures.
        Handles HDG spatial schemes by updating the solution for each element block using wEquation.
        Updates solution vectors for each DIRK stage and after the final stage.

    Key Concepts:
    - Solution vectors (udg, wdg, etc.) represent the state of the simulation at each time step.
    - DIRK and BDF are implicit time integration schemes used for advancing the solution in time.
    - The functions use utility routines (ArrayAXPBY, ArrayCopy, ArrayExtract, etc.) for efficient vector operations.
    - Special handling is provided for wave problems and HDG spatial discretizations.

    Arguments:
    - solstruct &sol: Structure containing solution vectors.
    - sysstruct &sys: Structure containing system-level temporary vectors.
    - appstruct &app: Structure containing application-specific data.
    - resstruct &res: Structure containing residual vectors.
    - tempstruct &tmp: Structure containing temporary storage for element-wise operations.
    - commonstruct &common: Structure containing common simulation parameters (e.g., time step, scheme type).
    - Int backend: Identifier for computational backend (e.g., CPU/GPU).

*/
#ifndef __UPDATESOLUTION
#define __UPDATESOLUTION

void UpdateSolutionDIRK(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{                                   
    Int N = common.ndof1;
    Int N2 = common.npe*common.nc*common.ne2;            
    
    // update sys.u
    ArrayExtract(sys.u, sol.udg, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1);                                                  
    
    // update the solution at each DIRK stage
    if (common.wave==1)
        ArrayAXPBY(sys.utmp, sol.udg, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N2);                    
    else
        ArrayAXPBY(sys.utmp, sys.u, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N);                    
    
    // after the last DIRK stage
    if (common.currentstage == common.tstages-1) {
       // copy utmp to udg
        if (common.wave==1)
            ArrayCopy(sol.udg, sys.utmp, N2);
        else
            ArrayInsert(sol.udg, sys.utmp, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1);                                                  
    }   

    // update the solution w at each DIRK stage
    if (common.ncw>0) {
        N2 = common.npe*common.ncw*common.ne2;            
        ArrayAXPBY(sys.wtmp, sol.wdg, sys.wtmp, common.DIRKcoeff_c[common.currentstage], one, N2);                
        // after the last DIRK stage
        if (common.currentstage == common.tstages-1) 
            ArrayCopy(sol.wdg, sys.wtmp, N2);                                 
    }    
}

void UpdateSolutionBDF(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{       
    Int N = common.ndof1;
    
    // solve dw/dt = u for wave problems
    if (common.wave==1) {          
        dstype dt = common.dt[common.currentstep];
        
        // update the source term
        switch (common.torder) {
            case 1:
                ArrayAXPB(sol.wsrc, sys.wprev1, -common.BDFcoeff_c[1]/dt, zero, N);
                break;
            case 2:                                    
                ArrayAXPBY(sol.wsrc, sys.wprev2, sys.wprev1, -common.BDFcoeff_c[1]/dt, -common.BDFcoeff_c[2]/dt, N);
                break;
            case 3:                            
                ArrayAdd3Vectors(sol.wsrc, sys.wprev3, sys.wprev2, sys.wprev1, -common.BDFcoeff_c[1]/dt, 
                        -common.BDFcoeff_c[2]/dt, -common.BDFcoeff_c[3]/dt, N);                
                break;
            default:
                error("BDF scheme not implemented yet.");
        }        
        dstype scalar = common.BDFcoeff_c[0]/dt;                
        ArrayAXPBY(sol.wdg, sys.u, sol.wsrc, one/scalar,one/scalar, N);                
    }
}

void UpdateSolution(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{           
    if (common.temporalScheme==0) // DIRK
        UpdateSolutionDIRK(sol, sys, common, backend);
    else // BDF
        UpdateSolutionBDF(sol, sys, common, backend);
}


void UpdateSolution(solstruct &sol, sysstruct &sys, appstruct &app, resstruct &res, tempstruct &tmp, commonstruct &common, Int backend)
{                                   
    Int N = common.ndof1;
    Int N2 = common.npe*common.nc*common.ne2;                        
    
    // update the solution at each DIRK stage
    if (common.wave==1) {
        ArrayAXPBY(sys.utmp, sol.udg, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N2);                    
    }
    else {
        ArrayExtract(res.Rq, sol.udg, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1);                                                  
        ArrayAXPBY(sys.utmp, res.Rq, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N);                    
    }

    // after the last DIRK stage
    if (common.currentstage == common.tstages-1) {
       // copy utmp to udg
        if (common.wave==1)
            ArrayCopy(sol.udg, sys.utmp, N2);
        else
            ArrayInsert(sol.udg, sys.utmp, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1);                                                  
    }   

    // update the solution w at each DIRK stage
    if (common.ncw>0) {
        if (common.spatialScheme > 0)  { // HDG
            for (Int j=0; j<common.nbe; j++) {         
                Int e1 = common.eblks[3*j]-1;
                Int e2 = common.eblks[3*j+1];
                Int ns = e2-e1;        
                Int ng = common.npe*ns;
                Int ncw = common.ncw;
                Int ncx = common.ncx;
                Int nc = common.nc;
                Int nco = common.nco;
                dstype* wdg = &tmp.tempn[0];
                dstype* xdg = &tmp.tempn[ng*ncw];
                dstype* udg = &tmp.tempn[ng*(ncw+ncx)];
                dstype* odg = &tmp.tempn[ng*(ncw+ncx+nc)];
                dstype* sdg = &tmp.tempn[ng*(ncw+ncx+nc+nco)];
                GetElemNodes(wdg, sol.wdg, common.npe, ncw, 0, ncw, e1, e2);
                GetElemNodes(xdg, sol.xdg, common.npe, ncx, 0, ncx, e1, e2);
                GetElemNodes(udg, sol.udg, common.npe, nc, 0, nc, e1, e2);
                GetElemNodes(odg, sol.odg, common.npe, nco, 0, nco, e1, e2);
                GetElemNodes(sdg, sol.wsrc, common.npe, ncw, 0, ncw, e1, e2);
                wEquation(wdg, xdg, udg, odg, sdg, tmp.tempg, app, common, ng, common.backend);
                PutElemNodes(sol.wdg, wdg, common.npe, ncw, 0, ncw, e1, e2);
            }   
        }

        N2 = common.npe*common.ncw*common.ne2;            
        ArrayAXPBY(sys.wtmp, sol.wdg, sys.wtmp, common.DIRKcoeff_c[common.currentstage], one, N2);                
        // after the last DIRK stage
        if (common.currentstage == common.tstages-1) 
            ArrayCopy(sol.wdg, sys.wtmp, N2);                                 
    }    
}

#endif
