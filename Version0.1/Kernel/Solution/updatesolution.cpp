#ifndef __UPDATESOLUTION
#define __UPDATESOLUTION

void UpdateSolutionDIRK(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{                                   
    Int N = common.ndof1;
    Int N2 = common.npe*common.nc*common.ne2;            
    
    // update sys.u
    ArrayExtract(sys.u, sol.udg, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1, backend);                                                  
    
    // update the solution at each DIRK stage
    if (common.wave==1)
        ArrayAXPBY(sys.utmp, sol.udg, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N2, backend);                    
    else
        ArrayAXPBY(sys.utmp, sys.u, sys.utmp, common.DIRKcoeff_c[common.currentstage], one, N, backend);                    
    
    // after the last DIRK stage
    if (common.currentstage == common.tstages-1) {
       // copy utmp to udg
        if (common.wave==1)
            ArrayCopy(sol.udg, sys.utmp, N2, backend);
        else
            ArrayInsert(sol.udg, sys.utmp, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, 0, common.ne1, backend);                                                  
    }   

   // solve dw/dt = u for wave problems              
    if (common.wave==1) { 
        Int dirkStage = common.tstages;
        dstype *dirkd = &common.DIRKcoeff_d[0];    
        dstype dt = common.dt[common.currentstep];    
        dstype scalar;

        // update the source term
        switch (common.currentstage) {
            case 0:
                ArrayAXPB(sys.wsrc, sys.w, dirkd[0]/dt, zero, N, backend);
                break;
            case 1:                
                scalar = (dirkd[1]+dirkd[dirkStage+1])/dirkd[0];                
                ArrayAXPBY(sys.wsrc, sys.w, sys.wsrc, -dirkd[1]/dt, scalar, N, backend);
                                
                break;
            case 2:                
                dstype pcalar;
                scalar = (dirkd[2]/dirkd[1]);
                pcalar = dirkd[2]+dirkd[dirkStage+2]+dirkd[2*dirkStage+2] - scalar*(dirkd[1]+dirkd[dirkStage+1]);                
                ArrayAdd3Vectors(sys.wsrc, sys.w, sys.wsrc, sys.wprev, -dirkd[dirkStage+2]/dt, 
                        scalar, pcalar/dt, N, backend);                
                break;
            default:
                error("DIRK scheme not implemented yet.");
        }                 
        // dw/dt = u
        scalar = dirkd[common.currentstage*dirkStage+common.currentstage]/dt;        
        ArrayAXPBY(sys.w, sys.u, sys.wsrc, one/scalar,one/scalar, N, backend);                        
          
        // update the solution
        ArrayAXPBY(sys.wtmp, sys.w, sys.wtmp, common.DIRKcoeff_c[common.currentstage], one, N, backend);                
        
        // after the last DIRK stage
        if (common.currentstage == dirkStage-1) 
            // copy wtmp to w
            ArrayCopy(sys.w, sys.wtmp, N, backend);                                 
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
                ArrayAXPB(sys.wsrc, sys.wprev1, -common.BDFcoeff_c[1]/dt, zero, N, backend);
                break;
            case 2:                                    
                ArrayAXPBY(sys.wsrc, sys.wprev2, sys.wprev1, -common.BDFcoeff_c[1]/dt, -common.BDFcoeff_c[2]/dt, N, backend);
                break;
            case 3:                            
                ArrayAdd3Vectors(sys.wsrc, sys.wprev3, sys.wprev2, sys.wprev1, -common.BDFcoeff_c[1]/dt, 
                        -common.BDFcoeff_c[2]/dt, -common.BDFcoeff_c[3]/dt, N, backend);                
                break;
            default:
                error("BDF scheme not implemented yet.");
        }        
        dstype scalar = common.BDFcoeff_c[0]/dt;                
        ArrayAXPBY(sys.w, sys.u, sys.wsrc, one/scalar,one/scalar, N, backend);                
    }
}

void UpdateSolution(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{           
    if (common.temporalScheme==0) // DIRK
        UpdateSolutionDIRK(sol, sys, common, backend);
    else // BDF
        UpdateSolutionBDF(sol, sys, common, backend);
}

#endif
