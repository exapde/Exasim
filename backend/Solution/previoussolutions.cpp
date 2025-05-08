#ifndef __PREVIOUSSOLUTIONS
#define __PREVIOUSSOLUTIONS

void PreviousSolutions(solstruct &sol, sysstruct &sys, commonstruct &common, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)    
    Int ncs = common.ncs;// number of compoments of (s)        
    Int npe = common.npe; // number of nodes on master element    
    //Int ne = common.ne1; // number of elements in this subdomain         
    Int ne2 = common.ne2; // number of elements in this subdomain       
    //Int N = common.ndof1;
    Int N2 = npe*common.ncw*ne2;  
        
    if (common.temporalScheme==0) {//DIRK     
        dstype scalar = one;
        for (Int j=0; j<common.tstages; j++)
            scalar = scalar - common.DIRKcoeff_c[j];
        
        /* get solution from the previous timestep */   
        ArrayExtract(sys.udgprev, sol.udg, npe, nc, ne2, 0, npe, 0, ncs, 0, ne2);   
                        
        /* utmp = (1-sum(c))*udg */  
        if (common.wave==1)
            ArrayAXPB(sys.utmp, sol.udg, scalar, zero, npe*nc*ne2);
        else         
            ArrayAXPB(sys.utmp, sys.udgprev, scalar, zero, npe*ncu*ne2);        
                        
//         if (common.wave==1) { 
//             /* get solution from the previous timestep */   
//             ArrayCopy(sys.wprev, sol.wdg, N);   
//             /* wtmp = (1-sum(c))*w */  
//             ArrayAXPB(sys.wtmp, sol.wdg, scalar, zero, N);        
//         }
        
        if (common.ncw>0) {             
            /* get solution from the previous timestep */   
            ArrayCopy(sys.wprev, sol.wdg, N2);   
            /* wtmp = (1-sum(c))*w */  
            ArrayAXPB(sys.wtmp, sol.wdg, scalar, zero, N2);        
        }
    }
    else { //BDF
        if (common.torder==1) {  // BDF1
            /* get solution from the one previous timestep */            
            ArrayExtract(sys.udgprev1, sol.udg, npe, nc, ne2, 0, npe, 0, ncs, 0, ne2);     
//             if (common.wave==1)
//                 ArrayCopy(sys.wprev1, sol.wdg, N);
            if (common.ncw>0)
                ArrayCopy(sys.wprev1, sol.wdg, N2);
        }
        else if (common.torder==2) { // BDF2          
            /* get solution from the two previous timesteps */            
            ArrayCopy(sys.udgprev1, sys.udgprev2, npe*ncs*ne2);             
            /* get solution from the previous timestep */  
            ArrayExtract(sys.udgprev2, sol.udg, npe, nc, ne2, 0, npe, 0, ncs, 0, ne2);    
            
//             if (common.wave==1) {
//                 /* get solution from the two previous timesteps */ 
//                 ArrayCopy(sys.wprev1, sys.wprev2, N);            
//                 /* get solution from the previous timestep */  
//                 ArrayCopy(sys.wprev2, sol.wdg, N); 
//             }
            if (common.ncw>0) {
                /* get solution from the two previous timesteps */ 
                ArrayCopy(sys.wprev1, sys.wprev2, N2);            
                /* get solution from the previous timestep */  
                ArrayCopy(sys.wprev2, sol.wdg, N2); 
            }
        }                    
        else if (common.torder==3) { //BDF3
            /* get solution from the three previous timesteps */            
            ArrayCopy(sys.udgprev1, sys.udgprev2, npe*ncs*ne2); 
            /* get solution from the two previous timesteps */
            ArrayCopy(sys.udgprev2, sys.udgprev3, npe*ncs*ne2); 
            /* get solution from the one previous timestep */            
            ArrayExtract(sys.udgprev3, sol.udg, npe, nc, ne2, 0, npe, 0, ncs, 0, ne2);  
//             if (common.wave==1) {
//                 /* get solution from the three previous timesteps */    
//                 ArrayCopy(sys.wprev1, sys.wprev2, N);             
//                 /* get solution from the two previous timesteps */
//                 ArrayCopy(sys.wprev2, sys.wprev3, N); 
//                 /* get solution from the one previous timestep */     
//                 ArrayCopy(sys.wprev3, sol.wdg, N); 
//             }
            if (common.ncw>0) {
                /* get solution from the three previous timesteps */    
                ArrayCopy(sys.wprev1, sys.wprev2, N2);             
                /* get solution from the two previous timesteps */
                ArrayCopy(sys.wprev2, sys.wprev3, N2); 
                /* get solution from the one previous timestep */     
                ArrayCopy(sys.wprev3, sol.wdg, N2); 
            }
        }                    
    }
}

#endif


