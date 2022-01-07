#ifndef __MATVEC
#define __MATVEC

void MatVec(dstype *w, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master,
      meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, dstype *v, 
      dstype *u, dstype *Ru, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int N = npe*ncu*ne;
    Int order = common.matvecOrder;
    dstype epsilon = common.matvecTol;
//    dstype nrm;
    
    if (order==1) {
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);
                            
        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u))/epsilon    
        ArrayAXPBY(w, res.Ru, Ru, 1.0/epsilon, -1.0/epsilon, N, backend);
        
//         nrm = PNORM(common.cublasHandle, N, res.Ru, backend);
//         if (common.mpiRank==0) {
//             printf("The norm R is %g\n",nrm);  
//         }
//         nrm = PNORM(common.cublasHandle, N, Ru, backend);
//         if (common.mpiRank==0) {
//             printf("The norm R is %g\n",nrm);  
//         }
//         nrm = PNORM(common.cublasHandle, N, v, backend);
//         if (common.mpiRank==0) {
//             printf("The norm v is %g\n",nrm);  
//         }
//         nrm = PNORM(common.cublasHandle, N, sol.udg, backend);
//         if (common.mpiRank==0) {
//             printf("The norm udg is %g\n",nrm);  
//         }
//         nrm = PNORM(common.cublasHandle, N, w, backend);
//         if (common.mpiRank==0) {
//             printf("The norm w is %g\n",nrm);  
//         }                    
//         if (nrm==0) {
//             nrm = PNORM(disc.common.cublasHandle, N, w, backend);
//             if (disc.common.mpiRank==0) {
//                 printf("The norm w is %g\n",nrm);  
//             }            
//         }
        
    }
    else if (order==2) {
        // calculate w = u - epsilon*v
        ArrayAXPBY(w, u, v, 1.0, -epsilon, N, backend);

        // insert (u-epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u-epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // copy res.Ru to Ru
        ArrayCopy(Ru, res.Ru, N, backend);
        
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N, backend);

        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne, backend);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);
        
        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u-epsilon*v))/(2*epsilon)    
        ArrayAXPBY(w, res.Ru, Ru, 0.5/epsilon, -0.5/epsilon, N, backend);
    }
    else
        error("Matrix-vector multiplication order is not implemented");
}


#endif


