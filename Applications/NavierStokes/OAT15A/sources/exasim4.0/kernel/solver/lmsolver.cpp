#ifndef __LMSOLVER
#define __LMSOLVER

Int LevenbergMarquardt(sysstruct &sys, CDiscretization& disc, CPreconditioner& prec, Int backend)
{    
    Int N = disc.common.ndof1, it = 0;
    Int RBdim = disc.common.RBcurrentdim;
    dstype lambda=0.0;
    dstype tol = 1e-6, epsilon = 1e-6;    
    dstype nrmr, nrmrold;
    
    //dstype *RBcoef = &prec.precond.RBcoef[0];    
    dstype *RBcoef = &disc.tmp.tempn[0];
    dstype *RBcoef_tmp = &disc.tmp.tempn[RBdim];
    dstype *H = &disc.tmp.tempg[0];
    dstype *H_tmp = &disc.tmp.tempg[RBdim*RBdim];
    
    // evaluate the residual vector R(u) 
    disc.evalResidual(sys.r, sys.u, backend);    

    // residual norm
    nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    
    if (disc.common.mpiRank==0)
        cout<<"LM Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
    
    // b = u
    ArrayCopy(sys.b, sys.u, N, backend);
    
    /* RBcoef = W^T u */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, prec.precond.W, N, sys.u, inc1, &zero, 
            RBcoef, inc1, H_tmp, backend);

    // u = W*RBcoef (projection of u onto the RB space)
    PGEMNV(disc.common.cublasHandle, N, RBdim, &one, prec.precond.W, N, 
            RBcoef, inc1, &zero, sys.u, inc1, backend);    
        
//     if (disc.common.mpiRank==0)
//     printArray2D(RBcoef,1,RBdim,backend);
    
    // use GNLM algorithm to solve the least-squares problem: min ||R(u(c))^T R(u(c))||
    for (it=1; it<2; it++) {                        
        nrmrold = nrmr;
        
        // dR/da = dR/du*du/da = dR/du*W = dR/du*w_i = (R(u+epsilon*w_i) - R(u))/epsilon
        /* Compute jacobian matrix using finite differences */
        for (Int i=0; i<RBdim; i++) {        
            // calculate x = W*(RBcoef + epsilon e_i) = u + epsion W e_i = u + epsilon w_i            
            ArrayAXPBY(sys.x, sys.u, &prec.precond.W[i*N], one, epsilon, N, backend);
            
            // evaluate the residual vector R(x) and set it to sys.x
            disc.evalResidual(sys.x, sys.x, backend);

            // calculate v_i = (R(x)-R(u))/epsilon    
            ArrayAXPBY(&sys.v[i*N], sys.x, sys.r, 1.0/epsilon, -1.0/epsilon, N, backend);
        }        
        
        /* RBcoef_tmp = V^T r(u) */
        PGEMTV(disc.common.cublasHandle, N, RBdim, &one, sys.v, N, sys.r, inc1, 
                &zero, RBcoef_tmp, inc1, H_tmp, backend);

//         if (disc.common.mpiRank==0)
//         printArray2D(RBcoef_tmp,1,RBdim,backend);
        
        //cout<<NORM(disc.common.cublasHandle, RBdim, RBcoef_tmp, backend)<<endl;                
        
        // check convergence        
        //if (NORM(disc.common.cublasHandle, RBdim, RBcoef_tmp, backend) <= tol)        
        //    return it;
        
        /* H = V^T V */
        PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, sys.v, N, sys.v, N, &zero, 
                H, RBdim, H_tmp, backend);    

        /* Diagonal(H) = (1.0 + lambda)*Diagonal(H) */
        ArrayMultiplyScalarDiagonal(H, one+lambda, RBdim, backend);       
        
//         if (disc.common.mpiRank==0)
//         printArray2D(H,RBdim,RBdim,backend);
        
        /* Solve the linear system */
        Inverse(disc.common.cublasHandle, H, H_tmp, prec.precond.ipiv, RBdim, 1, backend);       
        
//         if (disc.common.mpiRank==0)
//         printArray2D(H,RBdim,RBdim,backend);
        
        // RBcoef = inverse(H)*RBcoef_tmp
        PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, H, RBdim, 
                RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
                        
        // compute u = u - W*RBcoef
        PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, prec.precond.W, N, RBcoef, 
                inc1, &one, sys.u, inc1, backend);         
        
//         if (disc.common.mpiRank==0)
//         printArray2D(RBcoef,1,RBdim,backend);
        
        // evaluate the residual vector R(u) 
        disc.evalResidual(sys.r, sys.u, backend);    
        
        // residual norm
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        
        if (disc.common.mpiRank==0)
            cout<<"LM Iteration: "<<it<<",  Residual Norm: "<<nrmr<<endl;                           
        
        // update LM parameter 
        if (nrmr<nrmrold) { // 
            lambda = lambda/10.0;
            // b = u
            ArrayCopy(sys.b, sys.u, N, backend);
        }
        else {            
            lambda = lambda*10.0;
            // u = b
            ArrayCopy(sys.u, sys.b, N, backend);
        }        
    }
        
    //error("here");
    
    return it;
}

#endif

