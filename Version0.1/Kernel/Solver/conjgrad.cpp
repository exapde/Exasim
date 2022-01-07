#ifndef __CONJGRADSOLVER
#define __CONJGRADSOLVER

// Solve the linear system J(u) x = -R(u) using the conjugate gradient method
Int ConjGrad(sysstruct &sys, CDiscretization &disc, CPreconditioner& prec, Int backend)
{    
    Int maxit, it;
    //Int nc = disc.common.nc; // number of compoments of (u, q, p)
    Int ncu = disc.common.ncu;// number of compoments of (u)    
    Int npe = disc.common.npe; // number of nodes on master element    
    Int ne = disc.common.ne1; // number of elements in this subdomain 
    Int N = npe*ncu*ne;
    
    dstype alpha, itertol, rNormOld, rNormNew, nrmb, nrmr;    
    itertol = disc.common.linearSolverTol;    
    maxit = disc.common.linearSolverMaxIter;    
                
    // Set initial guess to zero
    ArraySetValue(sys.x, zero, N, backend);
    
    // evaluate the residual R(u) and set it to sys.b
    disc.evalResidual(sys.b, sys.u, backend);
    
    // copy sys.b to sys.v
    //ArrayCopy(sys.v, sys.b, N, backend);  
    
    // norm of the RHS vector b
    //nrmb = PNORM(disc.common.cublasHandle, N, sys.b, backend); 
    if (disc.common.RBcurrentdim>=0) {        
        // copy sys.b to sys.v
        ArrayCopy(sys.v, sys.b, N, backend);  
        
        // apply preconditioner to v
        prec.ApplyPreconditioner(sys.v, sys, disc, backend);
    
        // norm of the RHS vector b
        nrmb = PNORM(disc.common.cublasHandle, N, sys.v, backend);             
    }
    else
        nrmb = PNORM(disc.common.cublasHandle, N, sys.b, backend); 
    
    if  (nrmb == 0)
        nrmb = one; 
    
    //compute J(u)*x
    disc.evalMatVec(&sys.v[N], &sys.x[0], &sys.u[0], sys.b, backend);
    
    // compute (PTCparam*PTCmatrix + J(u))*x   
    if (disc.common.PTCparam>0)
        ArrayAXYPBZ(&sys.v[N], sys.PTCmatrix, sys.x, &sys.v[N], disc.common.PTCparam, one, N, backend);
    
    //compute r = -Ru - J(u)*x
    ArrayAXPBY(&sys.r[0], &sys.b[0], &sys.v[N], minusone, minusone, N, backend);    
            
    // copy the residual to sys.v = p   
    ArrayCopy(sys.v, sys.r, N, backend);    
    
    // apply preconditioner to v = p
    if (disc.common.RBcurrentdim>=0)
        prec.ApplyPreconditioner(sys.v, sys, disc, backend);
    
    // copy p to y   
    ArrayCopy(&sys.v[2*N], sys.v, N, backend);    
    
    // calculate rNormOld = (r^T y)
    PDOT(disc.common.cublasHandle, N, &sys.r[0], inc1, &sys.v[2*N], inc1, &rNormOld, backend); 
    
    // residual norm
    nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);    
    disc.common.linearSolverRelError = nrmr/nrmb;        
    if (disc.common.linearSolverRelError < itertol)
        return 0;
        
    for (it=1; it<maxit; it++)
    {        
        // do matrix-vector multiplication Ap = J(u)*p, where p = &sys.v[0]
        disc.evalMatVec(&sys.v[N], &sys.v[0], &sys.u[0], &sys.b[0], backend);
      
        // compute (PTCparam*PTCmatrix + J(u))*p  
        if (disc.common.PTCparam>0)
            ArrayAXYPBZ(&sys.v[N], sys.PTCmatrix, &sys.v[0], &sys.v[N], disc.common.PTCparam, one, N, backend);
                        
        // calculate alpha = (p^T * A * p), v[0] = p and v[N] = A*p
        PDOT(disc.common.cublasHandle, N, &sys.v[0], inc1, &sys.v[N], inc1, &alpha, backend); 
                
        // calculate alpha = (r^T * y)/ (p^T * A * p)
        alpha = rNormOld/alpha;
        
        // x = x + alpha*p
        ArrayAXPBY(&sys.x[0], &sys.x[0], &sys.v[0], one, alpha, N, backend);
        
        // r = r - alpha*A*p
        ArrayAXPBY(&sys.r[0], &sys.r[0], &sys.v[N], one, -alpha, N, backend);    
        
        // copy r to y
        ArrayCopy(&sys.v[2*N], sys.r, N, backend);    
        
        // apply preconditioner to y
        if (disc.common.RBcurrentdim>=0)
            prec.ApplyPreconditioner(&sys.v[2*N], sys, disc, backend);   
        
        // calculate rNormNew = (r^T y)
        PDOT(disc.common.cublasHandle, N, &sys.r[0], inc1, &sys.v[2*N], inc1, &rNormNew, backend); 
        
        // residual norm
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        disc.common.linearSolverRelError = nrmr/nrmb;        
        
        // check convergence
        if (disc.common.linearSolverRelError<itertol)
            break;
        
        // p = y + (rNormNew/rNormOld)*p
        ArrayAXPBY(&sys.v[0], &sys.v[2*N], &sys.v[0], one, rNormNew/rNormOld, N, backend);
        
        // update rNormOld
        rNormOld = rNormNew;
        
//         if (disc.common.mpiRank==0)
//            cout<<"CG Iteration: "<<it<<",  Relative Error: "<<disc.common.linearSolverRelError<<endl;
    }        
    
    if (disc.common.linearSolverRelError > itertol) {
        if (disc.common.mpiRank==0) {        
            printf("Warning: Conjugate gradient does not converge to the tolerance %g within % d iterations\n",itertol,maxit);
            printf("Warning: The current relative error is %g \n",disc.common.linearSolverRelError);
        }
    }
    
    return it;
}

#endif
