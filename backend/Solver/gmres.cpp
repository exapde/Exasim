/*
================================================================================
GMRES Solver Implementation
================================================================================

This file implements the Generalized Minimal Residual (GMRES) iterative solver 
with support for polynomial preconditioning and restart functionality. The solver 
is designed for use with CUDA/cuBLAS and supports both classical and modified 
Gram-Schmidt orthogonalization methods.

--------------------------------------------------------------------------------
Functions:
--------------------------------------------------------------------------------

1. void CGS(cublasHandle_t handle, dstype *V, dstype *H, dstype *temp, Int N, Int m, Int backend)
    - Performs Classical Gram-Schmidt orthogonalization on Krylov vectors.
    - Computes the projection and normalization of the new vector.

2. void ApplyPoly(dstype *w, CDiscretization &disc, CPreconditioner& prec,
                        sysstruct &sys, dstype *q, dstype *p, int N, int backend)
    - Applies a polynomial preconditioner to a vector using Ritz values.
    - Handles both real and complex Ritz values for improved convergence.

3. void UpdateSolution(cublasHandle_t handle, dstype *x, dstype *y, dstype *H, dstype *s, dstype *V,
                             Int i, Int N, Int n, Int backend)
    - Updates the solution vector using the computed Krylov basis and coefficients.

4. Int GMRES(sysstruct &sys, CDiscretization &disc, CPreconditioner& prec, Int backend)
    - Main GMRES solver routine.
    - Handles initialization, polynomial preconditioning, orthogonalization, 
      convergence checks, and solution updates.
    - Supports restart and timing of key operations.

5. void ApplyPoly(dstype *w, CDiscretization &disc, CPreconditioner& prec,
                        sysstruct &sys, dstype *q, dstype *p, int N, Int spatialScheme, int backend)
    - Overloaded version of ApplyPoly with spatialScheme parameter for flexibility.

6. Int GMRES(sysstruct &sys, CDiscretization &disc, CPreconditioner& prec, Int N, Int spatialScheme, Int backend)
    - Overloaded GMRES routine supporting additional spatial scheme parameter.
    - Provides detailed timing for matrix-vector products, preconditioning, 
      orthogonalization, and solution updates.

--------------------------------------------------------------------------------
Key Concepts:
--------------------------------------------------------------------------------

- Polynomial Preconditioning: Improves convergence by transforming the system 
  using a polynomial of the matrix, based on computed Ritz values.
- Restart: Limits the number of Krylov vectors to control memory usage and 
  improve robustness.
- Orthogonalization: Supports both Classical (CGS) and Modified Gram-Schmidt (MGS).
- Timing: Measures and reports the time spent in major solver components.
- Convergence: Checks relative error against a tolerance and reports warnings 
  if not achieved.

--------------------------------------------------------------------------------
Parameters:
--------------------------------------------------------------------------------

- sysstruct: Structure containing system vectors and temporary memory.
- CDiscretization: Discretization object providing matrix-vector operations.
- CPreconditioner: Preconditioner object for applying preconditioning.
- N: Problem size (number of unknowns).
- spatialScheme: Indicates the spatial discretization scheme.
- backend: Specifies the computational backend (e.g., CPU, GPU).
- cublasHandle_t: cuBLAS handle for GPU operations.

--------------------------------------------------------------------------------
Usage:
--------------------------------------------------------------------------------

Call GMRES with appropriate system, discretization, and preconditioner objects.
Configure solver parameters (max iterations, restart, tolerance, orthogonalization method)
via the 'common' member of CDiscretization.

================================================================================
*/
#ifndef __GMRESSOLVER
#define __GMRESSOLVER 

void CGS(cublasHandle_t handle, dstype *V, dstype *H, dstype *temp, Int N, Int m, Int backend)
{
    PGEMTV(handle, N, m, &one, V, N, &V[m*N], inc1, &zero, H, inc1, temp, backend);
    PGEMNV(handle, N, m, &minusone, V, N, H, inc1, &one, &V[m*N], inc1, backend);
    PDOT(handle, N, &V[m*N], inc1, &V[m*N], inc1, &H[m], backend);
    H[m] = sqrt(H[m]);
    ArrayMultiplyScalar(&V[m*N], one/H[m], N);
}

template <typename Model>
void ApplyPoly(dstype *w, CDiscretization<Model> &disc, CPreconditioner<Model>& prec,
        sysstruct &sys, dstype *q, dstype *p, int N, int backend)
{
    int m = disc.common.ppdegree;
    dstype *sr = &sys.lam[2*m];
    dstype *si = &sys.lam[3*m];    

    ArrayCopy(q, w, N);
    ArraySetValue(w, 0.0, N);
    
    dstype a, b, a2b2;
    int i = 0;
    
    while (i<(m-1)) {
        if (si[i] == 0) {
           ArrayAXPBY(w, w, q, 1.0, 1.0/sr[i], N);                
           disc.evalMatVec(p, q, sys.u, sys.b, backend);      
           prec.ApplyPreconditioner(p, sys, disc, backend);
           ArrayAXPBY(q, q, p, 1.0, -1.0/sr[i], N);                   
           i = i + 1;            
        }        
        else {
            a = sr[i];
            b = si[i];
            a2b2 = a*a + b*b;
            disc.evalMatVec(p, q, sys.u, sys.b, backend);   
            prec.ApplyPreconditioner(p, sys, disc, backend);
            ArrayAXPBY(p, q, p, 2*a, -1.0, N);      
            ArrayAXPBY(w, w, p, 1.0, 1.0/a2b2, N);                
           if ( i < (m - 2) ) {
               disc.evalMatVec(p, p, sys.u, sys.b, backend);      
               prec.ApplyPreconditioner(p, sys, disc, backend);
               ArrayAXPBY(q, q, p, 1.0, -1.0/a2b2, N);   
           }
           i += 2; 
        }
    }
    if (si[m-1] == 0) 
        ArrayAXPBY(w, w, q, 1.0, 1.0/sr[m-1], N);                  
}

void UpdateSolution(cublasHandle_t handle, dstype *x, dstype *y, dstype *H, dstype *s, dstype *V,
        Int i, Int N, Int n, Int backend)
{
    cpuBackSolve(y, H, s, i, n);
    PGEMNV(handle, N, i+1, &one, V, N, y, inc1, &one, x, inc1, backend);
}

template <typename Model>
Int GMRES(sysstruct &sys, CDiscretization<Model> &disc, CPreconditioner<Model>& prec, Int backend)
{
    INIT_TIMING;    
    
    Int maxit, nrest, orthogMethod, n1, i, k, j = 0;
    Int ncu = disc.common.ncu;
    Int npe = disc.common.npe;
    Int ne = disc.common.ne1;
    Int N = npe*ncu*ne;
    dstype nrmb, nrmr, tol, scalar;
    tol = std::min(0.01,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
    maxit = disc.common.linearSolverMaxIter;
    nrest = disc.common.gmresRestart;
    orthogMethod = disc.common.gmresOrthogMethod;
    n1 = nrest + 1;
                
    dstype *s, *y, *cs, *sn, *H;
    s = &sys.tempmem[0];
    y = &sys.tempmem[n1];
    cs = &sys.tempmem[2*n1];
    sn = &sys.tempmem[3*n1];
    H = &sys.tempmem[4*n1];
        
    // calculate Ritz values of polynomial preconditioner
    if  (disc.common.ppdegree > 1) {
        getPoly(disc, prec, sys, sys.lam, sys.randvect, sys.ipiv, N, disc.common.ppdegree, backend);        
        dstype lmin=1.0e10, lmax=-1.0e10, smax=0.0; 
        int m = disc.common.ppdegree;
        for (int i=0; i<m; i++) {
            lmax = (sys.lam[2*m+i] > lmax) ? sys.lam[2*m+i] : lmax;
            lmin = (sys.lam[2*m+i] < lmin) ? sys.lam[2*m+i] : lmin;
            smax = (fabs(sys.lam[3*m+i]) > smax) ? fabs(sys.lam[3*m+i]) : smax;
        }        
        if (smax == 0.0) {
            for (int i=1; i<=m; i++) 
                sys.lam[i-1] = lmin + (lmax-lmin)*(0.5 - 0.5*cos((2*i-1)*M_PI/(2*m))/cos(M_PI/(2*m)));            
            LejaSort(&sys.lam[2*m], &sys.lam[3*m], sys.lam, &sys.lam[3*m], &sys.lam[4*m], m);
        }        
    }    
    
    // compute b
    //disc.evalResidual(sys.b, sys.u, backend);            
    nrmb = PNORM(disc.common.cublasHandle, N, sys.b, backend);    
    
    scalar = PNORM(disc.common.cublasHandle, N, sys.x, backend);    
    if (scalar>1e-12) {
        // r = A*x
        disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);

        // compute the new RHS vector: r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N);    

        // norm of the new RHS vector
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    }
    else {
        ArrayAXPBY(sys.r, sys.b, sys.b, minusone, zero, N);
        nrmr = nrmb;  
    }

    if (disc.common.mpiRank==0)
        cout<<"Old RHS Norm: "<<nrmb<<",  New RHS Norm: "<<nrmr<<endl; 
    
    disc.common.linearSolverTolFactor = nrmb/nrmr;
    tol = std::min(0.1,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
        
    //disc.common.ppdegree = 0;
    
    // compute r = P*r
    if (disc.common.ppdegree>1) {
        prec.ApplyPreconditioner(sys.r, sys, disc, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);    
    }
    else if (disc.common.RBcurrentdim>=0) {
        prec.ApplyPreconditioner(sys.r, sys, disc, backend);
        //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
    }
    
    // compute ||r||
    nrmb = PNORM(disc.common.cublasHandle, N, sys.r, backend);
    nrmr = nrmb;
        
    j = 0;
    while (j < maxit) {
        // v = r/||r||
        ArrayAXPB(sys.v, sys.r, one/nrmr, zero, N);
        
        // initialize s
        s[0] = nrmr;
        for (k=1; k<n1; k++)
            s[k] = zero;
        
        for (i = 0; i < nrest && j < maxit; i++, j++) {
            Int m = i + 1;
                        
            // compute v[m] = A*v[i]
            START_TIMING;                        
            disc.evalMatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, backend);            
            END_TIMING_DISC(0);    
                                                            
            START_TIMING;                                    
            // compute v[m] = P*v[m]
            if (disc.common.ppdegree>1) {
                prec.ApplyPreconditioner(&sys.v[m*N], sys, disc, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
                ApplyPoly(&sys.v[m*N], disc, prec, sys, sys.q, sys.p, N, backend);
            }
            else if (disc.common.RBcurrentdim>=0) {
                prec.ApplyPreconditioner(&sys.v[m*N], sys, disc, backend);   
                //ApplyComponentNorm(disc, sys.normcu, &sys.v[m*N], disc.res.Ru, ncu, ncu, npe, ne, backend);
            }
            END_TIMING_DISC(1);    
                                    
            // orthogonalize Krylov vectors
            
            START_TIMING;    
            if (orthogMethod == 0)
                MGS(disc.common.cublasHandle, sys.v, &H[n1*i], N, m, backend);
            else
                CGS(disc.common.cublasHandle, sys.v, &H[n1*i], y, N, m, backend);
            END_TIMING_DISC(2);    
                        
            START_TIMING;    
            // Apply Givens rotation to compute s
            cpuApplyGivensRotation(&H[n1*i], s, cs, sn, i);           
                        
            // compute relative error
            disc.common.linearSolverRelError = fabs(s[i+1])/nrmb;
            
            // check convergence and update solution: x = x + v*s
            if (disc.common.linearSolverRelError < tol) {                
                UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, i, N, n1, backend);
                return j;
            }                        
            END_TIMING_DISC(3);    
        }        
        
        // update solution: x = x + v*s
        UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, nrest-1, N, n1, backend);
               
        // compute r = A*x
        disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, backend);
                
        // r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, minusone, minusone, N);
        
        // r = P*r = P*(-b-A*x)
        if (disc.common.ppdegree>1) {
            prec.ApplyPreconditioner(sys.r, sys, disc, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
            ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, backend);
        }
        else if (disc.common.RBcurrentdim>=0) {
            prec.ApplyPreconditioner(sys.r, sys, disc, backend);            
            //ApplyComponentNorm(disc, sys.normcu, sys.r, disc.res.Ru, ncu, ncu, npe, ne, backend);
        }
        
        // compute relative error
        nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        disc.common.linearSolverRelError = nrmr/nrmb;
        
        // check convergence
        if (disc.common.linearSolverRelError < tol) {
            return j;
        }
    }       
    
    if (disc.common.linearSolverRelError > tol) {
        if (disc.common.mpiRank==0) {
            printf("Warning: GMRES(%d) does not converge to the tolerance %g within % d iterations\n",nrest,tol,maxit);
            printf("Warning: The current relative error is %g \n",disc.common.linearSolverRelError);
        }
    }
        
    return j;
}

template <typename Model>
void ApplyPoly(dstype *w, CDiscretization<Model> &disc, CPreconditioner<Model>& prec,
        sysstruct &sys, dstype *q, dstype *p, int N, Int spatialScheme, int backend)
{
    int m = disc.common.ppdegree;
    dstype *sr = &sys.lam[2*m];
    dstype *si = &sys.lam[3*m];    

    ArrayCopy(q, w, N);
    ArraySetValue(w, 0.0, N);
    
    dstype a, b, a2b2;
    int i = 0;
    
    while (i<(m-1)) {
        if (si[i] == 0) {
           ArrayAXPBY(w, w, q, 1.0, 1.0/sr[i], N);                
           disc.evalMatVec(p, q, sys.u, sys.b, spatialScheme, backend);      
           prec.ApplyPreconditioner(p, sys, disc, spatialScheme, backend);
           ArrayAXPBY(q, q, p, 1.0, -1.0/sr[i], N);                   
           i = i + 1;            
        }        
        else {
            a = sr[i];
            b = si[i];
            a2b2 = a*a + b*b;
            disc.evalMatVec(p, q, sys.u, sys.b, spatialScheme, backend);   
            prec.ApplyPreconditioner(p, sys, disc, spatialScheme, backend);
            ArrayAXPBY(p, q, p, 2*a, -1.0, N);      
            ArrayAXPBY(w, w, p, 1.0, 1.0/a2b2, N);                
           if ( i < (m - 2) ) {
               disc.evalMatVec(p, p, sys.u, sys.b, spatialScheme, backend);      
               prec.ApplyPreconditioner(p, sys, disc, spatialScheme, backend);
               ArrayAXPBY(q, q, p, 1.0, -1.0/a2b2, N);   
           }
           i += 2; 
        }
    }
    if (si[m-1] == 0) 
        ArrayAXPBY(w, w, q, 1.0, 1.0/sr[m-1], N);                  
}

template <typename Model>
Int GMRES(sysstruct &sys, CDiscretization<Model> &disc, CPreconditioner<Model>& prec, Int N, Int spatialScheme, Int backend)
{
    INIT_TIMING;    
    
    Int maxit, nrest, orthogMethod, n1, i, k, j = 0;
    dstype nrmb, nrmr, tol, scalar;
    maxit = disc.common.linearSolverMaxIter;
    nrest = disc.common.gmresRestart;
    orthogMethod = disc.common.gmresOrthogMethod;
    n1 = nrest + 1;
                
    dstype *s, *y, *cs, *sn, *H;
    s = &sys.tempmem[0];
    y = &sys.tempmem[n1];
    cs = &sys.tempmem[2*n1];
    sn = &sys.tempmem[3*n1];
    H = &sys.tempmem[4*n1];
    
    auto begin = chrono::high_resolution_clock::now(); 
    auto end = chrono::high_resolution_clock::now();
    double tm[10];
    for (int i=0; i<10; i++) tm[i]=0.0;
    
    // calculate Ritz values of polynomial preconditioner
    if  (disc.common.ppdegree > 1) {
        getPoly(disc, prec, sys, sys.lam, sys.randvect, sys.ipiv, N, disc.common.ppdegree, spatialScheme, backend);        
        dstype lmin=1.0e10, lmax=-1.0e10, smax=0.0; 
        int m = disc.common.ppdegree;
        for (int i=0; i<m; i++) {
            lmax = (sys.lam[2*m+i] > lmax) ? sys.lam[2*m+i] : lmax;
            lmin = (sys.lam[2*m+i] < lmin) ? sys.lam[2*m+i] : lmin;
            smax = (fabs(sys.lam[3*m+i]) > smax) ? fabs(sys.lam[3*m+i]) : smax;
        }        
        if (smax == 0.0) {
            for (int i=1; i<=m; i++) 
                sys.lam[i-1] = lmin + (lmax-lmin)*(0.5 - 0.5*cos((2*i-1)*M_PI/(2*m))/cos(M_PI/(2*m)));            
            LejaSort(&sys.lam[2*m], &sys.lam[3*m], sys.lam, &sys.lam[3*m], &sys.lam[4*m], m);
        }        
    }    
    
    nrmb = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.b, backend);       
    if (nrmb < disc.common.nonlinearSolverTol) {
        ArraySetValue(sys.x, zero, N);
        return 0;
    }
    scalar = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.x, backend);      
    
    dstype alpha = spatialScheme == 0 ? minusone : one;
    if (scalar>1e-12) {
        // r = A*x
        disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, spatialScheme, backend);
        
        // compute the new RHS vector: r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, alpha, minusone, N);    
        
        // norm of the new RHS vector
        nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.r, backend);
        // if (disc.common.mpiProcs>1 && disc.common.spatialScheme==1) {
        //   nrm = PNORM(disc.common.cublasHandle, disc.common.ncu*disc.common.npf*disc.common.ninterfacefaces, sys.b, backend);
        //   nrmr = sqrt(nrmr*nrmr - 0.5*nrm*nrm);
        // }                
    }
    else {
        ArrayAXPBY(sys.r, sys.b, sys.b, alpha, zero, N);
        nrmr = nrmb;  
    }

    if (disc.common.mpiRank==0)
        cout<<"Old RHS Norm: "<<nrmb<<",  New RHS Norm: "<<nrmr<<endl; 
    
    disc.common.linearSolverTolFactor = nrmb/nrmr;
    tol = min(0.1,disc.common.linearSolverTol*disc.common.linearSolverTolFactor);
            
    // compute r = P*r
    if (disc.common.ppdegree>1) {
        prec.ApplyPreconditioner(sys.r, sys, disc, spatialScheme, backend);
        ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, spatialScheme, backend);    
    }
    else if (disc.common.RBcurrentdim>=0) {
        prec.ApplyPreconditioner(sys.r, sys, disc, spatialScheme, backend);
    }
    
    // compute ||r||
    nrmb = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.r, backend);    
    nrmr = nrmb;
                
    //printf("%d %d %d %d %g\n", N, disc.common.nf, disc.common.ne, disc.common.ninterfacefaces, nrmr);
    
    j = 0;
    while (j < maxit) {
        // v = r/||r||
        ArrayAXPB(sys.v, sys.r, one/nrmr, zero, N);
        
        // initialize s
        s[0] = nrmr;
        for (k=1; k<n1; k++)
            s[k] = zero;
        
        for (i = 0; i < nrest && j < maxit; i++, j++) {
            Int m = i + 1;
            
            begin = chrono::high_resolution_clock::now();   

            // compute v[m] = A*v[i]
            disc.evalMatVec(&sys.v[m*N], &sys.v[i*N], sys.u, sys.b, spatialScheme, backend);                        
            
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
            
            end = chrono::high_resolution_clock::now();   
            tm[0] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
            
            begin = chrono::high_resolution_clock::now();   
            
            // compute v[m] = P*v[m]
            if (disc.common.ppdegree>1) {
                prec.ApplyPreconditioner(&sys.v[m*N], sys, disc, spatialScheme, backend);   
                ApplyPoly(&sys.v[m*N], disc, prec, sys, sys.q, sys.p, N, spatialScheme, backend);
            }
            else if (disc.common.RBcurrentdim>=0) {
                prec.ApplyPreconditioner(&sys.v[m*N], sys, disc, spatialScheme, backend);   
            }
            
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
            
            end = chrono::high_resolution_clock::now();   
            tm[1] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
                        
            // orthogonalize Krylov vectors
            begin = chrono::high_resolution_clock::now();   
            //CGS(disc.common.cublasHandle, sys.v, &H[n1*i], y, N, m, backend);
            if (orthogMethod == 0)
                MGS(disc.common.cublasHandle, sys.v, &H[n1*i], N, m, disc.common.ndofuhatinterface, backend);
            else
                CGS(disc.common.cublasHandle, sys.v, &H[n1*i], y, N, m, backend);
            
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
            
            end = chrono::high_resolution_clock::now();   
            tm[2] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
            
            begin = chrono::high_resolution_clock::now();   
            // Apply Givens rotation to compute s
            cpuApplyGivensRotation(&H[n1*i], s, cs, sn, i);
            
            // compute relative error
            disc.common.linearSolverRelError = fabs(s[i+1])/nrmb;           

#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
            
            //cout<<i<<"  "<<j<<"  "<<disc.common.linearSolverRelError<<endl;                         
            end = chrono::high_resolution_clock::now();   
            tm[3] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;                                
                        
            // check convergence and update solution: x = x + v*y
            if (disc.common.linearSolverRelError < tol) {                
              UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, i, N, n1, backend);                
              if (disc.common.mpiRank==0) {
                printf("Matrix-vector product time: %g miliseconds\n", tm[0]);
                printf("Applying preconditioner time: %g miliseconds\n", tm[1]);
                printf("Orthgonalization time: %g miliseconds\n", tm[2]);    
                printf("Solution update time: %g miliseconds\n", tm[3]);      
              }
              return j;
            }                        
        }        
                
        // update solution: x = x + v*s
        UpdateSolution(disc.common.cublasHandle, sys.x, y, H, s, sys.v, i-1, N, n1, backend);
               
        // compute r = A*x
        disc.evalMatVec(sys.r, sys.x, sys.u, sys.b, spatialScheme, backend);
                
        // r = -b - A*x
        ArrayAXPBY(sys.r, sys.b, sys.r, alpha, minusone, N);
        
        // r = P*r = P*(-b-A*x)
        if (disc.common.ppdegree>1) {
            prec.ApplyPreconditioner(sys.r, sys, disc, spatialScheme, backend);            
            ApplyPoly(sys.r, disc, prec, sys, sys.q, sys.p, N, spatialScheme, backend);
        }
        else if (disc.common.RBcurrentdim>=0) {
            prec.ApplyPreconditioner(sys.r, sys, disc, spatialScheme, backend);            
        }
        
        // compute relative error
        nrmr = PNORM(disc.common.cublasHandle, N, disc.common.ndofuhatinterface, sys.r, backend);
        disc.common.linearSolverRelError = nrmr/nrmb;
        
        // check convergence
        if (disc.common.linearSolverRelError < tol) {
          if (disc.common.mpiRank==0) {
            printf("Matrix-vector product time: %g miliseconds\n", tm[0]);
            printf("Applying preconditioner time: %g miliseconds\n", tm[1]);
            printf("Orthgonalization time: %g miliseconds\n", tm[2]);    
            printf("Solution update time: %g miliseconds\n", tm[3]);          
          }
            return j;
        }
    }       
    
    if (disc.common.linearSolverRelError > tol) {
        if (disc.common.mpiRank==0) {
            printf("Matrix-vector product time: %g miliseconds\n", tm[0]);
            printf("Applying preconditioner time: %g miliseconds\n", tm[1]);
            printf("Orthgonalization time: %g miliseconds\n", tm[2]);    
            printf("Solution update time: %g miliseconds\n", tm[3]);          
            printf("Warning: GMRES(%d) does not converge to the tolerance %g within % d iterations\n",nrest,tol,maxit);
            printf("Warning: The current relative error is %g \n",disc.common.linearSolverRelError);
        }
    }
        
    return j;
}

#endif
