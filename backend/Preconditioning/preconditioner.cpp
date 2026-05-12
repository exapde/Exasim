/*
    preconditioner.cpp

    This file implements the CPreconditioner class and related functions for preconditioning linear systems
    in the Exasim backend. The preconditioner is used to accelerate the convergence of iterative solvers
    by transforming the original system into a more favorable form.

    Main Components:
    ----------------
    - CPreconditioner::CPreconditioner(CDiscretization& disc, Int backend)
        Constructor: Initializes the preconditioner structure and prints info if in debug mode.

    - CPreconditioner::~CPreconditioner()
        Destructor: Frees memory allocated for the preconditioner.

    - CPreconditioner::ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend)
        Computes the initial guess and sets up the preconditioner using reduced basis techniques.
        Involves matrix-vector multiplications and small matrix inversions.

    - CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
        Applies the preconditioner to the input vector x, storing the result in disc.res.Ru.

    - CPreconditioner::ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int N, Int spatialScheme, Int backend)
        Variant of the initial guess and preconditioner computation supporting different spatial schemes.

    - ApplyBlockILU0(double* x, double* A, double* b, double *B, double *C, commonstruct& common)
        Implements block ILU(0) factorization and solve for super-element additive Schwarz preconditioner.
        Performs forward and backward solves with block matrices.

    - CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int spatialScheme, Int backend)
        Applies the preconditioner based on the spatial scheme and preconditioner type:
            - Block Jacobi
            - Elemental additive Schwarz
            - Super-element additive Schwarz with block ILU(0)

    Dependencies:
    -------------
    - setprecondstruct.cpp: Sets up the preconditioner structure.
    - applymatrix.cpp: Applies matrix operations.
    - Various BLAS-like routines (PGEMTV, PGEMTM, PGEMNV, PGEMNMStridedBached).
    - Utility functions (ArrayCopy, ArraySetValue, SubtractColumns, etc.).

    Notes:
    ------
    - The code supports both CPU and GPU backends.
    - Uses Kokkos for synchronization and parallelism.
    - Handles both sequential and non-sequential block solves for efficiency.
    - Designed for high-performance scientific computing applications.

*/
#ifndef __PRECONDITIONER
#define __PRECONDITIONER

#include "preconditioner.h"
#include "setprecondstruct.cpp"
#include "applymatrix.cpp"


// constructor
CPreconditioner::CPreconditioner(CDiscretization& disc, Int backend)
{
    mpiRank = disc.common.mpiRank;
    setprecondstruct(precond, disc, backend);    
    if ((disc.common.mpiRank==0) && (disc.common.debugMode==1)) precond.printinfo();
    if (disc.common.mpiRank == 0) printf("finish CPreconditioner constructor... \n");    
}

// destructor
CPreconditioner::~CPreconditioner()
{            
    precond.freememory(precond.backend);
    if (mpiRank==0) printf("CPreconditioner destructor: precond memory is freed successfully.\n");
}

void CPreconditioner::ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend)
{       
    // P = B + V*W^T  
    // P*W = B*W + V*W^T*W = A*W -> V = (A-B)*W
    
    // inverse(P) = inverse(B + V*W^T) = invB - U*inv(I + W^T*U)*W^T*invB = (I - U*H*W^T)*C, 
    // U = invB*V = C*V = C*(A-B)*W = C*A*W - W = Q - W, with Q = C*A*W
    // H = I + W^T*U = I + W^T C*V = I + W^T C*(A-B)*W = I + W^T invB*(A-B)*W = W^T*C*A*W = W^T*Q
    
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
    dstype *RBcoef = &disc.tmp.tempn[0];
    dstype *RBcoef_tmp = &disc.tmp.tempn[RBdim];
    dstype *H = &disc.tmp.tempg[0];
    dstype *H_tmp = &disc.tmp.tempg[RBdim*RBdim];
    
    // compute V = J(u)*W        
    //Int i = RBdim-1;
    int i = disc.common.RBremovedind - 1;
    if (disc.common.RBremovedind==0) i = RBdim-1;
    disc.evalMatVec(&precond.U[i*N], &precond.W[i*N], sys.u, sys.b, backend);
    
    //disc.evalMatVec(view_1d Jv, view_1d v, view_1d u, view_1d Ru, Int backend);
    
    /* RBcoef_tmp = V^T b */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.U, N, sys.b, inc1, 
            &zero, RBcoef_tmp, inc1, H_tmp, backend);
    
    /* H = V^T V */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.U, N, precond.U, N, &zero, 
            H, RBdim, H_tmp, backend);    
    
    /* Solve the linear system */
//     Inverse(disc.common.cublasHandle, H, H_tmp, precond.ipiv, RBdim, 1, backend);       
//     PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, H, RBdim, 
//             RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
    
    // since the matrix is very small, it is faster to invert it in CPU 
    ArrayCopy(disc.common.cublasHandle, sys.tempmem, H, RBdim*RBdim, backend);         
    Kokkos::fence();
    cpuComputeInverse(sys.tempmem, &sys.tempmem[RBdim*RBdim], sys.ipiv, RBdim);

    // RBcoef = inverse(H)*RBcoef_tmp
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, sys.tempmem, RBdim, 
            RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
        
    // compute x = -W*RBcoef
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, precond.W, N, RBcoef, 
            inc1, &zero, sys.x, inc1, backend);                                                 
}

void CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
{        
    Int N = disc.common.ndof1;        
    
    ArrayCopy(disc.common.cublasHandle, disc.res.Ru, x, N, backend);
    ApplyMatrix(disc.common.cublasHandle, x, disc.res.Minv, disc.res.Ru, disc.common.npe, disc.common.ncu, 
        disc.common.ne1, disc.common.precMatrixType, disc.common.curvedMesh, backend);                
}

void CPreconditioner::ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int N, Int spatialScheme, Int backend)
{     
    Int RBdim = disc.common.RBcurrentdim;
    dstype *RBcoef = &disc.tmp.tempn[0];
    dstype *RBcoef_tmp = &disc.tmp.tempn[RBdim];
    dstype *H = &disc.tmp.tempg[0];
    dstype *H_tmp = &disc.tmp.tempg[RBdim*RBdim];
    
    // compute V = J(u)*W        
    //Int i = RBdim-1;
    int i = disc.common.RBremovedind - 1;
    if (disc.common.RBremovedind==0) i = RBdim-1;
    disc.evalMatVec(&precond.U[i*N], &precond.W[i*N], sys.u, sys.b, spatialScheme, backend);
    
    //disc.evalMatVec(view_1d Jv, view_1d v, view_1d u, view_1d Ru, Int backend);
    
    /* RBcoef_tmp = V^T b */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.U, N, sys.b, inc1, 
            &zero, RBcoef_tmp, inc1, H_tmp, backend);
    
    /* H = V^T V */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.U, N, precond.U, N, &zero, 
            H, RBdim, H_tmp, backend);    
    
    /* Solve the linear system */
    Inverse(disc.common.cublasHandle, H, H_tmp, precond.ipiv, RBdim, 1, backend);       
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, H, RBdim, 
            RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
    
    // since the matrix is very small, it is faster to invert it in CPU 
    // ArrayCopy(disc.common.cublasHandle, sys.tempmem, H, RBdim*RBdim, backend);         
    // Kokkos::fence();
    // cpuComputeInverse(sys.tempmem, &sys.tempmem[RBdim*RBdim], sys.ipiv, RBdim);

    // // RBcoef = inverse(H)*RBcoef_tmp
    // PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, sys.tempmem, RBdim, 
    //         RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
        
    // compute x = alpha * W * RBcoef
    dstype alpha = (spatialScheme == 0) ? minusone : one;
    PGEMNV(disc.common.cublasHandle, N, RBdim, &alpha, precond.W, N, RBcoef, 
            inc1, &zero, sys.x, inc1, backend);                                                 
}

void ApplyBlockILU0(double* x, double* A, double* b, double *B, double *C, commonstruct& common) 
{    
    Int nfe = common.nfe; 
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int ncf = ncu*npf;    
    Int nse = common.nse;
    Int nfse = common.nfse;
    Int Q = ncf*nse;
    Int N = ncf*ncf*nse;    
    Int M = 2 * (nfe - 1);
    
    // Forward solve: L*y = b (unit diagonal)
    for (int i = 0; i < nfse; ++i) {
//         double *yi = &b[Q*i];
//         int k = common.Lnum_ji[0 + i*2];
//         for (int l = 0; l < k; ++l) {
//             int ptr = common.Lind_ji[l + 0*M + i*2*M];
//             int j   = common.Lind_ji[l + 1*M + i*2*M];                
//             PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, minusone, &A[ptr*N], ncf, &b[Q*j], ncf, one, yi, ncf, nse, common.backend);         
//         }
      
        double *yi = &b[Q*i];                        
        int k = common.Lnum_ji[0 + i*2];
        int flag_seq = common.Lnum_ji[1 + i*2];
        
        if (flag_seq == 1 && k > 1) {
                      
            for (int l = 0; l < k; ++l) {
              int ji = common.Lind_ji[l + 1*M + i*2*M];
              ArrayCopy(&B[Q*l], &b[Q*ji], Q);
            }                        
            
            int p1 = common.Lind_ji[0 + 0*M + i*2*M];
            PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, one, &A[N*p1], ncf, B, ncf, zero, C, ncf, nse*k, common.backend);         
            
            SubtractColumns(yi, C, Q, k);                        
        } else {          
            // Non-sequential: direct pointer access to y(:,:,j)
            for (int l = 0; l < k; ++l) {
                int ptr = common.Lind_ji[l + 0*M + i*2*M];
                int j   = common.Lind_ji[l + 1*M + i*2*M];                
                PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, minusone, &A[ptr*N], ncf, &b[Q*j], ncf, one, yi, ncf, nse, common.backend);         
            }
        }
    }

    // Backward solve: U*x = y
    for (int i = nfse-1; i >= 0; --i) {
//         double *yi = &b[Q*i];
//         int k = common.Unum_ji[0 + i*3];
//         int rstart = common.Unum_ji[1 + i*3];
//         for (int l = 0; l < k; ++l) {
//             int ptr = common.Uind_ji[l + 0*M + i*2*M];
//             int j   = common.Uind_ji[l + 1*M + i*2*M];                                
//             PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, minusone, &A[ptr*N], ncf, &x[Q*j], ncf, one, yi, ncf, nse, common.backend);         
//         }
      
        double *yi = &b[Q*i];
        int k = common.Unum_ji[0 + i*3];
        int rstart = common.Unum_ji[1 + i*3];
        int flag_seq = common.Unum_ji[2 + i*3];

        if (flag_seq == 1 && k > 1) {
                        
            for (int l = 0; l < k; ++l) {
              int ji = common.Uind_ji[l + 1*M + i*2*M];
              ArrayCopy(&B[Q*l], &x[Q*ji], Q);
            }                        
            
            int p1 = common.Uind_ji[0 + 0*M + i*2*M];
            PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, one, &A[N*p1], ncf, B, ncf, zero, C, ncf, nse*k, common.backend);         
            
            SubtractColumns(yi, C, Q, k);                  
        } else {
            // Non-sequential: direct pointer access to x(:,:,j)
            for (int l = 0; l < k; ++l) {
                int ptr = common.Uind_ji[l + 0*M + i*2*M];
                int j   = common.Uind_ji[l + 1*M + i*2*M];                                
                PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, minusone, &A[ptr*N], ncf, &x[Q*j], ncf, one, yi, ncf, nse, common.backend);         
            }
        }
        
        // Final diagonal solve for node i: x(:,n,i) = A(:,:,n,rstart) * yi(:,n) for all n
        PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, one, &A[rstart*N], ncf, yi, ncf, zero, &x[i*Q], ncf, nse, common.backend);         
    }
}

void CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int spatialScheme, Int backend)
{                
    if (spatialScheme==0) {
      Int N = disc.common.ndof1;        
      ArrayCopy(disc.common.cublasHandle, disc.res.Ru, x, N, backend);
      ApplyMatrix(disc.common.cublasHandle, x, disc.res.Minv, disc.res.Ru, disc.common.npe, disc.common.ncu, 
          disc.common.ne1, disc.common.precMatrixType, disc.common.curvedMesh, backend);                
    }
    else if (spatialScheme==1) {      
      if (disc.common.preconditioner==0) { // Block Jacobi preconditioner
        Int nf = disc.common.nf; // number of faces in this subdomain
        Int ncu = disc.common.ncu;// number of compoments of (u)
        Int npf = disc.common.npf; // number of nodes on master face           
        Int ncf = ncu*npf;  

        ArrayCopy(disc.common.cublasHandle, disc.res.Rh, x, ncf*nf, backend);

        // (ncf)  * (ncf) * nf x (ncf) * nf -> (ncf) * nf
        PGEMNMStridedBached(disc.common.cublasHandle, ncf, 1, ncf, one, disc.res.K, ncf, disc.res.Rh, ncf, zero, x, ncf, nf, backend);         
      }
      else if (disc.common.preconditioner==1) { // Elemental additive Schwarz preconditioner        
        hdgMatVec(x, disc.res.K, x, disc.res.Rh, disc.res.Rq, disc.res, disc.app, disc.mesh, disc.common, disc.tmp, disc.common.cublasHandle, backend);
      }
      else if (disc.common.preconditioner==2) { // super-element additive Schwarz preconditioner with BLIU0
        Int nf = disc.common.nf; // number of faces in this subdomain
        Int ncu = disc.common.ncu;// number of compoments of (u)
        Int npf = disc.common.npf; // number of nodes on master face           
        Int ncf = ncu*npf;          
        Int nse = disc.common.nse;
        Int nfse = disc.common.nfse;
        
        //writearray2file(disc.common.fileout + "b.bin", x, ncf*nf, backend);  
        
        GetCollumnAtIndex(disc.res.Rq, x, disc.mesh.face, ncf, nse*nfse);
        
//         writearray2file(disc.common.fileout + "b1.bin", disc.res.Rq, ncf*nse*nfse, backend);  
//         
//         writearray2file(disc.common.fileout + "K.bin", disc.res.K, ncf*ncf*nse*disc.common.nnz, backend);  
//         
        ApplyBlockILU0(disc.tmp.tempn, disc.res.K, disc.res.Rq, disc.tmp.tempg, &disc.tmp.tempn[ncf*nse*nfse], disc.common); 
        
//        writearray2file(disc.common.fileout + "b2.bin", disc.tmp.tempn, ncf*nse*nfse, backend);  
        
        ArraySetValue(x, zero, ncf*nf);
        PutCollumnAtIndexAtomicAdd(x, disc.tmp.tempn, disc.mesh.face, ncf, nse*nfse);
        
//         writearray2file(disc.common.fileout + "x.bin", x, ncf*nf, backend);  
//         
//         error("here");
        
      }
    }
}

#endif        

