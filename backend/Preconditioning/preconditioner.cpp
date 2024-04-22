#ifndef __PRECONDITIONER
#define __PRECONDITIONER

#include "preconditioner.h"
#include "setprecondstruct.cpp"
#include "applymatrix.cpp"


// constructor
CPreconditioner::CPreconditioner(CDiscretization& disc, Int backend)
{
    setprecondstruct(precond, disc, backend);    
}

// destructor
CPreconditioner::~CPreconditioner()
{        
    precond.freememory(precond.cpuMemory);
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

void CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int spatialScheme, Int backend)
{                
    if (spatialScheme==0) {
      Int N = disc.common.ndof1;        
      ArrayCopy(disc.common.cublasHandle, disc.res.Ru, x, N, backend);
      ApplyMatrix(disc.common.cublasHandle, x, disc.res.Minv, disc.res.Ru, disc.common.npe, disc.common.ncu, 
          disc.common.ne1, disc.common.precMatrixType, disc.common.curvedMesh, backend);                
    }
    else if (spatialScheme==1) {
      Int nf = disc.common.nf; // number of faces in this subdomain
      Int ncu = disc.common.ncu;// number of compoments of (u)
      Int npf = disc.common.npf; // number of nodes on master face           
      Int ncf = ncu*npf;  

      ArrayCopy(disc.common.cublasHandle, disc.res.Rh, x, ncf*nf, backend);
      // (ncf)  * (ncf) * nf x (ncf) * nf -> (ncf) * nf
      PGEMNMStridedBached(disc.common.cublasHandle, ncf, 1, ncf, one, disc.res.K, ncf, disc.res.Rh, ncf, zero, x, ncf, nf, backend); 
    }
}

#endif        

