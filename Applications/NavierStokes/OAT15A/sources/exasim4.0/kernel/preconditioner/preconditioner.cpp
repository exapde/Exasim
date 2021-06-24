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

void ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, precondstruct& precond, Int backend)
{
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
        
    // compute  x = C*x
    ArrayCopy(precond.y, x, N, backend);
    ApplyMatrix(disc.common.cublasHandle, x, precond.Cmat, precond.y, 
        disc.common.npe, disc.common.ncu, disc.common.ne1, 
        disc.common.precMatrixType, disc.common.curvedMesh, backend);            
    
    // disc.tmp.tempg = W^T * x
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, x, inc1, &zero, 
            disc.tmp.tempg, inc1, disc.tmp.tempn, backend);

    // disc.tmp.tempn = H*disc.tmp.tempg
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
            disc.tmp.tempg, inc1, &zero, disc.tmp.tempn, inc1, backend);    

    // compute x = x - U*disc.tmp.tempn
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, precond.U, N, disc.tmp.tempn, 
            inc1, &one, x, inc1, backend);                                     
}

void CheckLowRankPreconditioner(sysstruct& sys, CDiscretization& disc, precondstruct& precond, Int backend)
{       
    Int N = disc.common.ndof;
    Int RBdim = disc.common.RBcurrentdim;
    
    // evaluate the residual R(u) and set it to sys.b
    disc.evalResidual(sys.b, sys.u, backend);
       
    // P = B + V*W^T  
    // P*W = B*W + V*W^T*W = A*W -> V = (A-B)*W
    // P = B + (A-B)*W*W^T
    // P*W = B*W + (A-B)*W*(W^T*W) = A*W

    // inverse(P) = inverse(B + V*W^T) = invB - U*inv(I + W^T*U)*W^T*invB = (I - U*H*W^T)*C, 
    // U = invB*V = C*V = C*(A-B)*W = C*A*W - W
    // H = I + W^T*U = I + W^T C*V = I + W^T C*(A-B)*W = I + W^T invB*(A-B)*W = W^T*C*A*W
    for (Int i=0; i<RBdim; i++) {     
        // compute J(u)*W
        disc.evalMatVec(sys.r, &precond.W[i*N], sys.u, sys.b, backend);

        // compute x = (PTCparam*PTCmatrix + J(u))*W   
        if (disc.common.PTCparam>0) {
        ApplyMatrix(disc.common.cublasHandle, precond.y, sys.PTCmat, &precond.W[i*N], 
            disc.common.npe, disc.common.ncu, disc.common.ne1, disc.common.ptcMatrixType, disc.common.curvedMesh, backend);                
        ArrayAXPBY(sys.r, precond.y, sys.r, disc.common.PTCparam, one, N, backend);
        }
        
        // r = inv(P)*x = inv(P)*A*w 
        ApplyPreconditioner(sys.r, sys, disc, precond, backend);
        
        // r = r - w
        ArrayAXPBY(sys.r, sys.r, &precond.W[i*N], one, minusone, N, backend);    
        
        dstype nrmr = PNORM(disc.common.cublasHandle, N, sys.r, backend);
        
        cout<<"Check Preconditioner Inverse: "<<i<<",  Error norm: "<<nrmr<<endl;          
    }        
}

void CPreconditioner::ConstructPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend)
{       
    // P = B + V*W^T  
    // P*W = B*W + V*W^T*W = A*W -> V = (A-B)*W
    
    // inverse(P) = inverse(B + V*W^T) = invB - U*inv(I + W^T*U)*W^T*invB = (I - U*H*W^T)*C, 
    // U = invB*V = C*V = C*(A-B)*W = C*A*W - W = Q - W, with Q = C*A*W
    // H = I + W^T*U = I + W^T C*V = I + W^T C*(A-B)*W = I + W^T invB*(A-B)*W = W^T*C*A*W = W^T*Q
    
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
    
    // evaluate the residual R(u) and set it to sys.b
    disc.evalResidual(sys.b, sys.u, backend);
    
    for (Int i=0; i<RBdim; i++) {     
        // compute J(u)*W
        disc.evalMatVec(&precond.U[i*N], &precond.W[i*N], sys.u, sys.b, backend);

        // compute (PTCparam*PTCmatrix + J(u))*W   
        //ArrayAXYPBZ(&precond.U[i*N], sys.PTCmatrix, &precond.W[i*N], &precond.U[i*N], 
        //        disc.common.PTCparam, one, N, backend);
        
        if (disc.common.PTCparam>0) {
        ApplyMatrix(disc.common.cublasHandle, precond.y, sys.PTCmat, &precond.W[i*N], 
            disc.common.npe, disc.common.ncu, disc.common.ne1, disc.common.ptcMatrixType, disc.common.curvedMesh, backend);                
        ArrayAXPBY(precond.y, precond.y, &precond.U[i*N], 
                disc.common.PTCparam, one, N, backend);
        }
        else 
            ArrayCopy(precond.y, &precond.U[i*N], N, backend);
            
        // U = C*(PTCparam*PTCmatrix + J(u))*W   
        ApplyMatrix(disc.common.cublasHandle, &precond.U[i*N], precond.Cmat, precond.y, 
            disc.common.npe, disc.common.ncu, disc.common.ne1, 
            disc.common.precMatrixType, disc.common.curvedMesh, backend);                
    }
    
    /* H = W^T U = W^T C*(PTCparam*PTCmatrix + J(u))*W  */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.W, N, precond.U, N, &zero, 
            precond.H, RBdim, disc.tmp.tempg, backend);    
    
    /* Compute inverse of H */
    Inverse(disc.common.cublasHandle, precond.H, disc.tmp.tempg, precond.ipiv, RBdim, 1, backend);       
    
    // U = U - W
    ArrayAXPBY(precond.U, precond.U, precond.W, one, minusone, N*RBdim, backend);         
    
    //CheckLowRankPreconditioner(sys, disc, precond, backend);
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
    
    // evaluate the residual R(u) and set it to sys.b
    disc.evalResidual(sys.b, sys.u, backend);
    
    // residual norm
    //dstype oldnrm = PNORM(disc.common.cublasHandle, N, sys.b, backend);    
        
//     for (Int i=0; i<RBdim; i++) {     
//         // compute J(u)*W        
//         disc.evalMatVec(&sys.v[i*N], &precond.W[i*N], sys.u, sys.b, backend);
//                 
//         // U = C*(J(u))*W   
//         ApplyMatrix(disc.common.cublasHandle, &precond.U[i*N], precond.Cmat, &sys.v[i*N], 
//             disc.common.npe, disc.common.ncu, disc.common.ne1, 
//             disc.common.precMatrixType, disc.common.curvedMesh, backend);                
//     }        
    
    // compute U = J(u)*W        
    //Int i = RBdim-1;
    int i = disc.common.RBremovedind -1;
    if (disc.common.RBremovedind==0)
         i = RBdim-1;
    disc.evalMatVec(&precond.U[i*N], &precond.W[i*N], sys.u, sys.b, backend);
    
    /* RBcoef_tmp = U^T b */
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.U, N, sys.b, inc1, 
            &zero, RBcoef_tmp, inc1, H_tmp, backend);
    
    /* H = U^T U */
    PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.U, N, precond.U, N, &zero, 
            precond.H, RBdim, H_tmp, backend);    
    
    /* Solve the linear system */
    Inverse(disc.common.cublasHandle, precond.H, H_tmp, precond.ipiv, RBdim, 1, backend);       

    // RBcoef = inverse(H)*RBcoef_tmp
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
            RBcoef_tmp, inc1, &zero, RBcoef, inc1, backend);    
        
    // compute x = -W*RBcoef
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, precond.W, N, RBcoef, 
            inc1, &zero, sys.x, inc1, backend);                                                 

//     // H = W^T R
//     PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.W, N, precond.R, N, &zero, 
//             precond.H, RBdim, disc.tmp.tempg, backend);    
//     
//     /* Compute inverse of H */
//     Inverse(disc.common.cublasHandle, precond.H, disc.tmp.tempg, precond.ipiv, RBdim, 1, backend);       
//     
//     // V = W - V
//     ArrayAXPBY(precond.V, precond.W, precond.V, one, minusone, N*RBdim, backend);                 
    
//     /* H = W^T U = W^T C*(PTCparam*PTCmatrix + J(u))*W  */
//     PGEMTM(disc.common.cublasHandle, RBdim, RBdim, N, &one, precond.W, N, precond.U, N, &zero, 
//             precond.H, RBdim, disc.tmp.tempg, backend);    
//     
//     /* Compute inverse of H */
//     Inverse(disc.common.cublasHandle, precond.H, disc.tmp.tempg, precond.ipiv, RBdim, 1, backend);       
//     
//     // U = U - W
//     ArrayAXPBY(precond.U, precond.U, precond.W, one, minusone, N*RBdim, backend);                 
}

void CPreconditioner::ApplyLowRankPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
{
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
        
    // compute  x = C*x
    ArrayCopy(precond.y, x, N, backend);
    ApplyMatrix(disc.common.cublasHandle, x, precond.Cmat, precond.y, 
        disc.common.npe, disc.common.ncu, disc.common.ne1, 
        disc.common.precMatrixType, disc.common.curvedMesh, backend);            
    
    if (RBdim>0) {
    // disc.tmp.tempg = W^T * x
    //PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, x, inc1, &zero, 
    //        disc.tmp.tempg, inc1, disc.tmp.tempn, backend);
    PGEMTV2(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, x, inc1, &zero, 
            disc.tmp.tempg, inc1, disc.tmp.tempn, backend);    
    
    // disc.tmp.tempn = H*disc.tmp.tempg
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
            disc.tmp.tempg, inc1, &zero, disc.tmp.tempn, inc1, backend);    

    // compute x = x - U*disc.tmp.tempn
    PGEMNV(disc.common.cublasHandle, N, RBdim, &minusone, precond.U, N, disc.tmp.tempn, 
            inc1, &one, x, inc1, backend);   
    }
}

void CPreconditioner::ApplyBlockJacobiPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
{
    Int N = disc.common.ndof1;        
    Int RBdim = disc.common.RBcurrentdim;
    
    // compute  x = C*x
    ArrayCopy(precond.y, x, N, backend);
    ApplyMatrix(disc.common.cublasHandle, x, precond.Cmat, precond.y, 
        disc.common.npe, disc.common.ncu, disc.common.ne1, 
        disc.common.precMatrixType, disc.common.curvedMesh, backend);           

    
//     if (RBdim>0) {
//     /* RBcoef_tmp = U^T b */
//     PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.U, N, precond.y, inc1, 
//             &zero, disc.tmp.tempg, inc1, disc.tmp.tempn, backend);
// 
//     // RBcoef = inverse(H)*RBcoef_tmp
//     PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
//             disc.tmp.tempg, inc1, &zero, disc.tmp.tempn, inc1, backend);    
//     
//     dstype scalar=-1.0;
//     // compute x = x-W*RBcoef
//     PGEMNV(disc.common.cublasHandle, N, RBdim, &scalar, precond.W, N, disc.tmp.tempn, 
//             inc1, &one, x, inc1, backend);                                                 
//     }
    
//     if (RBdim>0) {
//     // disc.tmp.tempg = W^T * x
//     PGEMTV2(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, precond.y, inc1, &zero, 
//             disc.tmp.tempg, inc1, disc.tmp.tempn, backend);    
//         
//     // disc.tmp.tempn = H*disc.tmp.tempg
//     PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
//             disc.tmp.tempg, inc1, &zero, disc.tmp.tempn, inc1, backend);    
// 
//     // compute x = x + V*disc.tmp.tempn
//     PGEMNV(disc.common.cublasHandle, N, RBdim, &one, precond.V, N, disc.tmp.tempn, 
//             inc1, &one, x, inc1, backend);       
//     }
}

void CPreconditioner::ApplyReducedBasisPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
{    
    Int N = disc.common.ndof1;
    Int RBdim = disc.common.RBcurrentdim;
    dstype alpha = 0.5;
    
    // disc.tmp.tempg = W^T * x
    PGEMTV(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, x, inc1, &zero, 
            disc.tmp.tempg, inc1, disc.tmp.tempn, backend);

    // disc.tmp.tempn = H*disc.tmp.tempg = inv(W^T * A * W)*disc.tmp.tempg 
    PGEMNV(disc.common.cublasHandle, RBdim, RBdim, &one, precond.H, RBdim, 
            disc.tmp.tempg, inc1, &zero, disc.tmp.tempn, inc1, backend);    
    
    // compute z = W*disc.tmp.tempn  (RB approximation to A*z = x)
    PGEMNV(disc.common.cublasHandle, N, RBdim, &one, precond.W, N, disc.tmp.tempn, 
            inc1, &zero, precond.z, inc1, backend);           
            
    // compute y = J(u)*z
    disc.evalMatVec(precond.y, precond.z, sys.u, sys.b, backend);
        
    // compute y = (PTCparam*PTCmatrix + J(u))*z   
    ArrayAXYPBZ(precond.y, sys.PTCmatrix, precond.z, precond.y, 
            disc.common.PTCparam, one, N, backend);
//     ApplyMatrix(disc.common.cublasHandle, precond.y, sys.PTCmat, precond.z, 
//         disc.common.npe, disc.common.ncu, disc.common.ne1, disc.common.ptcMatrixType, disc.common.curvedMesh, backend);                
//     ArrayAXPBY(precond.y, precond.z, precond.y, 
//             disc.common.PTCparam, one, N, backend);
            
    // RB residual: y = x - A*z = x - y
    ArrayAXPBY(precond.y, x, precond.y, (one+alpha), minusone, N, backend);            
    
    // apply LR preconditioner to the RB residual
    this->ApplyLowRankPreconditioner(precond.y, sys, disc, backend);    
    
    // x = z (RB approximation) + y (preconditioned RB residual)
    ArrayAXPBY(x, precond.z, precond.y, (one-alpha), one, N, backend);                
}

void CPreconditioner::ApplyPreconditioner(dstype* x, sysstruct& sys, CDiscretization& disc, Int backend)
{    
    if (disc.common.preconditioner==0)        
        //this->ApplyLowRankPreconditioner(x, sys, disc, backend);
        this->ApplyBlockJacobiPreconditioner(x, sys, disc, backend);
    else
        this->ApplyReducedBasisPreconditioner(x, sys, disc, backend);
}

#endif        

