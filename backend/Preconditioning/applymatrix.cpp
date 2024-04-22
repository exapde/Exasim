#ifndef __APPLYMATRIX
#define __APPLYMATRIX

void ApplyMassInv(cublasHandle_t handle, dstype *MinvR, dstype *Minv, dstype *R, 
     Int curvedmesh, Int npe, Int ncu, Int ne, Int backend)
{    
    Int N = npe*ncu*ne;
    Int M = npe*ncu;
                
    if (curvedmesh==1) { // curved mesh
        if (backend<=1) {
            for(int i=0; i<ne; i++)
                Gauss2Node(handle, &MinvR[npe*ncu*i], &R[npe*ncu*i], &Minv[npe*npe*i], npe, npe, ncu, backend); 
        }
        else {
            ArrayGemmBatch(MinvR, Minv, R, npe, ncu, npe, ne);            
        }
    }
    else { // straight mesh
        ApplyJacInv(R, &Minv[npe*npe], M, N);        
        Gauss2Node(handle, MinvR, R, Minv, npe, npe, ncu*ne, backend); 
    }                 
}

void ApplyMatrix(cublasHandle_t handle, dstype *y, dstype *MassInv, dstype *x, 
        Int npe, Int ncu, Int ne, Int matrixtype, Int curvedmesh, Int backend)
{    
    Int N = npe*ncu*ne;
       
    if (matrixtype==0) // identity matrix      
    {
        ArrayCopy(y, x, N);
    }
    else 
      ApplyMassInv(handle, y, MassInv, x, curvedmesh, npe, ncu, ne, backend);       
}

//hdgApplyBlockJacobi(dstype *w, dstype *BE, dstype *v, commonstruct &common, cublasHandle_t handle, Int backend)

// void ApplyMass(cublasHandle_t handle, dstype *MR, dstype *Mass, dstype *R, 
//      Int curvedmesh, Int npe, Int ncu, Int ne, Int backend)
// {        
//     Int N = npe*ncu*ne;
//     Int M = npe*ncu;
//                 
//     if (curvedmesh==1) { // curved mesh
//         if (backend==2) {
//             ArrayGemmBatch(MR, Mass, R, npe, ncu, npe, ne, backend);
//         }
//         else {
//             for(Int i=0; i<ne; i++)
//                 Gauss2Node(handle, &MR[npe*ncu*i], &R[npe*ncu*i], &Mass[npe*npe*i], npe, npe, ncu, backend); 
//         }
//     }
//     else { // straight mesh
//         ApplyJac(R, R, &Mass[npe*npe], M, N, backend);        
//         Gauss2Node(handle, MR, R, Mass, npe, npe, ncu*ne, backend); 
//     }         
//     
// }
// 
// void ApplyMassInv(cublasHandle_t handle, dstype *MinvR, dstype *Minv, dstype *R, 
//      Int curvedmesh, Int npe, Int ncu, Int ne, Int backend)
// {    
//     Int N = npe*ncu*ne;
//     Int M = npe*ncu;
//                 
//     if (curvedmesh==1) { // curved mesh
//         if (backend==2) {
//             ArrayGemmBatch(MinvR, Minv, R, npe, ncu, npe, ne, backend);
//         }
//         else {
//             for(int i=0; i<ne; i++)
//                 Gauss2Node(handle, &MinvR[npe*ncu*i], &R[npe*ncu*i], &Minv[npe*npe*i], npe, npe, ncu, backend); 
//         }
//     }
//     else { // straight mesh
//         ApplyJacInv(R, R, &Minv[npe*npe], M, N, backend);        
//         Gauss2Node(handle, MinvR, R, Minv, npe, npe, ncu*ne, backend); 
//     }                 
// }
// 
// void ApplyMatrix(cublasHandle_t handle, dstype *y, dstype *Matrix, dstype *x, 
//         Int npe, Int ncu, Int ne, Int matrixtype, Int curvedmesh, Int backend)
// {
//     // compute: y = Matrix*x 
//     
//     Int N = npe*ncu*ne;
//        
//     if (matrixtype==0) // diagonal matrix      
//     {
//         ArrayAXY(y, Matrix, x, one, N, backend);   
//     }
//     else if (matrixtype==1) // mass matrix
//     { 
//         ApplyMass(handle, y, Matrix, x, curvedmesh, npe, ncu, ne, backend);   
//     }
//     else if (matrixtype==2) // inverse mass matrix
//     { 
//         ApplyMassInv(handle, y, Matrix, x, curvedmesh, npe, ncu, ne, backend);   
//     }
// }


#endif        