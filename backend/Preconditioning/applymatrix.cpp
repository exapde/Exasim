/*

    This file provides functions for applying matrix operations in the context of preconditioning,
    specifically for mass matrix inversion and general matrix application on finite element meshes.

    Functions:

    void ApplyMassInv(
        cublasHandle_t handle,
        dstype *MinvR,
        dstype *Minv,
        dstype *R,
        Int curvedmesh,
        Int npe,
        Int ncu,
        Int ne,
        Int backend
    )
    ---------------------------------------------------------------------------
    Applies the inverse mass matrix to a vector R, storing the result in MinvR.
    Handles both curved and straight mesh cases, and supports different backends.
    For curved meshes, applies either a batched GEMM or a looped Gauss2Node operation.
    For straight meshes, applies a Jacobian inverse and then Gauss2Node.
    Parameters:
        - handle: cuBLAS handle for GPU operations.
        - MinvR: Output array for the result.
        - Minv: Inverse mass matrix (or Jacobian inverse for straight mesh).
        - R: Input vector.
        - curvedmesh: Flag indicating curved (1) or straight (0) mesh.
        - npe: Number of points per element.
        - ncu: Number of components per unknown.
        - ne: Number of elements.
        - backend: Backend type (CPU/GPU).

    void ApplyMatrix(
        cublasHandle_t handle,
        dstype *y,
        dstype *MassInv,
        dstype *x,
        Int npe,
        Int ncu,
        Int ne,
        Int matrixtype,
        Int curvedmesh,
        Int backend
    )
    ---------------------------------------------------------------------------
    Applies a matrix operation to vector x, storing the result in y.
    If matrixtype is 0, performs an identity operation (copy).
    Otherwise, applies the mass matrix inverse using ApplyMassInv.
    Parameters:
        - handle: cuBLAS handle for GPU operations.
        - y: Output array for the result.
        - MassInv: Mass inverse matrix.
        - x: Input vector.
        - npe: Number of points per element.
        - ncu: Number of components per unknown.
        - ne: Number of elements.
        - matrixtype: Type of matrix operation (0 for identity, otherwise mass inverse).
        - curvedmesh: Flag indicating curved (1) or straight (0) mesh.
        - backend: Backend type (CPU/GPU).
*/
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

#endif        