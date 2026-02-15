/*
    geometry.cpp

    This file contains functions for computing geometric quantities at element and face Gauss points
    for finite element discretizations in 1D, 2D, and 3D. The geometric data is used for numerical integration
    and flux computations in high-order methods.

    Main Functions:
    ---------------
    - ElemGeomBlock: Computes geometric quantities (Jacobian, coordinates, etc.) at Gauss points for a block of elements.
    - ElemGeom: Allocates memory and computes geometric quantities for all elements in the domain.
    - FaceGeomBlock: Computes geometric quantities (Jacobian, normal vectors, etc.) at Gauss points for a block of faces.
    - FaceGeom: Allocates memory and computes geometric quantities for all faces in the domain.
    - ElemFaceGeomBlock: Computes geometric quantities at Gauss points on all faces of a block of elements.
    - ElemFaceGeom: Allocates memory and computes geometric quantities for all element faces in the domain.

    Key Concepts:
    -------------
    - Uses master element shape functions and mappings to compute coordinates and derivatives at Gauss points.
    - Handles different spatial dimensions (nd = 1, 2, 3) with specialized routines for each case.
    - Supports block-wise computation for parallelization and efficient memory access.
    - Uses temporary storage for intermediate results and copies final results to solution structures.
    - Integrates with CUDA/cuBLAS for GPU acceleration via the cublasHandle_t parameter.

    Parameters:
    -----------
    - solstruct, masterstruct, meshstruct, tempstruct, commonstruct: Data structures holding solution, master element,
      mesh, temporary arrays, and common parameters.
    - cublasHandle_t handle: cuBLAS handle for GPU computations.
    - Int backend: Specifies the computational backend (CPU/GPU).
    - Int nd, npe, nge, ncx, e1, e2, f1, f2, nfe, ne, ngf, npf: Various mesh and discretization parameters.

    Notes:
    ------
    - The functions assume that memory allocation and initialization of input structures are handled elsewhere.
    - ArrayCopy and TemplateMalloc are utility functions for memory management.
    - Specialized routines (ElemGeom1D, ElemGeom2D, ElemGeom3D, FaceGeom1D, FaceGeom2D, FaceGeom3D) are called
      depending on the spatial dimension.
    - Normal vectors on boundaries are adjusted to ensure correct orientation.
*/
#ifndef __GEOMETRY
#define __GEOMETRY

void ElemGeomBlock(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, 
        commonstruct &common, cublasHandle_t handle, Int nd, 
        Int npe, Int nge, Int ncx, Int e1, Int e2, Int backend)
{            
    Int ne = e2-e1;
    //Int nn =  npe*ne; 
    Int nga = nge*ne;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // Xx
    Int n2 = nga*(ncx+nd*nd);                   // jac
    Int n3 = nga*(ncx+nd*nd+1);                 // Jg, ug
    
    GetElemNodes(tmp.tempn, sol.xdg, npe, ncx, 0, ncx, e1, e2);    
    Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapegt, nge, npe, ne*ncx, backend);
    
    if (nd==1) {
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ne*nd, backend);        
        ElemGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
    }
    else if (nd==2){
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ne*nd, backend);                
        Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapegt[2*nge*npe], nge, npe, ne*nd, backend);        
        ElemGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n1+2*nga], &tmp.tempg[n1+nga], &tmp.tempg[n1+3*nga],
        &tmp.tempg[n3], &tmp.tempg[n3+nga], &tmp.tempg[n3+2*nga], &tmp.tempg[n3+3*nga], nga);                    
    }
    else if (nd==3) {
       Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ne*nd, backend);        
       Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapegt[2*nge*npe], nge, npe, ne*nd, backend);        
       Node2Gauss(handle, &tmp.tempg[n3+2*nga*nd], tmp.tempn, &master.shapegt[3*nge*npe], nge, npe, ne*nd, backend);        
       ElemGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n1+3*nga], &tmp.tempg[n1+6*nga], 
                &tmp.tempg[n1+nga], &tmp.tempg[n1+4*nga], &tmp.tempg[n1+7*nga], 
                &tmp.tempg[n1+2*nga], &tmp.tempg[n1+5*nga], &tmp.tempg[n1+8*nga], 
                &tmp.tempg[n3], &tmp.tempg[n3+nga], &tmp.tempg[n3+2*nga], 
                &tmp.tempg[n3+3*nga], &tmp.tempg[n3+4*nga], &tmp.tempg[n3+5*nga],
                &tmp.tempg[n3+6*nga], &tmp.tempg[n3+7*nga], &tmp.tempg[n3+8*nga], nga);
    }
}

void ElemGeom(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, 
        commonstruct &common, cublasHandle_t handle, Int backend)
{    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element    
    Int ne = common.ne; // number of elements in this subdomain 
    
    TemplateMalloc(&sol.elemg, nge*ne*(ncx+nd*nd+1), backend);  
    sol.szelemg = nge*ne*(ncx+nd*nd+1);
    for (Int j=0; j<common.nbe; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];    
        ElemGeomBlock(sol, master, mesh, tmp, common, common.cublasHandle, nd, npe, nge, ncx, e1, e2, backend);                
        ArrayCopy(&sol.elemg[nge*e1*(ncx+nd*nd+1)], &tmp.tempg[0], nge*(e2-e1)*(ncx+nd*nd+1));   
    }                     
}

void FaceGeomBlock(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, 
        commonstruct &common, cublasHandle_t handle, Int f1, Int f2, Int backend)
{            
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face              

    Int nf = f2-f1;
    //Int nn =  npf*nf; 
    Int nga = ngf*nf;   
    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg

    GetFaceNodes(tmp.tempn, sol.xdg, mesh.facecon, npf, ncx, npe, ncx, f1, f2, 1);            
    Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfgt, ngf, npf, nf*ncx, backend);    
    // (jac, nlg) = (tmp.tempg[n2], tmp.tempg[n1]) at gauss points on face
    if (nd==1) {
        FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);    
        // Change the direction of the normal vector on the left boundary to make it outward
        FixNormal1D(&tmp.tempg[n1], &mesh.facecon[2*f1], nga);      
    }
    else if (nd==2){
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfgt[ngf*npf], ngf, npf, nf*nd, backend);                
        FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
    }
    else if (nd==3) {
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfgt[ngf*npf], ngf, npf, nf*nd, backend);                     
        Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfgt[2*ngf*npf], ngf, npf, nf*nd, backend);                
        FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
    }
}

void FaceGeom(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        cublasHandle_t handle, Int backend)
{           
    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int ngf = common.ngf; // number of gauss points on master face                      
    Int nbf = common.nbf;
    Int nf = common.fblks[3*(nbf-1)+1];    
    
    TemplateMalloc(&sol.faceg, ngf*nf*(ncx+nd+1), backend);        
    sol.szfaceg = ngf*nf*(ncx+nd+1);
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];       
        //Int ib = common.fblks[3*j+2];    
        FaceGeomBlock(sol, master, mesh, tmp, common, handle, f1, f2, backend);                
        ArrayCopy(&sol.faceg[ngf*f1*(ncx+nd+1)], tmp.tempg, ngf*(f2-f1)*(ncx+nd+1));  
    }        
}

void ElemFaceGeomBlock(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, 
        commonstruct &common, cublasHandle_t handle, Int e1, Int e2, Int backend)
{            
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master face           
    Int ngf = common.ngf; // number of gauss poInts on master face              
    Int nfe = common.nfe; // number of faces per element

    Int ne = e2-e1;
    Int nga = ngf*nfe*ne;   

    Int n0 = 0;                                 // xg
    Int n1 = nga*ncx;                           // nlg
    Int n2 = nga*(ncx+nd);                      // jac
    Int n3 = nga*(ncx+nd+1);                    // Jg, uhg

    GetElementFaceNodes(tmp.tempn, sol.xdg, mesh.perm, npf*nfe, ncx, npe, ncx, e1, e2);            
    Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapfgt, ngf, npf, nfe*ne*ncx, backend);    
    
    // (jac, nlg) = (tmp.tempg[n2], tmp.tempg[n1]) at gauss points on face
    if (nd==1) {
        FaceGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);      
        FixNormal1D(&tmp.tempg[n1], nga);    
    }
    else if (nd==2){
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfgt[ngf*npf], ngf, npf, nfe*ne*nd, backend);                
        FaceGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
    }
    else if (nd==3) {
        Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapfgt[ngf*npf], ngf, npf, nfe*ne*nd, backend);                     
        Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapfgt[2*ngf*npf], ngf, npf, nfe*ne*nd, backend);                
        FaceGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
    }
}

void ElemFaceGeom(solstruct &sol, masterstruct &master, meshstruct &mesh, tempstruct &tmp, 
        commonstruct &common, cublasHandle_t handle, Int backend)
{    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension        
    Int ngf = common.ngf; // number of gauss poInts on master face    
    Int nfe = common.nfe; // number of faces per element
    Int ne = common.ne; // number of elements in this subdomain 
    
    TemplateMalloc(&sol.elemfaceg, ngf*nfe*ne*(ncx+nd+1), backend);  // fixed bug here
    sol.szelemfaceg = ngf*nfe*ne*(ncx+nd+1);
    for (Int j=0; j<common.nbe; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];    
        ElemFaceGeomBlock(sol, master, mesh, tmp, common, common.cublasHandle, e1, e2, backend);                
        ArrayCopy(&sol.elemfaceg[ngf*nfe*e1*(ncx+nd+1)], &tmp.tempg[0], ngf*nfe*(e2-e1)*(ncx+nd+1));  // fixed bug here       
    }                     
}

#endif


