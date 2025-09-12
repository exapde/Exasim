/*
    massinv.cpp

    This file contains functions for computing and applying the inverse of the mass matrix in finite element discretizations.

    Functions:

    1. void ComputeMinv(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
                        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
        - Computes the inverse of the mass matrix (Minv) for either straight or curved meshes.
        - Allocates memory for mass and inverse matrices, workspace, and pivot arrays.
        - For straight meshes, computes a single master element mass matrix and its inverse.
        - For curved meshes, computes element-wise mass matrices and their inverses.
        - Uses geometric and shape information to assemble matrices.
        - Handles 1D, 2D, and 3D cases with appropriate geometric computations.
        - Frees allocated memory at the end.

    2. void ApplyMinv(dstype* MinvR, dstype* Minv, dstype* R, dstype scalar, Int curvedmesh, Int npe, Int ncr, Int e1, Int e2)
        - Applies the inverse mass matrix to a residual vector R, storing the result in MinvR.
        - Handles both curved and straight mesh cases.
        - Performs batched matrix multiplication or direct application depending on mesh type.
        - Scales the result by a given scalar.

    Notes:
    - The functions rely on several utility routines for memory management, matrix operations, and geometric computations.
    - The code is designed to work with both CPU and GPU backends, as indicated by the 'backend' parameter.
    - Assumes existence of supporting data structures and utility functions (e.g., TemplateMalloc, ArrayCopy, Inverse, etc.).
*/
#ifndef __MASSINV
#define __MASSINV

void ComputeMinv(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle,  Int backend)
{   
    dstype *work=NULL;  
    Int *ipiv=NULL;
    
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int nd = common.nd;     // spatial dimension    
    Int npe = common.npe; // number of nodes on master element
    Int nge = common.nge; // number of gauss points on master element    
    Int ne = common.ne; // number of elements in this subdomain 
    Int nbe = common.nbe; // number of blocks for elements   
    Int neb = common.neb; // maximum number of elements per block
            
    if (common.curvedMesh==0) { //straight mesh                        
        TemplateMalloc(&res.Mass, npe*npe+ne, backend);        
        TemplateMalloc(&res.Minv, npe*npe+ne, backend);
        TemplateMalloc(&work, npe*npe, backend);        
        TemplateMalloc(&ipiv, npe+1, backend);         
        res.szMass = npe*npe+ne;
        res.szMinv = npe*npe+ne;
        
        // mass inverse for the master element
        Gauss2Node(handle, &res.Minv[0], master.shapegt, master.shapegw, nge, npe, npe, backend);        
        ArrayCopy(res.Mass, res.Minv, npe*npe);
        Inverse(handle, &res.Minv[0], work, ipiv, npe, 1, backend);        
    }
    else { // curved mesh        
        TemplateMalloc(&res.Mass, npe*npe*ne, backend);
        TemplateMalloc(&res.Minv, npe*npe*ne, backend);
        TemplateMalloc(&work, max(nge*npe*neb,npe*npe*neb), backend);      
        TemplateMalloc(&ipiv, npe+1, backend);      
        res.szMass = npe*npe*ne;
        res.szMinv = npe*npe*ne;
    }         
        
    Int mm = 0;
    for (Int j=0; j<nbe; j++) // for each block of elements
    {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];
        Int ns = e2-e1;        
        Int nga = nge*ns;
        //Int nn =  npe*ns;         
        Int n0 = 0;                                 // xg
        Int n1 = nga*ncx;                           // Xx
        Int n2 = nga*(ncx+nd*nd);                   // jac
        Int n3 = nga*(ncx+nd*nd+1);                 // Jg, ug

        GetElemNodes(tmp.tempn, sol.xdg, npe, ncx, 0, ncx, e1, e2);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapegt, nge, npe, ns*ncx, backend);

        if (nd==1) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ns*nd, backend);               
            ElemGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ns*nd, backend);        
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapegt[2*nge*npe], nge, npe, ns*nd, backend);        
            ElemGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n1+2*nga], &tmp.tempg[n1+nga], &tmp.tempg[n1+3*nga],
            &tmp.tempg[n3], &tmp.tempg[n3+nga], &tmp.tempg[n3+2*nga], &tmp.tempg[n3+3*nga], nga);
        }
        else if (nd==3) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ns*nd, backend);        
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapegt[2*nge*npe], nge, npe, ns*nd, backend);        
            Node2Gauss(handle, &tmp.tempg[n3+2*nga*nd], tmp.tempn, &master.shapegt[3*nge*npe], nge, npe, ns*nd, backend);        
            ElemGeom3D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n1+3*nga], &tmp.tempg[n1+6*nga], 
            &tmp.tempg[n1+nga], &tmp.tempg[n1+4*nga], &tmp.tempg[n1+7*nga], 
            &tmp.tempg[n1+2*nga], &tmp.tempg[n1+5*nga], &tmp.tempg[n1+8*nga], 
            &tmp.tempg[n3], &tmp.tempg[n3+nga], &tmp.tempg[n3+2*nga], 
            &tmp.tempg[n3+3*nga], &tmp.tempg[n3+4*nga], &tmp.tempg[n3+5*nga],
            &tmp.tempg[n3+6*nga], &tmp.tempg[n3+7*nga], &tmp.tempg[n3+8*nga], nga);
        }
                
        if (common.curvedMesh==0) {
            ArrayExtract(&res.Minv[npe*npe+mm], &tmp.tempg[n2], 
                            nge, ns, 1, 0, 1, 0, ns, 0, 1);
            
            ArrayCopy(&res.Mass[npe*npe+mm], &res.Minv[npe*npe+mm], ns);
        }
        else {       
            ShapJac(work, master.shapegt, &tmp.tempg[n2], nge, npe, ne); // fixed bug here                           
            Gauss2Node(handle, &res.Minv[npe*npe*mm], work, master.shapegw, nge, npe, npe*ns, backend);           
            ArrayCopy(&res.Mass[npe*npe*mm], &res.Minv[npe*npe*mm], npe*npe*ns);
            Inverse(handle, &res.Minv[npe*npe*mm], work, ipiv, npe, ns, backend);
        }
        mm = mm+ns;
    }                
    
    TemplateFree(work, backend);   
    TemplateFree(ipiv, backend);        
}

void ApplyMinv(dstype* MinvR, dstype* Minv, dstype* R, dstype scalar, Int curvedmesh, Int npe, Int ncr, Int e1, Int e2)
{    
    Int ns = e2-e1;
    Int N = npe*ncr*ns;        
    
    if (curvedmesh==1) {
      ArrayGemmBatch1(MinvR, &Minv[npe*npe*e1], R, npe, ncr, npe, ns);    
    }
    else {        
        ApplyJacInv(R, &Minv[npe*npe+e1], npe*ncr, npe*ncr*ns);        
        ArrayMatrixMultiplication1(MinvR, Minv, R, npe, ncr*ns, npe); 
    }

    ArrayMultiplyScalar(MinvR, scalar, N);        
}


#endif  

