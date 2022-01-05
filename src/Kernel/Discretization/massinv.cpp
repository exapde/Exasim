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
        if ((common.PTCparam>0) && (common.ptcMatrixType>0))
            TemplateMalloc(&res.Mass, npe*npe+ne, backend);        
        TemplateMalloc(&res.Minv, npe*npe+ne, backend);
        TemplateMalloc(&work, npe*npe, backend);        
        TemplateMalloc(&ipiv, npe+1, backend);           
        
        // mass inverse for the master element
        Gauss2Node(handle, &res.Minv[0], master.shapegt, master.shapegw, nge, npe, npe, backend);
        if ((common.PTCparam>0) && (common.ptcMatrixType>0))
            ArrayCopy(res.Mass, res.Minv, npe*npe, backend);
        Inverse(handle, &res.Minv[0], work, ipiv, npe, 1, backend);        
    }
    else { // curved mesh
        if ((common.PTCparam>0) && (common.ptcMatrixType>0))
            TemplateMalloc(&res.Mass, npe*npe*ne, backend);
        TemplateMalloc(&res.Minv, npe*npe*ne, backend);
        TemplateMalloc(&work, max(nge*npe*neb,npe*npe*neb), backend);      
        TemplateMalloc(&ipiv, npe+1, backend);           
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

        GetElemNodes(tmp.tempn, sol.xdg, npe, ncx, 0, ncx, e1, e2, backend);
        Node2Gauss(handle, &tmp.tempg[n0], tmp.tempn, master.shapegt, nge, npe, ns*ncx, backend);

        if (nd==1) {
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ns*nd, backend);               
            ElemGeom1D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n3], nga, backend);
        }
        else if (nd==2){
            Node2Gauss(handle, &tmp.tempg[n3], tmp.tempn, &master.shapegt[nge*npe], nge, npe, ns*nd, backend);        
            Node2Gauss(handle, &tmp.tempg[n3+nga*nd], tmp.tempn, &master.shapegt[2*nge*npe], nge, npe, ns*nd, backend);        
            ElemGeom2D(&tmp.tempg[n2], &tmp.tempg[n1], &tmp.tempg[n1+2*nga], &tmp.tempg[n1+nga], &tmp.tempg[n1+3*nga],
            &tmp.tempg[n3], &tmp.tempg[n3+nga], &tmp.tempg[n3+2*nga], &tmp.tempg[n3+3*nga], nga, backend);
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
            &tmp.tempg[n3+6*nga], &tmp.tempg[n3+7*nga], &tmp.tempg[n3+8*nga], nga, backend);
        }
        
        if (common.curvedMesh==0) {
            ArrayExtract(&res.Minv[npe*npe+mm], &tmp.tempg[n2], 
                            nge, ns, 1, 0, 1, 0, ns, 0, 1, backend);
            if ((common.PTCparam>0) && (common.ptcMatrixType>0))
                ArrayCopy(&res.Mass[npe*npe+mm], &res.Minv[npe*npe+mm], ns, backend);
        }
        else {       
            Int N = nge*npe*ns;
            Int M = nge*npe;
            ShapJac(work, master.shapegt, &tmp.tempg[n2], nge, M, N, backend);                            
            Gauss2Node(handle, &res.Minv[npe*npe*mm], work, master.shapegw, nge, npe, npe*ns, backend);
            if ((common.PTCparam>0) && (common.ptcMatrixType>0))
                ArrayCopy(&res.Mass[npe*npe*mm], &res.Minv[npe*npe*mm], npe*npe*ns, backend);
            Inverse(handle, &res.Minv[npe*npe*mm], work, ipiv, npe, ns, backend);
        }
        mm = mm+ns;
    }                
                
#ifdef HAVE_CUDA            
    if (backend==2) {   
        cudaFree(work);   
        cudaFree(ipiv);    
    }
#endif
    if (backend<2) {
        free(work);
        free(ipiv);
    }    
}

//ApplyMinv(handle, &res.Rq[N], res.Minv, &res.Rq[npe*ncq*e1], app.fc_q, scalar, common.curvedMesh, npe, ncq, e1, e2, backend);              
void ApplyMinv(cublasHandle_t handle, dstype *MinvR, dstype *Minv, dstype *R, dstype *fc_r, 
        dstype scalar, Int curvedmesh, Int npe, Int ncr, Int e1, Int e2, Int backend)
{    
    Int ns = e2-e1;
    Int N = npe*ncr*ns;
    Int M = npe*ncr;
    
    // apply factor fc_r to MinvR
    //ApplyFactor(MinvR, MinvR, fc_r, npe, M, N, backend);    
    ArrayMultiplyScalar(MinvR, scalar, N, backend);
        
    if (curvedmesh==1) { // curved mesh
        //ApplyFactor(R, R, fc_r, npe, M, N, backend);
        ArrayMultiplyScalar(R, scalar, N, backend);
        if (backend==2) {
            ArrayGemmBatch1(&MinvR[0], &Minv[npe*npe*e1], &R[0], npe, ncr, npe, ns, backend);
        //gpuGEMMStridedBatched1(handle, MinvR, &Minv[npe*npe*e1], R, npe, ncr, npe, ns);
        //gpuGEMMBatched1(handle, MinvR, &Minv[npe*npe*e1], R, npe, ncr, npe, ns);     
            // PGEMNMStridedBached(handle, npe, ncr, npe, cublasOne, &Minv[npe*npe*e1], npe,
            //                        R, npe, cublasOne, MinvR, npe, ns, backend);     
        //CHECK_CUBLAS(cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, npe, ncr, npe, 
        //    cublasOne, &Minv[npe*npe*e1], npe, npe*npe, R, npe, npe*ncr, cublasOne, MinvR, npe, npe*ncr, ns));            
//         gemmStridedBatched(cublasHandle_t handle, 
//                       cublasOperation_t transA, cublasOperation_t transB,
//                       int M, int N, int K, 
//                       const T* alpha,
//                       const T* A, int ldA, int strideA, 
//                       const T* B, int ldB, int strideB, 
//                       const T* beta,
//                       T* C, int ldC, int strideC,
//                       int batchCount)            
        }
        else {
            for(int i=0; i<ns; i++)
                Gauss2Node1(handle, &MinvR[npe*ncr*i], &R[npe*ncr*i], &Minv[npe*npe*(i+e1)], npe, npe, ncr, backend); 
        }
    }
    else { // straight mesh
        ApplyFactorJac(R, R, fc_r, &Minv[npe*npe+e1], npe, M, N, backend);                
        //ApplyJacInv(R, R, &Minv[npe*npe+e1], M, N, backend);        
        //ArrayMultiplyScalar(R, scalar, N, backend);
        Gauss2Node1(handle, MinvR, R, &Minv[0], npe, npe, ncr*ns, backend); 
    }         
}

void ApplyMinv(cublasHandle_t handle, dstype *MinvR, dstype *Minv, dstype *R, 
        Int curvedmesh, Int npe, Int ncr, Int e1, Int e2, Int backend)
{    
    Int ns = e2-e1;
    Int N = npe*ncr*ns;
    Int M = npe*ncr;
            
    if (curvedmesh==1) { // curved mesh
        if (backend==2) {
            ArrayGemmBatch(&MinvR[0], &Minv[npe*npe*e1], &R[0], npe, ncr, npe, ns, backend);
        //gpuGEMMStridedBatched1(handle, MinvR, &Minv[npe*npe*e1], R, npe, ncr, npe, ns);
        //gpuGEMMBatched1(handle, MinvR, &Minv[npe*npe*e1], R, npe, ncr, npe, ns);
        }
        else {
            for(int i=0; i<ns; i++)
                Gauss2Node(handle, &MinvR[npe*ncr*i], &R[npe*ncr*i], &Minv[npe*npe*(i+e1)], npe, npe, ncr, backend); 
        }
    }
    else { // straight mesh
        ApplyJacInv(R, R, &Minv[npe*npe+e1], M, N, backend);        
        Gauss2Node(handle, MinvR, R, &Minv[0], npe, npe, ncr*ns, backend); 
    }         
}

#endif  

