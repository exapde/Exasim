#ifndef __MATVEC
#define __MATVEC

#include "ioutilities.cpp"
void MatVec(dstype *w, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master,
      meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, dstype *v, 
      dstype *u, dstype *Ru, Int backend)
{   
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    //Int nd = common.nd;
    Int N = npe*ncu*ne;
    
    Int order = common.matvecOrder;
    dstype epsilon = common.matvecTol;
#ifdef HAVE_ENZYME
//TODO: there might be a cleaner way to do this...matvecAD v. matvecFD functions? 
    ArrayInsert(sol.dudg, v, npe, nc, ne, 0, npe, 0, ncu, 0, ne); //insert v into dudgs
    dResidual(sol, res, app, master, mesh, tmp, common, handle, backend);
    ArrayAXPBY(w, u, res.dRu, 0.0, 1, N);
#else
    if (order==1) {
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N);
                            
        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u))/epsilon    
        ArrayAXPBY(w, res.Ru, Ru, 1.0/epsilon, -1.0/epsilon, N);
    }
    else if (order==2) {
        // calculate w = u - epsilon*v
        ArrayAXPBY(w, u, v, 1.0, -epsilon, N);

        // insert (u-epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne);  

        // compute the residual R(u-epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);

        // copy res.Ru to Ru
        ArrayCopy(Ru, res.Ru, N);
        
        // calculate w = u + epsilon*v
        ArrayAXPBY(w, u, v, 1.0, epsilon, N);

        // insert (u+epsilon*v) into udg
        ArrayInsert(sol.udg, w, npe, nc, ne, 0, npe, 0, ncu, 0, ne);  

        // compute the residual R(u+epsilon*v)
        Residual(sol, res, app, master, mesh, tmp, common, handle, backend);
        
        // calculate w = J(u)*v = (R(u+epsilon*v)-R(u-epsilon*v))/(2*epsilon)    
        ArrayAXPBY(w, res.Ru, Ru, 0.5/epsilon, -0.5/epsilon, N);
    }
    else
        error("Matrix-vector multiplication order is not implemented");
#endif
}

void hdgAssembleRHS(dstype *R, dstype *Rh, meshstruct &mesh, commonstruct &common)
{   
    Int nf = common.nf; // number of faces in this subdomain
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int nfe = common.nfe; // number of faces in each element

    // send Rh to neighbor elements 

     // ncu * npf * nfe * ne -> ncu * npf * nf
    PutElementFaceNodes(R, Rh, mesh.f2e, npf, nfe, ncu, nf);
    PutElementFaceNodes(R, Rh, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    // receive Rh from neighbor elements 

    if (common.debugMode == 1) {
      writearray2file(common.fileout + "hdgAssembleRHS.bin", R, ncu * npf * nf, common.backend);
    }
}

void hdgBlockJacobi(dstype *BE, dstype *AE, dstype *BEtmp, int *ipiv, meshstruct &mesh, commonstruct &common, cublasHandle_t handle, Int backend)
{   
    Int nf = common.nf; // number of faces in this subdomain
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int nfe = common.nfe; // number of faces in each element
    Int ncf = ncu*npf;

     // ncf * nfe * ncf * nfe * ne -> ncf * ncf * nf
    blockJacobi(BE, AE, mesh.f2e, npf, nfe, ncu, nf);
    blockJacobi(BE, AE, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    
    Inverse(handle, BE, BEtmp, ipiv, ncf, nf, backend); 

    if (common.debugMode == 1) {
      writearray2file(common.fileout + "hdgBlockJacobi.bin", BE, ncf * ncf * nf, backend);
    }
}

void hdgApplyBlockJacobi(dstype *w, dstype *BE, dstype *v, commonstruct &common, cublasHandle_t handle, Int backend)
{   
    Int nf = common.nf; // number of faces in this subdomain
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int m = ncu*npf;  

    // (ncf)  * (ncf) * nf x (ncf) * nf -> (ncf) * nf
    PGEMNMStridedBached(handle, m, 1, m, one, BE, m, v, m, zero, w, m, nf, backend); 

    if (common.debugMode == 1) {
      writearray2file(common.fileout + "hdgApplyBlockJacobi.bin", w, m * nf, backend);
    }
}

void hdgGetDUDG(dstype *w, dstype *F, dstype *duh, dstype *ve, meshstruct &mesh, 
        commonstruct &common,  Int backend)
{   
    Int ne = common.ne1; // number of elements in this subdomain 
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int nfe = common.nfe; // number of faces in each element
    Int m = ncu*npf*nfe;  
    Int n = common.npe*ncu;

    // ncu * npf * nf -> ncu * npf * nfe * ne
    GetElementFaceNodes(ve, duh, mesh.elemcon, npf*nfe, ncu, 0, ne, 2);

    // (npe * ncu)  * (ncu * npf * nfe) * ne x (ncu * npf * nfe) * ne -> (npe * ncu) * ne
    PGEMNMStridedBached(common.cublasHandle, n, 1, m, one, F, n, ve, m, one, w, n, ne, backend); 

    if (common.debugMode == 1) {
      writearray2file(common.fileout + "hdgGetDUDG.bin", w, n * ne, backend);
    }
}

void hdgMatVec(dstype *w, dstype *AE, dstype *v, dstype *ve, dstype *we, meshstruct &mesh, 
        commonstruct &common, tempstruct &tmp, cublasHandle_t handle, Int backend)
{   
    Int ne1 = common.ne1; // number of elements in this subdomain 
    Int nf = common.nf; // number of faces in this subdomain
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int nfe = common.nfe; // number of faces in each element
    Int m = ncu*npf*nfe;  

    // ncu * npf * nf -> ncu * npf * nfe * ne
    GetElementFaceNodes(ve, v, mesh.elemcon, npf*nfe, ncu, 0, ne1, 2);

#ifdef HAVE_MPI     
    Int ne0 = common.ne0; // number of interior elements in this subdomain
    Int bsz = ncu*npf*nfe;  

    // perform matrix-vector products for interface elements
    PGEMNMStridedBached(handle, m, 1, m, one, &AE[m*m*ne0], m, &ve[m*ne0], m, zero, &we[m*ne0], m, ne1-ne0, backend); 
    
    // copy we to buffsend
    for (int n=0; n<common.nelemsend; n++)  {       
      ArrayCopy(&tmp.buffsend[bsz*n], &we[bsz*common.elemsend[n]], bsz);     
    }

    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    Int nrecv, precv = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    // perform matrix-vector products for interior elements
    PGEMNMStridedBached(handle, m, 1, m, one, AE, m, ve, m, zero, we, m, ne0, backend); 

    // copy buffrecv to we
    MPI_Waitall(request_counter, common.requests, common.statuses);    
    for (int n=0; n<common.nelemrecv; n++) {        
      ArrayCopy(&we[bsz*common.elemrecv[n]], &tmp.buffrecv[bsz*n], bsz);       
    }

    // assemble vector w from we using the FIRST elements in mesh.f2e
    PutElementFaceNodes(w, we, mesh.f2e, npf, nfe, ncu, nf);

    // assemble vector w from we using the SECOND elements in mesh.f2e
    PutElementFaceNodes(w, we, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);        
#else 
    // (ncu * npf * nfe)  * (ncu * npf * nfe) * ne x (ncu * npf * nfe) * ne -> (ncu * npf * nfe) * ne
    PGEMNMStridedBached(handle, m, 1, m, one, AE, m, ve, m, zero, we, m, ne1, backend); 
    
     // ncu * npf * nfe * ne -> ncu * npf * nf
    PutElementFaceNodes(w, we, mesh.f2e, npf, nfe, ncu, nf);
    PutElementFaceNodes(w, we, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    if (common.debugMode == 1) {
      writearray2file(common.fileout + "hdgMatVec.bin", w, ncu * npf * nf, backend);
    }
#endif
}

void hdgAssembleSerial(dstype *b, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,  cublasHandle_t handle, Int backend)
{
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int npe = common.npe; // number of nodes on master element
    Int nfe = common.nfe; // number of faces in each element
    Int ncq = common.ncq; // number of compoments of (q)
    Int ncf = ncu*npf;

    string filename;
    if (common.mpiProcs==1)
      filename = common.fileout;
    else
      filename = common.fileout + NumberToString(common.mpiRank);      

    // perform HDG descrization for interface elements
    for (Int j=0; j<common.nbe1; j++) {                  
      uEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);      
      uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);

      Int e1 = common.eblks[3*j]-1;
      Int e2 = common.eblks[3*j+1];            
      Int ne = e2-e1;

      if (common.mpiRank==0)  {        
        // print2darray(&res.Ru[npe*e1], npe, ne);   
        // print3darray(&res.D[0], npe, npe, ne);   
        // print3darray(&res.B[0], npe, npe, ne);   
        // print3darray(&res.B[npe*npe*ne], npe, npe, ne);   
        print3darray(&res.H[npf*nfe*npf*nfe*e1], npf*nfe, npf*nfe, ne);   
      }
    }                             

#ifdef  HAVE_MPI     
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    error("stop here");

    Int ne = common.ne;
    Int e1 = 0;
    Int n = npe*ncu; 
    Int m = npf*nfe*ncu;

    writearray2file(filename + "uEquationElem_Ru.bin", res.Ru, npe*ncu*ne, backend);
    writearray2file(filename + "uEquationElem_D.bin", res.D, npe*npe*ncu*ncu*ne, backend);  
    if (ncq > 0) writearray2file(filename + "uEquationElem_B.bin", res.B, npe*npe*ncu*ncq*ne, backend);      

    for (Int j=0; j<common.nbe; j++) {                  
      uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                             

    writearray2file(filename + "uEquationElemFace_Ru.bin", &res.Ru[npe*ncu*e1], npe*ne*ncu, backend);
    writearray2file(filename + "uEquationElemFace_D.bin", res.D, npe*npe*ncu*ncu*ne, backend);
    writearray2file(filename + "uEquationElemFace_B.bin", res.B, npe*npe*ncu*ncq*ne, backend);
    writearray2file(filename + "uEquationElemFace_F.bin", &res.F[npe*npf*nfe*ncu*ncu*e1], npe*npf*nfe*ncu*ncu*ne, backend);
    writearray2file(filename + "uEquationElemFace_Rh.bin", &res.Rh[npf*nfe*ncu*e1], npf*nfe*ne*ncu, backend);
    writearray2file(filename + "uEquationElemFace_K.bin", res.K, npf*nfe*npe*ne*ncu*ncu, backend);  
    writearray2file(filename + "uEquationElemFace_G.bin", res.G, npf*nfe*npe*ne*ncu*ncq, backend);  
    writearray2file(filename + "uEquationElemFace_H.bin", &res.H[npf*nfe*npf*nfe*ncu*ncu*e1], npf*nfe*npf*nfe*ne*ncu*ncu, backend);      

    for (Int j=0; j<common.nbe; j++) {                  
      uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                             

    dstype *DinvF = &res.F[n*m*e1];    
    dstype *Ru = &res.Ru[n*e1];
    dstype *DinvH =&res.H[m*m*e1];
    dstype *Rh = &res.Rh[m*e1];
    writearray2file(filename + "uEquationElemSchur_DinvF.bin", DinvF, n*m*ne, backend);
    writearray2file(filename + "uEquationElemSchur_Ru.bin", Ru, n*ne, backend);
    writearray2file(filename + "uEquationElemSchur_DinvH.bin", DinvH, m*m*ne, backend);  
    writearray2file(filename + "uEquationElemSchur_Rh.bin", Rh, m*ne, backend);      

    // assemble RHS vector b from res.Rh using the FIRST elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, npf, nfe, ncu, common.nf);    
    // assemble RHS vector b from res.Rh using the SECOND elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    // assemble block Jacobi matrix from res.H using the FIRST elements in mesh.f2e
    blockJacobi(res.K, res.H, mesh.f2e, npf, nfe, ncu, common.nf);
    // assemble block Jacobi matrix from res.H using the SECOND elements in mesh.f2e
    blockJacobi(res.K, res.H, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    // inverse block Jacobi matrix
    Inverse(handle, res.K, tmp.tempg, res.ipiv, ncf, common.nf, backend); 


#ifdef  HAVE_MPI     
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    error("stop here");
}

#ifdef  HAVE_MPI     
void hdgAssembleLinearSystemMPI(dstype *b, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,  cublasHandle_t handle, Int backend)
{
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int npe = common.npe; // number of nodes on master element
    Int nfe = common.nfe; // number of faces in each element
    Int ncf = ncu*npf;
    Int szR = ncu*npf*nfe;  
    Int szH = szR*szR;
    Int bsz = szH + szR; // send and receive H and Rh together   

    // perform HDG descrization for interface elements
    for (Int j=common.nbe0; j<common.nbe1; j++) {     // fixed bug here             
      uEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                             
        
    // copy H and Rh to buffsend
    for (int n=0; n<common.nelemsend; n++)  {       
      ArrayCopy(&tmp.buffsend[bsz*n], &res.H[szH*common.elemsend[n]], szH);     
      ArrayCopy(&tmp.buffsend[bsz*n + szH], &res.Rh[szR*common.elemsend[n]], szR);         
    }

    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    Int nrecv, precv = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    // perform HDG descrization for interior elements
    for (Int j=0; j<common.nbe0; j++) {                  
      uEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                        

    // copy buffrecv to H and Rh 
    MPI_Waitall(request_counter, common.requests, common.statuses);    
    for (int n=0; n<common.nelemrecv; n++) {        
      ArrayCopy(&res.H[szH*common.elemrecv[n]], &tmp.buffrecv[bsz*n], szH); 
      ArrayCopy(&res.Rh[szR*common.elemrecv[n]], &tmp.buffrecv[bsz*n + szH], szR);            
    }

    ArraySetValue(b, 0.0, ncu*npf*common.nf); // fix bug here

    // assemble RHS vector b from res.Rh using the FIRST elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, npf, nfe, ncu, common.nf);

    // assemble RHS vector b from res.Rh using the SECOND elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    // assemble block Jacobi matrix from res.H using the FIRST elements in mesh.f2e
    blockJacobi(res.K, res.H, mesh.f2e, npf, nfe, ncu, common.nf);

    // assemble block Jacobi matrix from res.H using the SECOND elements in mesh.f2e
    blockJacobi(res.K, res.H, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    

    // inverse block Jacobi matrix
    Inverse(handle, res.K, tmp.tempg, res.ipiv, ncf, common.nf, backend); 

    // string filename = common.fileout + NumberToString(common.mpiRank);    
    // writearray2file(filename + "hdgAssembleRHS.bin", b, ncu * npf * common.nf, common.backend);    

//     if ((common.mpiRank==0) || (common.mpiRank==0))  {
//       int m = npf*nfe*ncu; 
//       int n = npe*ncu;
//       // print2darray(&res.Ru[0], n, common.ne1);       
//       print2darray(&res.Rh[0], m, common.ne);         
//       // print3darray(&res.F[0], n, m, common.ne1);    
//       // print3darray(&res.H[0], m, m, common.ne1);     
//       cout<<common.nf<<endl; 
//       print2darray(b, ncu*npf, common.nf);             
//     }
  
//     for (Int j=0; j<common.mpiProcs; j++) {       
//       if (common.mpiRank==j) {
//         cout<<common.mpiRank<<endl;                       
//         print2darray(b, ncu*npf, common.nf);             
// // #ifdef  HAVE_MPI     
// //         MPI_Barrier(MPI_COMM_WORLD);
// #endif
//       }
//     }

    // cout<<common.mpiRank<<endl;                       
//     print2darray(b, ncu*npf, common.nf);             
//     // cout<<PNORM(common.cublasHandle, ncu*npf*nfe*common.nf, b, backend)<<endl;

// #ifdef  HAVE_MPI     
//     MPI_Barrier(MPI_COMM_WORLD);
// #endif

//     error("stop here");

}

void hdgAssembleResidualMPI(dstype *b, solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,  cublasHandle_t handle, Int backend)
{
    Int ncu = common.ncu;// number of compoments of (u)
    Int npf = common.npf; // number of nodes on master face           
    Int npe = common.npe; // number of nodes on master element
    Int nfe = common.nfe; // number of faces in each element
    Int ncf = ncu*npf;
    Int szR = ncu*npf*nfe;  
    Int bsz = szR; // send and receive H and Rh together   

    // perform HDG descrization for interface elements
    for (Int j=common.nbe0; j<common.nbe1; j++) {     // fixed bug here             
      RuEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      RuEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);        
    }                             
        
    // copy H and Rh to buffsend
    for (int n=0; n<common.nelemsend; n++)  {       
      ArrayCopy(&tmp.buffsend[bsz*n], &res.Rh[szR*common.elemsend[n]], szR);         
    }

    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    Int nrecv, precv = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    // perform HDG descrization for interior elements
    for (Int j=0; j<common.nbe0; j++) {                  
      RuEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      RuEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);        
    }                        

    // copy buffrecv to Rh 
    MPI_Waitall(request_counter, common.requests, common.statuses);    
    for (int n=0; n<common.nelemrecv; n++) {        
      ArrayCopy(&res.Rh[szR*common.elemrecv[n]], &tmp.buffrecv[bsz*n], szR);            
    }

    // assemble RHS vector b from res.Rh using the FIRST elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, npf, nfe, ncu, common.nf);

    // assemble RHS vector b from res.Rh using the SECOND elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    
}
#endif

#endif


