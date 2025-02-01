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

void hdgMatVec(dstype *w, dstype *AE, dstype *v, dstype *ve, dstype *we, resstruct &res, appstruct &app, 
        meshstruct &mesh, commonstruct &common, tempstruct &tmp, cublasHandle_t handle, Int backend)
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
//     for (int n=0; n<common.nelemsend; n++)  {       
//       ArrayCopy(&tmp.buffsend[bsz*n], &we[bsz*common.elemsend[n]], bsz);     
//     }
    GetCollumnAtIndex(tmp.buffsend, we, mesh.elemsend, bsz, common.nelemsend);   

#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
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

    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      int neI = common.nintfaces;
      
      // perform  matrix-vector products for interface elements
      // (ncu12 * npf) * (ncu * npf * nfe ) * neI x (ncu * npf * nfe ) * neI  = (ncu12 * npf) * neI    
      PGEMNMStridedBached(handle, szRi, 1, m, one, res.Hi, szRi, &ve[m*ne0], m, zero, res.Ri, szRi, neI, backend); 
      
//       for (int n=0; n<common.nfacesend; n++)  {       
//         ArrayCopy(&tmp.bufffacesend[szRi*n], &res.Ri[szRi*common.facesend[n]], szRi);     
//       }
      GetCollumnAtIndex(tmp.bufffacesend, res.Ri, mesh.facesend, szRi, common.nfacesend);   
       
#ifdef HAVE_HIP
      hipDeviceSynchronize();
#endif
      
      /* non-blocking send */
      psend = 0; 
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nsend = common.facesendpts[n]*szRi;
          if (nsend>0) {
              MPI_Isend(&tmp.bufffacesend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              psend += nsend;
              request_counter += 1;
          }
      }
      
      /* non-blocking receive */
      precv = 0;
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nrecv = common.facerecvpts[n]*szRi;
          if (nrecv>0) {
              MPI_Irecv(&tmp.bufffacerecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              precv += nrecv;
              request_counter += 1;
          }
      }      
    }
    
    // perform matrix-vector products for interior elements
    PGEMNMStridedBached(handle, m, 1, m, one, AE, m, ve, m, zero, we, m, ne0, backend); 

    // copy buffrecv to we
    MPI_Waitall(request_counter, common.requests, common.statuses);    
//     for (int n=0; n<common.nelemrecv; n++) {        
//       ArrayCopy(&we[bsz*common.elemrecv[n]], &tmp.buffrecv[bsz*n], bsz);       
//     }
    PutCollumnAtIndex(we, tmp.buffrecv, mesh.elemrecv, bsz, common.nelemrecv);   

    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      
//       for (int n=0; n<common.nfacerecv; n++) {        
//         ArrayCopy(&res.Ri[szRi*common.facerecv[n]], &tmp.bufffacerecv[szRi*n], szRi); 
//       }      
      PutCollumnAtIndex(res.Ri, tmp.bufffacerecv, mesh.facerecv, szRi, common.nfacerecv);   
      
      PutBoudaryNodes(we, res.Ri, mesh.intfaces, mesh.faceperm, app.interfacefluxmap, nfe, npf, ncu12, ncu, common.nintfaces); 
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

    //printf("hdgAssembleLinearSystemMPI: %d %d %d %d\n", common.nbe0, common.nbe1, common.nelemsend, common.nnbsd);
    
    // perform HDG descrization for interface elements
    for (Int j=common.nbe0; j<common.nbe1; j++) {     // fixed bug here             
      uEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      uEquationSchurBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
    }                             
        
    // copy H and Rh to buffsend
//     for (int n=0; n<common.nelemsend; n++)  {       
//       ArrayCopy(&tmp.buffsend[bsz*n], &res.H[szH*common.elemsend[n]], szH);     
//       ArrayCopy(&tmp.buffsend[bsz*n + szH], &res.Rh[szR*common.elemsend[n]], szR);         
//     }
    GetCollumnAtIndex(tmp.buffsend, res.H, mesh.elemsend, 0, bsz, szH, common.nelemsend);
    GetCollumnAtIndex(tmp.buffsend, res.Rh, mesh.elemsend, szH, bsz, szR, common.nelemsend);
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
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

    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      
//       if (common.mpiRank==1)
//         print2darray(res.Ri, szRi, common.nfacesend);
      
//       for (int n=0; n<common.nfacesend; n++)  {       
//         ArrayCopy(&tmp.bufffacesend[szRi*n], &res.Ri[szRi*common.facesend[n]], szRi);     
//       }
      GetCollumnAtIndex(tmp.bufffacesend, res.Ri, mesh.facesend, szRi, common.nfacesend);   

      //printf("hdgAssembleLinearSystemMPI: %d %d %d %d %d %d %d %d %d %d %d\n", common.mpiRank, common.nnbintf, common.nbintf[0], common.facesendpts[0], common.facerecvpts[0], common.nfacesend, common.nelemsend, szRi, common.coupledcondition, common.coupledboundarycondition, app.interfacefluxmap[0]);
      
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
      
      /* non-blocking send */
      psend = 0; 
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nsend = common.facesendpts[n]*szRi;
          if (nsend>0) {
              MPI_Isend(&tmp.bufffacesend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              psend += nsend;
              request_counter += 1;
          }
      }
      
      /* non-blocking receive */
      precv = 0;
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nrecv = common.facerecvpts[n]*szRi;
          if (nrecv>0) {
              MPI_Irecv(&tmp.bufffacerecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              precv += nrecv;
              request_counter += 1;
          }
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
//     for (int n=0; n<common.nelemrecv; n++) {        
//       ArrayCopy(&res.H[szH*common.elemrecv[n]], &tmp.buffrecv[bsz*n], szH); 
//       ArrayCopy(&res.Rh[szR*common.elemrecv[n]], &tmp.buffrecv[bsz*n + szH], szR);            
//     }
    PutCollumnAtIndex(res.H, tmp.buffrecv, mesh.elemrecv, 0, bsz, szH, common.nelemrecv);
    PutCollumnAtIndex(res.Rh, tmp.buffrecv, mesh.elemrecv, szH, bsz, szR, common.nelemrecv);
    
    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      
//       for (int n=0; n<common.nfacerecv; n++) {        
//         ArrayCopy(&res.Ri[szRi*common.facerecv[n]], &tmp.bufffacerecv[szRi*n], szRi); 
//       }      
      PutCollumnAtIndex(res.Ri, tmp.bufffacerecv, mesh.facerecv, szRi, common.nfacerecv);  
      
      PutBoudaryNodes(res.Rh, res.Ri, mesh.intfaces, mesh.faceperm, app.interfacefluxmap, nfe, npf, ncu12, ncu, common.nintfaces); 
      
//       if (common.mpiRank==1)
//         print2darray(res.Ri, szRi, common.nfacerecv);
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
//     for (int n=0; n<common.nelemsend; n++)  {       
//       ArrayCopy(&tmp.buffsend[bsz*n], &res.Rh[szR*common.elemsend[n]], szR);         
//     }
    GetCollumnAtIndex(tmp.buffsend, res.Rh, mesh.elemsend, szR, common.nelemsend);   
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
        
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

    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      
//       for (int n=0; n<common.nfacesend; n++)  {       
//         ArrayCopy(&tmp.bufffacesend[szRi*n], &res.Ri[szRi*common.facesend[n]], szRi);     
//       }            
      GetCollumnAtIndex(tmp.bufffacesend, res.Ri, mesh.facesend, szRi, common.nfacesend);   
      
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
      /* non-blocking send */
      psend = 0; 
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nsend = common.facesendpts[n]*szRi;
          if (nsend>0) {
              MPI_Isend(&tmp.bufffacesend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              psend += nsend;
              request_counter += 1;
          }
      }
      
      /* non-blocking receive */
      precv = 0;
      for (int n=0; n<common.nnbintf; n++) {
          neighbor = common.nbintf[n];
          nrecv = common.facerecvpts[n]*szRi;
          if (nrecv>0) {
              MPI_Irecv(&tmp.bufffacerecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                    MPI_COMM_WORLD, &common.requests[request_counter]);
              precv += nrecv;
              request_counter += 1;
          }
      }      
    }
    
    // perform HDG descrization for interior elements
    for (Int j=0; j<common.nbe0; j++) {                  
      RuEquationElemBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);
      RuEquationElemFaceBlock(sol, res, app, master, mesh, tmp, common, handle, j, backend);        
    }                        
    
    // copy buffrecv to Rh 
    MPI_Waitall(request_counter, common.requests, common.statuses);    
//     for (int n=0; n<common.nelemrecv; n++) {        
//       ArrayCopy(&res.Rh[szR*common.elemrecv[n]], &tmp.buffrecv[bsz*n], szR);            
//     }
    PutCollumnAtIndex(res.Rh, tmp.buffrecv, mesh.elemrecv, szR, common.nelemrecv);         
    
    if ((common.nnbintf > 0) && (common.nfacesend > 0) && (common.coupledcondition>0)) {
      int ncu12 = common.szinterfacefluxmap;
      int szRi = ncu12*npf;  
      
//       for (int n=0; n<common.nfacerecv; n++) {        
//         ArrayCopy(&res.Ri[szRi*common.facerecv[n]], &tmp.bufffacerecv[szRi*n], szRi); 
//       }      
      PutCollumnAtIndex(res.Ri, tmp.bufffacerecv, mesh.facerecv, szRi, common.nfacerecv);  

      PutBoudaryNodes(res.Rh, res.Ri, mesh.intfaces, mesh.faceperm, app.interfacefluxmap, nfe, npf, ncu12, ncu, common.nintfaces);       
    }
        
    // assemble RHS vector b from res.Rh using the FIRST elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, npf, nfe, ncu, common.nf);

    // assemble RHS vector b from res.Rh using the SECOND elements in mesh.f2e
    PutElementFaceNodes(b, res.Rh, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    
    
}

#endif

#endif


