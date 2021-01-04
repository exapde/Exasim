#ifndef __RESIDUAL
#define __RESIDUAL

#include "geometry.cpp"
#include "massinv.cpp"
#include "qresidual.cpp"
#include "wresidual.cpp"
#include "uresidual.cpp"
#include "getuhat.cpp"

void GetQ(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nbe1, Int nbe2, Int nbf1, Int nbf2, Int backend)
{    
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncu = common.ncu;// number of compoments of (u)
    Int ncq = common.ncq;// number of compoments of (q)
    Int npe = common.npe; // number of nodes on master element
    Int npf = common.npf; // number of nodes on master element
    Int ne = common.ne2; // number of elements in this subdomain 
    Int N = npe*ncq*ne;
    
    INIT_TIMING;    
    
    START_TIMING;
    // Element integrals
    RqElem(sol, res, app, master, mesh, tmp, common, handle, nbe1, nbe2, backend);
    END_TIMING(15);    
        
    START_TIMING;
    // Face integrals
    RqFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
    END_TIMING(16);       
        
    // elements in the range [e1, e2)
    Int e1 = common.eblks[3*nbe1]-1;    
    Int e2 = common.eblks[3*(nbe2-1)+1];
       
    //cout<<res.Rh[common.nf*common.npf*ncq-1]<<endl;
    //print3darray(res.Rh, common.npf, ncq, 5);
    // assemble face residual vector into element residual vector
    PutFaceNodes(res.Rq, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, npf, npe, ncq, e1, e2, 0, backend);
    
    if (common.wave==1)
        // get the source term due to the time derivative (for wave problem)  
        ArrayExtract(&res.Rq[N], sol.sdg, npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2, backend);  
    else
        // set it to zero
        ArraySetValue(&res.Rq[N], zero, npe*ncq*(e2-e1), backend);
        
    dstype scalar = one;
    if (common.wave==1)
        scalar = one/common.dtfactor;

    START_TIMING;
    // Apply the mass matrix inverse and the factor fc_q to the residual
    ApplyMinv(handle, &res.Rq[N], res.Minv, &res.Rq[npe*ncq*e1], app.fc_q, scalar, common.curvedMesh, npe, ncq, e1, e2, backend);              
    END_TIMING(17);   
                
    START_TIMING;
    // append q into udg
    ArrayInsert(sol.udg, &res.Rq[N], npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2, backend);           
    END_TIMING(18);       
}

void GetW(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nbe1, Int nbe2, Int nbf1, Int nbf2, Int backend)
{        
    Int nc = common.nc; // number of compoments of (u, q, p)
    Int ncw = common.ncw;// number of compoments of (w)
    Int nco = common.nco;// number of compoments of (o)
    Int ncx = common.ncx;// number of compoments of (xdg)        
    Int npe = common.npe; // number of nodes on master element    
    for (Int j=nbe1; j<nbe2; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];
                
        if (common.wave==1) {
            // dw/dt = u
            dstype scalar = one/common.dtfactor;
            ArrayExtract(tmp.tempn, sol.udg, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, e1, e2, backend);                                                  
            ArrayAXPBY(&sol.wdg[npe*ncw*e1], tmp.tempn, &sol.wsrc[npe*ncw*e1], scalar, scalar, npe*ncw*(e2-e1), backend);                        
        }
        else {
            // calculate the source term Sourcew(xdg, udg, odg, wdg)
            SourcewDriver(tmp.tempn, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1], 
                    &sol.wdg[npe*ncw*e1], mesh, master, app, sol, tmp, common, npe, e1, e2, backend);            
            if (common.dae_steps==0) {
                // 1.0/(alpha*dirkd/dt + beta)
                dstype scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta);
                
                //cout<<common.dae_alpha<<" "<<common.dtfactor<<" "<<scalar<<" "<<app.physicsparam[0]<<endl;

                // calculate w = (1/(alpha*dirkd/dt + beta))*(alpha*wsrc + Sourcew(xdg, udg, odg, wdg))  
                ArrayAXPBY(&sol.wdg[npe*ncw*e1], &sol.wsrc[npe*ncw*e1], tmp.tempn, common.dae_alpha*scalar, scalar, npe*ncw*(e2-e1), backend);                    
            }
            else {
                // calculate tmp = alpha*wsrc + Sourcew(xdg, udg, odg, wdg)) 
                ArrayAXPBY(tmp.tempn, &sol.wsrc[npe*ncw*e1], tmp.tempn, common.dae_alpha, one, npe*ncw*(e2-e1), backend);                    

                // 1.0/(alpha*dirkd/dt + beta)
                dstype scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta + common.dae_gamma);
                
                //cout<<scalar<<endl;
                
                // calculate w = (1/(alpha*dirkd/dt + beta + gamma))*(gamma*walg + alpha*wsrc + Sourcew(xdg, udg, odg, wdg))  
                ArrayAXPBY(&sol.wdg[npe*ncw*e1], &sol.wdual[npe*ncw*e1], tmp.tempn, common.dae_gamma*scalar, scalar, npe*ncw*(e2-e1), backend);                                    
            }                
        }                    
    }    
}

void RuResidual(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
   Int nbe1u, Int nbe2u, Int nbe1q, Int nbe2q, Int nbf1, Int nbf2, Int backend)
{    
    // compute uhat
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
    
    // compute q
    if (common.ncq>0)
        GetQ(sol, res, app, master, mesh, tmp, common, handle, nbe1q, nbe2q, nbf1, nbf2, backend);                
    
//     // compute w
    if (common.ncw>0)
         GetW(sol, res, app, master, mesh, tmp, common, handle, nbe1q, nbe2q, nbf1, nbf2, backend);                
        
    // Element integrals
    RuElem(sol, res, app, master, mesh, tmp, common, handle, nbe1u, nbe2u, backend);    
                    
    // Face integrals
    RuFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
        
    Int e1 = common.eblks[3*nbe1u]-1;    
    Int e2 = common.eblks[3*(nbe2u-1)+1];    
    
    // assemble face residual vector into element residual vector
    PutFaceNodes(res.Ru, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0, backend);
        
    // change sign 
    //ArrayMultiplyScalar(res.Ru, minusone, common.ndof1, backend);        
}

#ifdef  HAVE_MPI

void GetQMPI(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{  
    // non-blocking send and receive solutions on exterior and outer elements to neighbors
    
    Int bsz = common.npe*common.ncu;
    Int nudg = common.npe*common.nc;
    Int n;
    
    /* copy some portion of u to buffsend */
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend, backend);
    //for (n=0; n<common.nelemsend; n++)         
    //    ArrayCopy(&tmp.buffsend[bsz*n], &sol.udg[nudg*common.elemsend[n]], bsz, backend);            
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
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
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
            
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
       
    // calculate q for interior elements
    if (common.ncq>0)         
        GetQ(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, 0, common.nbf, backend);        
    
    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);

    /* copy buffrecv to udg */
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv, backend);
    //for (n=0; n<common.nelemrecv; n++) 
    //    ArrayCopy(&sol.udg[nudg*common.elemrecv[n]], &tmp.buffrecv[bsz*n], bsz, backend);        
    
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
    
    // calculate q for interface and exterior elements
    if (common.ncq>0)         
        GetQ(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe2, 0, common.nbf, backend);                
}

void RuResidualMPI(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{      
    // non-blocking send and receive solutions on exterior and outer elements to neighbors    
    Int bsz = common.npe*common.ncu;
    Int nudg = common.npe*common.nc;
    Int n;
   
    INIT_TIMING;        
    
    START_TIMING;
    /* copy some portion of u to buffsend */
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend, backend);
    //for (n=0; n<common.nelemsend; n++)         
    //    ArrayCopy(&tmp.buffsend[bsz*n], &sol.udg[nudg*common.elemsend[n]], bsz, backend);                    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

    END_TIMING(13);    
    
    START_TIMING;    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
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
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
    END_TIMING(6);    
                    
    START_TIMING; 
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);            
    END_TIMING(7);    
    
    START_TIMING;  
    // calculate q for interior elements
    if (common.ncq>0)         
        GetQ(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, 0, common.nbf, backend);        
    if (common.ncw>0)         
        GetW(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, 0, common.nbf, backend);            
    END_TIMING(8);    
    
    START_TIMING; 
    // calculate Ru for interior elements
    RuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, backend);    
    END_TIMING(9);    
        
    START_TIMING; 
    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);
    END_TIMING(10);    

    START_TIMING; 
    /* copy buffrecv to udg */
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv, backend);
 
    END_TIMING(14);    
    
    START_TIMING; 
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
    END_TIMING(7);    

    START_TIMING; 
    // calculate q for interface and exterior elements
    if (common.ncq>0)         
        GetQ(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe2, 0, common.nbf, backend);        
    if (common.ncw>0)         
        GetW(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe2, 0, common.nbf, backend);        
    
    // calculate Ru for interface elements
    RuElem(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe1, backend);    
    END_TIMING(11);    
     
    START_TIMING; 
    // calculate Ru for all faces
    RuFace(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
    END_TIMING(12);    
        
    // assemble face residual vector into element residual vector
    Int e1 = common.eblks[3*0]-1;    
    Int e2 = common.eblks[3*(common.nbe1-1)+1];       
    PutFaceNodes(res.Ru, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0, backend);
    
    // change sign 
    //ArrayMultiplyScalar(res.Ru, minusone, common.ndof1, backend);      
}

void RuResidualMPI1(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{      
    // non-blocking send and receive solutions on exterior and outer elements to neighbors    
    Int bsz = common.npe*common.ncu;
    Int nudg = common.npe*common.nc;
    Int n;
    
    /* copy some portion of u to buffsend */
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend, backend);
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
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
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
                    
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);                    
    
    // calculate q for interior elements
    if (common.ncq>0)         
        GetQ(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, 0, common.nbf, backend);        
    
    // calculate Ru for interior elements
    RuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, backend);    
        
    // calculate Ru for interior faces
    RuFace(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf0, backend);
    
    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);

    /* copy buffrecv to udg */
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv, backend);
    
    // compute uhat for all faces
    GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);

    // calculate q for interface and exterior elements
    if (common.ncq>0)
        GetQ(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe2, common.nbf0, common.nbf, backend);        
                
    // calculate Ru for interface elements
    RuElem(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe1, backend);    
     
    // calculate Ru for all other faces
    RuFace(sol, res, app, master, mesh, tmp, common, handle, common.nbf0, common.nbf, backend);
        
    // change sign 
    //ArrayMultiplyScalar(res.Ru, minusone, common.ndof1, backend);      
}

#endif

void Residual(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{    
    if (common.mpiProcs>1) { // mpi processes
#ifdef  HAVE_MPI        
        RuResidualMPI(sol, res, app, master, mesh, tmp, common, handle, backend);
#endif        
    }
    else {
        RuResidual(sol, res, app, master, mesh, tmp, common, handle,
                0, common.nbe, 0, common.nbe, 0, common.nbf, backend);
    }        
            
    // change sign for matrix-vector product
    ArrayMultiplyScalar(res.Ru, minusone, common.ndof1, backend);      
        
    //common.dtfactor
    if (common.tdep==1) 
        ArrayMultiplyScalar(res.Ru, one/common.dtfactor, common.ndof1, backend);                
    
    if (common.debugMode==1) {
        writearray2file(common.fileout + "_uh.bin", sol.uh, common.npf*common.ncu*common.nf, backend);
        writearray2file(common.fileout + "_udg.bin", sol.udg, common.npe*common.nc*common.ne, backend);
        writearray2file(common.fileout + "_Ru.bin", res.Ru, common.npe*common.ncu*common.ne, backend);
    }    
}

void ComputeQ(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{
    if (common.mpiProcs>1) {
#ifdef  HAVE_MPI        
        GetQMPI(sol, res, app, master, mesh, tmp, common, handle, backend);
#endif                
    }
    else {
        // compute uhat
        GetUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
        GetQ(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe, 0, common.nbf, backend);    
    }        
    
    if (common.debugMode==1) {
        writearray2file(common.fileout + "_udg.bin", sol.udg, common.npe*common.nc*common.ne, backend);
    }
}

#endif

