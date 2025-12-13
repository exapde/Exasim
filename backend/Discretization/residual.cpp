/*
    residual.cpp

    This file contains the implementation of residual and related computations for the Exasim backend Discretization module.
    The main functionalities provided include:

    - DG2CGAVField: Converts a discontinuous Galerkin (DG) artificial viscosity field to a continuous Galerkin (CG) field, with smoothing iterations.
    - GetQ: Computes the q variable (e.g., gradients or fluxes) for the solution, including element and face integrals, and applies mass matrix inversion.
    - GetW: Computes the w variable (e.g., auxiliary or thermodynamic variables) using either explicit time-stepping or Newton's method for nonlinear systems.
    - GetAv: Computes and smooths the artificial viscosity field, and interpolates it to Gauss points for elements and faces.
    - RuResidual: Assembles the residual vector for the main unknowns, including computation of uhat, q, w, and AV fields, and both element and face integrals.
    - MPI variants (GetQMPI, RuResidualMPI, RuResidualMPI1): Support for distributed memory parallelism using MPI, with non-blocking communication for solution exchange between subdomains.
    - Residual: Top-level function to compute the residual, handling both serial and MPI-parallel cases, and applying time-stepping and debug output.
    - Enzyme AD support (GetdQ, GetdAv, dRuResidual, dRuResidualMPI, dResidual): Functions for computing directional derivatives of the residual using automatic differentiation (Enzyme), including MPI support.
    - ComputeQ: Utility to compute q variable, supporting both serial and MPI-parallel execution.

    The file relies on several helper functions and drivers for element/face integrals, matrix operations, and communication primitives.
    It supports both CUDA/HIP GPU acceleration and MPI parallelism, with conditional compilation for these features.

    Data structures:
    - solstruct: Solution variables (udg, wdg, odg, etc.)
    - resstruct: Residual vectors and intermediate storage
    - appstruct: Application-specific parameters
    - masterstruct: Master element data (shape functions, etc.)
    - meshstruct: Mesh connectivity and indexing
    - tempstruct: Temporary arrays for computation
    - commonstruct: Common parameters (problem size, flags, MPI info, etc.)

    Key features:
    - Modular design for element/face integrals and variable updates
    - Support for artificial viscosity smoothing and interpolation
    - MPI communication for distributed memory parallelism
    - Optional debug output to binary files
    - Automatic differentiation support for Jacobian-vector products

    Note: Many functions rely on external drivers and array utilities (e.g., ArrayExtract, ArrayInsert, PutFaceNodes, etc.) defined elsewhere in the codebase.
*/
#ifndef __RESIDUAL
#define __RESIDUAL

#include "geometry.cpp"
#include "massinv.cpp"
#include "qequation.cpp"
#include "wequation.cpp"
#include "uequation.cpp"
#include "qresidual.cpp"
#include "uresidual.cpp"
#include "getuhat.cpp"

void DG2CGAVField(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, dstype* avcg, dstype* avdg, Int backend)
{
    Int ncavdg = common.nco;
    Int ncavcg = common.nco;

    for(Int i = 0; i<common.ncAV; i++)
    { 
        //extract ith component of av field and store it in tmp.tempn
        ArrayExtract(tmp.tempn, avdg, common.npe, ncavdg, common.ne, 0, common.npe, i, i+1, 0, common.ne); 

        //make it a CG field and store in tmp.tempg
        ArrayDG2CG2(tmp.tempg, tmp.tempn, mesh.colent2elem, mesh.rowent2elem, common.ndofucg, common.npe);   

        // convert CG field to DG field   
        GetArrayAtIndex(tmp.tempn, tmp.tempg, mesh.cgelcon, common.npe*common.ne);      

        // insert cg field (tmp.tempn) into avcg
        ArrayInsert(avcg, tmp.tempn, common.npe, ncavcg, common.ne, 0, common.npe, i, i+1, 0, common.ne);    
    }
}

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

    // Element integrals
    RqElem(sol, res, app, master, mesh, tmp, common, handle, nbe1, nbe2, backend);
        
    START_TIMING;
    // Face integrals
    RqFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
    END_TIMING(16);       
        
    // elements in the range [e1, e2)
    Int e1 = common.eblks[3*nbe1]-1;    
    Int e2 = common.eblks[3*(nbe2-1)+1];
    Int f1 = common.fblks[3*nbf1]-1;
    Int f2 = common.fblks[3*(nbf2-1)+1];    

    //cout<<res.Rh[common.nf*common.npf*ncq-1]<<endl;
    //print3darray(res.Rh, common.npf, ncq, 5);
    // assemble face residual vector into element residual vector
    // PutFaceNodes(res.Rq, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
    //         mesh.ent2ind2, npf, npe, ncq, e1, e2, 0);    
    PutFaceNodes(res.Rq, &res.Rh[common.npf*common.ncq*f1],  mesh.facecon, npf, ncq, npe, ncq, f1, f2);
  
    if (common.wave==1)
        // get the source term due to the time derivative (for wave problem)  
        ArrayExtract(&res.Rq[N], sol.sdg, npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2);  
    else
        // set it to zero
        ArraySetValue(&res.Rq[N], zero, npe*ncq*(e2-e1));
        
    dstype scalar = one;
    if (common.wave==1)
        scalar = one/common.dtfactor;

    // START_TIMING;
    // // Apply the mass matrix inverse and the factor fc_q to the residual
    ApplyMinv(&res.Rq[N], res.Minv, &res.Rq[npe*ncq*e1], scalar, common.curvedMesh, npe, ncq, e1, e2);              
    // END_TIMING(17);   
                
    // START_TIMING;
    // // append q into udg
    ArrayInsert(sol.udg, &res.Rq[N], npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2);           
    // END_TIMING(18);       
}

void GetW(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
        Int nbe1, Int nbe2, Int nbf1, Int nbf2, Int backend)
{        
    if (common.subproblem==0) {
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
            ArrayExtract(tmp.tempn, sol.udg, common.npe, common.nc, common.ne1, 0, common.npe, 0, common.ncu, e1, e2);                                                  
            ArrayAXPBY(&sol.wdg[npe*ncw*e1], tmp.tempn, &sol.wsrc[npe*ncw*e1], scalar, scalar, npe*ncw*(e2-e1));                        
        }        
        else if (common.wave==0) {             
            if ((fabs(common.dae_alpha) < 1e-10) && (fabs(common.dae_beta) < 1e-10)) {
                // use Newton to solve the nonlinear system F(w, u) = 0 to obtain w for given u                
                for (int iter=0; iter<10; iter++) {
                  // evaluate nonlinear system F(w, u)
                  EosDriver(tmp.tempn, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1], 
                      &sol.wdg[npe*ncw*e1], mesh, master, app, sol, tmp, common, npe, e1, e2, backend);            
                  
                  int nn = npe*(e2-e1);                  
                  // check convergence
                  dstype nrm = NORM(common.cublasHandle, nn*ncw, tmp.tempn, backend);                   
                  if (nrm < 1e-6) break;                                       
                  
                  // compute jacobian matrix dF/dw
                  EosdwDriver(tmp.tempg, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1], 
                      &sol.wdg[npe*ncw*e1], mesh, master, app, sol, tmp, common, npe, e1, e2, backend);            
                                  
                  // compute the inverse of jacobian matrix
                  if (ncw==1) 
                    ArrayEosInverseMatrix11(tmp.tempg, npe, ncw, e2-e1);
                  else if (ncw==2)
                    ArrayEosInverseMatrix22(tmp.tempg, npe, ncw, e2-e1);
                  else if (ncw==3)
                    ArrayEosInverseMatrix33(tmp.tempg, npe, ncw, e2-e1);
                  else {
                    error("Equation of states functionality supports at most three dependent variables.");
                  }
                  
                  // perform dw = inverse(dF/dw) * F(w, u)
                  ArrayEosMatrixMultiplication(&tmp.tempn[nn*ncw], tmp.tempg, tmp.tempn, npe, ncw, e2-e1, 1);
                  
                  // update w = w - dw
                  ArrayAXPBY(&sol.wdg[npe*ncw*e1], &sol.wdg[npe*ncw*e1], &tmp.tempn[nn*ncw], one, minusone, nn*ncw);                    
                }                                
            }
            else {
                // alpha * dw/dt + beta w = sourcew(u,q,v)
                // calculate the source term Sourcew(xdg, udg, odg, wdg)
                SourcewDriver(tmp.tempn, &sol.xdg[npe*ncx*e1], &sol.udg[npe*nc*e1], &sol.odg[npe*nco*e1], 
                        &sol.wdg[npe*ncw*e1], mesh, master, app, sol, tmp, common, npe, e1, e2, backend);            
                if (common.dae_steps==0) { // alpha * dw/dt + beta w = sourcew(u,q,v)
                    // 1.0/(alpha*dirkd/dt + beta)
                    dstype scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta);

                    //cout<<common.dae_alpha<<" "<<common.dtfactor<<" "<<scalar<<" "<<app.physicsparam[0]<<endl;

                    // calculate w = (1/(alpha*dirkd/dt + beta))*(alpha*wsrc + Sourcew(xdg, udg, odg, wdg))  
                    ArrayAXPBY(&sol.wdg[npe*ncw*e1], &sol.wsrc[npe*ncw*e1], tmp.tempn, common.dae_alpha*scalar, scalar, npe*ncw*(e2-e1));                    
                }
                else {
                    // calculate tmp = alpha*wsrc + Sourcew(xdg, udg, odg, wdg)) 
                    ArrayAXPBY(tmp.tempn, &sol.wsrc[npe*ncw*e1], tmp.tempn, common.dae_alpha, one, npe*ncw*(e2-e1));                    

                    // 1.0/(alpha*dirkd/dt + beta + gamma)
                    dstype scalar = one/(common.dae_alpha*common.dtfactor + common.dae_beta + common.dae_gamma);

                    //cout<<scalar<<endl;

                    // calculate w = (1/(alpha*dirkd/dt + beta + gamma))*(gamma*walg + alpha*wsrc + Sourcew(xdg, udg, odg, wdg))  
                    ArrayAXPBY(&sol.wdg[npe*ncw*e1], &sol.wdual[npe*ncw*e1], tmp.tempn, common.dae_gamma*scalar, scalar, npe*ncw*(e2-e1));                                    
                }                
            }
        }
    }    
    }
}

void GetAv(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{
    AvfieldDriver(sol.odg, sol.xdg, sol.udg, sol.odg, sol.wdg, mesh, master, app, sol, tmp, common, backend); 
    // note that the args imply avfield can depend on odg...not sure this is true.
    // This is true actually; but we need to be careful with autodiff. Might need a seperate avfield...
    for (Int iav = 0; iav<common.AVsmoothingIter; iav++){
        DG2CGAVField(sol, res, app, master, mesh, tmp, common, sol.odg, sol.odg, backend);
    }

    for (Int j=0; j<common.nbe; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];                
        GetElemNodes(tmp.tempn, sol.odg, common.npe, common.nco, 
                0, common.nco, e1, e2);        
        Node2Gauss(common.cublasHandle, &sol.odgg[common.nge*common.nco*e1], 
            tmp.tempn, master.shapegt, common.nge, common.npe, (e2-e1)*common.nco, backend);        
    }         
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];            
        
        GetFaceNodes(tmp.tempn, sol.odg, mesh.facecon, common.npf, common.nco, 
                common.npe, common.nco, f1, f2, 1);          
        Node2Gauss(common.cublasHandle, &sol.og1[common.ngf*common.nco*f1], 
            tmp.tempn, master.shapfgt, common.ngf, common.npf, (f2-f1)*common.nco, backend);               
        
        GetFaceNodes(tmp.tempn, sol.odg, mesh.facecon, common.npf, common.nco, 
                common.npe, common.nco, f1, f2, 2);          
        Node2Gauss(common.cublasHandle, &sol.og2[common.ngf*common.nco*f1], 
            tmp.tempn, master.shapfgt, common.ngf, common.npf, (f2-f1)*common.nco, backend);               
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
    
    // compute w
    if (common.ncw>0)
         GetW(sol, res, app, master, mesh, tmp, common, handle, nbe1q, nbe2q, nbf1, nbf2, backend);                
    
    if (common.ncAV>0 && common.frozenAVflag == 0)
    /// FROM MEETING: MUST MOVE THE VOLUME RESIDUAL AND EVALUATE 0 to common.npe2 if unfrozen for mpi
        GetAv(sol, res, app, master, mesh, tmp, common, handle, backend);

    // Element integrals
    RuElem(sol, res, app, master, mesh, tmp, common, handle, nbe1u, nbe2u, backend);    
                    
    // Face integrals
    RuFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
        
    Int f1 = common.fblks[3*nbf1]-1;
    Int f2 = common.fblks[3*(nbf2-1)+1];        
    PutFaceNodes(res.Ru, &res.Rh[common.npf*common.ncu*f1],  mesh.facecon, common.npf, common.ncu, common.npe, common.ncu, f1, f2);            

    // Int e1 = common.eblks[3*nbe1u]-1;    
    // Int e2 = common.eblks[3*(nbe2u-1)+1];    
    
    // // assemble face residual vector into element residual vector
    // PutFaceNodes(res.Ru, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
    //         mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0);
                
    // change sign 
    // ArrayMultiplyScalar(res.Ru, minusone, common.ndof1, backend);        
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
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend);
    //for (n=0; n<common.nelemsend; n++)         
    //    ArrayCopy(&tmp.buffsend[bsz*n], &sol.udg[nudg*common.elemsend[n]], bsz, backend);            
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv);
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
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend);
    
    //for (n=0; n<common.nelemsend; n++)         
    //    ArrayCopy(&tmp.buffsend[bsz*n], &sol.udg[nudg*common.elemsend[n]], bsz, backend);                    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

#ifdef HAVE_HIP
    hipDeviceSynchronize();
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
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
    if (common.frozenAVflag == 1)
    // if AVfield is not part of residual, we can evaluate interior volume integrals before receving neighbor information
        RuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, backend);    
    END_TIMING(9);    
        
    START_TIMING; 
    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);

    END_TIMING(10);    

    START_TIMING; 
    /* copy buffrecv to udg */
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv);
 
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
    
    if (common.ncAV>0 && common.frozenAVflag == 0)
    /// FROM MEETING: MUST MOVE THE VOLUME RESIDUAL AND EVALUATE 0 to common.npe2 if unfrozen for mpi
        GetAv(sol, res, app, master, mesh, tmp, common, handle, backend);

    
    if (common.frozenAVflag == 1)
    { // calculate Ru for interface elements
        RuElem(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe1, backend); 
    } 
    else
    { // For unfrozen AV, we can only compute Ru for interior and interface after calculating AV 
        RuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe1, backend); 
    }
    
       
    END_TIMING(11);    
     
    START_TIMING; 
    // calculate Ru for all faces
    RuFace(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
    END_TIMING(12);    
        
    // assemble face residual vector into element residual vector
    Int f1 = common.fblks[3*0]-1;
    Int f2 = common.fblks[3*(common.nbf-1)+1];        
    PutFaceNodes(res.Ru, res.Rh,  mesh.facecon, common.npf, common.ncu, common.npe, common.ncu, f1, f2);            

    // Int e1 = common.eblks[3*0]-1;    
    // Int e2 = common.eblks[3*(common.nbe1-1)+1];       
    // PutFaceNodes(res.Ru, res.Rh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
    //         mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0);
    
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
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendind, bsz*common.nelemsend);
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
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
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv);
    
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
    ArrayMultiplyScalar(res.Ru, minusone, common.ndof1);      
        
    //common.dtfactor
    if (common.tdep==1) 
        ArrayMultiplyScalar(res.Ru, one/common.dtfactor, common.ndof1);                
    
    if (common.debugMode==1) {
        writearray2file(common.fileout + "_uh.bin", sol.uh, common.npf*common.ncu*common.nf, backend);
        writearray2file(common.fileout + "_udg.bin", sol.udg, common.npe*common.nc*common.ne, backend);
        writearray2file(common.fileout + "_Ru.bin", res.Ru, common.npe*common.ncu*common.ne, backend);
    }    
}

///////////////////////////////////////////////////////////////////////////////////////////
//// Calculate just dR(u)/du v, with u, q, uhat, R(u) already precalculated /////
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_ENZYME
void GetdQ(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
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

    // Element integrals
    dRqElem(sol, res, app, master, mesh, tmp, common, handle, nbe1, nbe2, backend);
        
    // Face integrals
    dRqFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
        
    // elements in the range [e1, e2)
    Int e1 = common.eblks[3*nbe1]-1;    
    Int e2 = common.eblks[3*(nbe2-1)+1];
       
    PutFaceNodes(res.dRq, res.dRh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, npf, npe, ncq, e1, e2, 0, backend);
    
    if (common.wave==1)
    { //TODO: not checked yet
        // get the source term due to the time derivative (for wave problem)  
        ArrayExtract(&res.dRq[N], sol.sdg, npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2);  
    }
    else
    {   // set it to zero
        ArraySetValue(&res.dRq[N], zero, npe*ncq*(e2-e1));
     }
    dstype scalar = one;
    if (common.wave==1)
        scalar = one/common.dtfactor;

    ApplyMinv(&res.dRq[N], res.Minv, &res.dRq[npe*ncq*e1], scalar, common.curvedMesh, npe, ncq, e1, e2);  

    ArrayInsert(sol.dudg, &res.dRq[N], npe, nc, ne, 0, npe, ncu, ncu+ncq, e1, e2);           
}

void GetdAv(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
        tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{
    AvfieldDriver(sol.odg, sol.dodg, sol.xdg, sol.udg, sol.dudg, sol.odg, sol.wdg, sol.dwdg, mesh, master, app, sol, tmp, common, backend); 

    for (Int iav = 0; iav<common.AVsmoothingIter; iav++){
        DG2CGAVField(sol, res, app, master, mesh, tmp, common, sol.dodg, sol.dodg, backend);
    }

    for (Int j=0; j<common.nbe; j++) {
        Int e1 = common.eblks[3*j]-1;
        Int e2 = common.eblks[3*j+1];                
        GetElemNodes(tmp.tempn, sol.dodg, common.npe, common.nco, 
                0, common.nco, e1, e2);        
        Node2Gauss(common.cublasHandle, &sol.dodgg[common.nge*common.nco*e1], 
            tmp.tempn, master.shapegt, common.nge, common.npe, (e2-e1)*common.nco, backend);        
    }         
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];            
        
        GetFaceNodes(tmp.tempn, sol.dodg, mesh.facecon, common.npf, common.nco, 
                common.npe, common.nco, f1, f2, 1);          
        Node2Gauss(common.cublasHandle, &sol.dog1[common.ngf*common.nco*f1], 
            tmp.tempn, master.shapfgt, common.ngf, common.npf, (f2-f1)*common.nco, backend);               
        
        GetFaceNodes(tmp.tempn, sol.dodg, mesh.facecon, common.npf, common.nco, 
                common.npe, common.nco, f1, f2, 2);          
        Node2Gauss(common.cublasHandle, &sol.dog2[common.ngf*common.nco*f1], 
            tmp.tempn, master.shapfgt, common.ngf, common.npf, (f2-f1)*common.nco, backend);               
    }       
}


void dRuResidual(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, 
   Int nbe1u, Int nbe2u, Int nbe1q, Int nbe2q, Int nbf1, Int nbf2, Int backend)
{   
    // compute (duhat/du v)
    // GetUhat(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);
    GetdUhat(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);

    // compute dq
    if (common.ncq>0)
    { 
        // GetQ(sol, res, app, master, mesh, tmp, common, handle, nbe1q, nbe2q, nbf1, nbf2, backend);
        GetdQ(sol, res, app, master, mesh, tmp, common, handle, nbe1q, nbe2q, nbf1, nbf2, backend);
    }               

    if (common.ncAV>0 && common.frozenAVflag == 0)
    /// FROM MEETING: MUST MOVE THE VOLUME RESIDUAL AND EVALUATE 0 to common.npe2 if unfrozen for mpi
        GetdAv(sol, res, app, master, mesh, tmp, common, handle, backend);


    //TODO: DAE functionality not implemented yet

    // Element integrals
    dRuElem(sol, res, app, master, mesh, tmp, common, handle, nbe1u, nbe2u, backend);    

    // Face integrals
    dRuFace(sol, res, app, master, mesh, tmp, common, handle, nbf1, nbf2, backend);

    Int e1 = common.eblks[3*nbe1u]-1;    
    Int e2 = common.eblks[3*(nbe2u-1)+1];    
    
    // insert face nodes into (dRu/du)v
    PutFaceNodes(res.dRu, res.dRh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0, backend);
}

#ifdef  HAVE_MPI

void dRuResidualMPI(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, 
   tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{      
    // non-blocking send and receive solutions on exterior and outer elements to neighbors    
    Int bsz = common.npe*common.ncu;
    Int nudg = common.npe*common.nc;
    Int n;
   
    // copy some portion of du to buffsend
    GetArrayAtIndex(tmp.buffsend, sol.dudg, mesh.elemsendind, bsz*common.nelemsend, backend);

#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif

#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.buffsend[psend], nsend, MPI_DOUBLE, neighbor, 0,
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    // non-blocking receive 
    Int nrecv, precv = 0;
    for (n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.buffrecv[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                   EXASIM_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
    // compute uhat for all faces
    GetdUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);            
    
    // calculate q for interior elements
    if (common.ncq>0)         
        GetdQ(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, 0, common.nbf, backend);        
    // TODO: DAE AD matvec

    // calculate Ru for interior elements
    if (common.frozenAVflag == 1)
    // if AVfield is not part of residual, we can evaluate interior volume integrals before receving neighbor information
        dRuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe0, backend);    
        
    // non-blocking receive solutions on exterior and outer elements from neighbors
    // wait until all send and receive operations are completely done
    MPI_Waitall(request_counter, common.requests, common.statuses);

    // copy buffrecv to udg
    PutArrayAtIndex(sol.dudg, tmp.buffrecv, mesh.elemrecvind, bsz*common.nelemrecv, backend);

    // compute uhat for all faces
    GetdUhat(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);

    // calculate q for interface and exterior elements
    if (common.ncq>0)         
        GetdQ(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe2, 0, common.nbf, backend);        
    // TODO: DAE AD Matvec    

    if (common.ncAV>0 && common.frozenAVflag == 0)
    /// FROM MEETING: MUST MOVE THE VOLUME RESIDUAL AND EVALUATE 0 to common.npe2 if unfrozen for mpi
        GetdAv(sol, res, app, master, mesh, tmp, common, handle, backend);
    
    if (common.frozenAVflag == 1)
    { // calculate Ru for interface elements
        dRuElem(sol, res, app, master, mesh, tmp, common, handle, common.nbe0, common.nbe1, backend); 
    } 
    else
    { // For unfrozen AV, we can only compute Ru for interior and interface after calculating AV 
        dRuElem(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbe1, backend); 
    }

    // calculate Ru for all faces
    dRuFace(sol, res, app, master, mesh, tmp, common, handle, 0, common.nbf, backend);
        
    // assemble face residual vector into element residual vector
    Int e1 = common.eblks[3*0]-1;    
    Int e2 = common.eblks[3*(common.nbe1-1)+1];       
    PutFaceNodes(res.dRu, res.dRh, mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, mesh.rowe2f2, mesh.cole2f2, 
            mesh.ent2ind2, common.npf, common.npe, common.ncu, e1, e2, 0, backend);
}


#endif

void dResidual(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common, cublasHandle_t handle, Int backend)
{    
    
    if (common.mpiProcs>1) { // mpi processes
#ifdef  HAVE_MPI        
        dRuResidualMPI(sol, res, app, master, mesh, tmp, common, handle, backend);
#endif        
    }
    else {
        dRuResidual(sol, res, app, master, mesh, tmp, common, handle,
                0, common.nbe, 0, common.nbe, 0, common.nbf, backend);
    } 
    // change sign for matrix-vector product
    ArrayMultiplyScalar(res.dRu, minusone, common.ndof1, backend);     
    if (common.tdep==1)
    { 
        ArrayMultiplyScalar(res.dRu, one/common.dtfactor, common.ndof1, backend); 
    }
    if (common.debugMode==1) {
        writearray2file(common.fileout + "enz_uh.bin", sol.uh, common.npf*common.ncu*common.nf, backend);
        writearray2file(common.fileout + "enz_udg.bin", sol.udg, common.npe*common.nc*common.ne, backend);
        writearray2file(common.fileout + "enz_Ru.bin", res.Ru, common.npe*common.ncu*common.ne, backend);
    }    
}

#endif

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


