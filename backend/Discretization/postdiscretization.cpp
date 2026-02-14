#ifndef __POSTDISCRETIZATION
#define __POSTDISCRETIZATION

#include "postdiscretization.h"
#include "ioutilities.cpp"

#include "../Model/KokkosDrivers.cpp"
#include "connectivity.cpp"
#include "readbinaryfiles.cpp"
#include "setstructs.cpp"
#include "residual.cpp"
#include "matvec.cpp"
#include "qoicalculation.cpp"

// Both CPU and GPU constructor
CDiscretization::CDiscretization(string filein, string fileout, string exasimpath, Int mpiprocs, Int mpirank, 
        Int fileoffset, Int omprank, Int backend, Int builtinmodelID, Int nsca, Int nvec, Int nten, Int nsurf, Int nvqoi) 
{
    common.backend = backend;
    common.exasimpath = exasimpath;
    //common.builtinmodelID = builtinmodelID;
    //app.builtinmodelID = builtinmodelID;

    if (backend>1) { // GPU
#ifdef HAVE_GPU        
        // host structs
        solstruct hsol;
        resstruct hres;
        appstruct happ;
        masterstruct hmaster; 
        meshstruct hmesh;
        tempstruct htmp;    
        commonstruct hcommon;     

        hcommon.backend = backend;        
        // allocate data for structs in CPU memory
        cpuInit(hsol, hres, happ, hmaster, hmesh, htmp, hcommon, filein, fileout, 
                mpiprocs, mpirank, fileoffset, omprank);                    
                
        // copy data from cpu memory to gpu memory
        gpuInit(sol, res, app, master, mesh, tmp, common, 
            hsol, hres, happ, hmaster, hmesh, htmp, hcommon);                
        app.read_uh = happ.read_uh;
        if (common.spatialScheme > 0)  { // HDG        
          TemplateMalloc(&mesh.bf, hcommon.nfe*hcommon.ne, 0);
          for (int i=0; i<hcommon.nfe*hcommon.ne; i++) mesh.bf[i] = hmesh.bf[i];   
        }

       // copy hsol.xcg to sol.xcg for paraview visualization
        sol.szxcg = hsol.szxcg;
        TemplateMalloc(&sol.xcg, sol.szxcg, 0);
        TemplateCopytoHost(sol.xcg, hsol.xcg, sol.szxcg, 0);
        if (common.mpiRank==0) printf("free CPU memory \n");
          
        // release CPU memory
        happ.freememory(1);        
        hmaster.freememory(1);        
        hmesh.freememory(1);        
        hsol.freememory(1);        
        htmp.freememory(1);        
        hres.freememory(1);        
        hcommon.freememory();             
#endif        
    }
    else  {// CPU
        cpuInit(sol, res, app, master, mesh, tmp, common, filein, fileout, 
                mpiprocs, mpirank, fileoffset, omprank);    
    }
    common.read_uh = app.read_uh;
    if (nsca > 0) common.nsca = nsca;
    if (nvec > 0) common.nvec = nvec;
    if (nten > 0) common.nten = nten;
    if (nsurf > 0) common.nsurf = nsurf;
    if (nvqoi > 0) common.nvqoi = nvqoi;     

    // compute the geometry quantities
    if (common.mpiRank==0) printf("start compGeometry... \n");
    compGeometry(backend);        
    if (common.mpiRank==0) printf("finish compGeometry... \n");        

    // compute the inverse of the mass matrix
    if (common.spatialScheme == 0) {
        if (common.mpiRank==0) printf("start compMassInverse... \n");
        compMassInverse(backend);    
        if (common.mpiRank==0) printf("finish compMassInverse... \n");        
    }
    
    // moved from InitSolution to here
    if ((common.ncq>0) && (common.wave==0) && (common.spatialScheme == 0)) evalQSer(backend); 
    
    if (common.spatialScheme > 0)  { // HDG
      Int neb = common.neb; // maximum number of elements per block
      Int npe = common.npe; // number of nodes on master element
      Int npf = common.npf; // number of nodes on master face
      Int nfe = common.nfe; // number of faces on master element
      Int ne = common.ne; // number of elements in this subdomain
      Int nf = common.nf; // number of faces in this subdomain
      Int ncx = common.ncx; // number of compoments of (xdg)
      Int nc = common.nc; // number of compoments of (u, q)
      Int ncu = common.ncu; // number of compoments of (u)
      Int ncq = common.ncq; // number of compoments of (q)      
      Int nbe = common.nbe; // number of blocks for elements
      int ncu12 = common.szinterfacefluxmap;
      
      if (common.mpiRank==0) 
        printf("Init HDG Discretization ... \n");        
      
      int nboufaces = 0; // number of boundary faces
      int maxbc = 0; // maximum number of boundary conditions
      for (int i=0; i<nfe*ne; i++) {
        if (mesh.bf[i] > 0) nboufaces++;
        maxbc = max(maxbc, mesh.bf[i]);
      }
      common.maxnbc = maxbc;      
                                      
      // compute uhat by getting u on faces
      if (!common.read_uh){
          if (common.mpiRank==0) 
              printf("===============================Constructing uh==========================\n");
          GetFaceNodes(sol.uh, sol.udg, mesh.f2e, mesh.perm, npf, ncu, npe, nc, nf);
      }

      if (common.ncq > 0) {                
        qEquation(sol, res, app, master, mesh, tmp, common, backend);      

        if (common.mpiRank==0) 
          printf("Finish qEquation ... \n");        

        // compute the flux q = -nabla u and store it in sol.udg
        if (common.wave == 0 && sol.szudg != npe*nc*ne) {
            hdgGetQ(sol.udg, sol.uh, sol, res, mesh, tmp, common, backend);
            if (common.mpiRank==0) printf("Finish hdgGetQ ... \n");     
        }        
      }
    }

    if (common.mpiRank==0) {
      if (common.debugMode==1) {
        common.printinfo();
        app.printinfo();
        res.printinfo();
        tmp.printinfo();
        sol.printinfo();
        mesh.printinfo();
        master.printinfo();
      }      
      printf("finish CDiscretization constructor... \n");        
    }
}

// destructor 
CDiscretization::~CDiscretization()
{        
    app.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: app memory is freed successfully.\n");
    master.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: master memory is freed successfully.\n");
    mesh.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: mesh memory is freeed successfully.\n");
    sol.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: sol memory is freed successfully.\n");
    tmp.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: tmp memory is freed successfully.\n");
    res.freememory(common.backend);
    if (common.mpiRank==0) printf("CDiscretization destructor: res memory is freed successfully.\n");
    common.freememory();
    if (common.mpiRank==0) printf("CDiscretization destructor: common memory is freed successfully.\n");

#ifdef HAVE_CUDA    
    if (common.backend==2) {
        CHECK(cudaEventDestroy(common.eventHandle));
        CHECK_CUBLAS(cublasDestroy(common.cublasHandle));
    }
#endif    
    
#ifdef HAVE_HIP    
    if (common.backend==3) {
        CHECK(hipEventDestroy(common.eventHandle));
        CHECK_HIPBLAS(hipblasDestroy(common.cublasHandle));
    }
#endif        
}

// Compute and store the geometry
void CDiscretization::compGeometry(Int backend) {
    if (common.mpiRank==0) printf("start ElemGeom... \n");
    ElemGeom(sol, master, mesh, tmp, common, common.cublasHandle, backend);   
    if (common.mpiRank==0) printf("Finish ElemGeom... \n");
    FaceGeom(sol, master, mesh, tmp, common, common.cublasHandle, backend);   

    if (common.spatialScheme>0)
      ElemFaceGeom(sol, master, mesh, tmp, common, common.cublasHandle, backend);   
}

// Compute and store the inverse of the mass matrix
void CDiscretization::compMassInverse(Int backend) {
    ComputeMinv(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
}

void CDiscretization::hdgAssembleLinearSystem(dstype *b, Int backend)
{
    int n = common.npe*common.ncu;
    int m = common.npf*common.nfe*common.ncu;
    int ne = common.ne1;

    ArraySetValue(res.H, zero, m*m*ne);
    ArraySetValue(res.Rh, zero, m*ne);
    ArraySetValue(res.Ru, zero, n*ne);
    ArraySetValue(res.F, zero, n*m*ne);    

#ifdef HAVE_MPI     
    hdgAssembleLinearSystemMPI(b, sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
#else    
    uEquationHDG(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
    hdgAssembleRHS(b, res.Rh, mesh, common);
#endif

    if (common.preconditioner==0) {
      // fix bug here: tmp.tempn is not enough memory to store ncu*npf*ncu*npf*nf 
      hdgBlockJacobi(res.K, res.H, res, mesh, tmp, common, common.cublasHandle, backend);      
    }
    else if (common.preconditioner==1) {
      hdgElementalAdditiveSchwarz(res.K, res.H, res, mesh, tmp, common, common.cublasHandle, backend);      
    }
    else if (common.preconditioner==2) {
      hdgBlockILU0(res.K, res.H, res, mesh, tmp, common, common.cublasHandle, backend);
    }
        
//     if (common.preconditioner==0) {
//       // assemble block Jacobi matrix from res.H using the FIRST elements in mesh.f2e
//       BlockJacobi(res.K, res.H, mesh.f2e, npf, nfe, ncu, common.nf);
// 
//       // assemble block Jacobi matrix from res.H using the SECOND elements in mesh.f2e
//       BlockJacobi(res.K, res.H, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf0);    
// 
//       // inverse block Jacobi matrix
//       //Inverse(handle, res.K, tmp.tempn, res.ipiv, ncf, common.nf, backend); // fix bug here
//       for (Int j=0; j<common.nbf; j++) {      
//         Int f1 = common.fblks[3*j]-1;
//         Int f2 = common.fblks[3*j+1];                  
//         Inverse(common.cublasHandle, &res.K[ncf*ncf*f1], tmp.tempn, res.ipiv, ncf, f2-f1, backend); 
//       }          
//     }
//     else if (common.preconditioner==1) { 
//       ArrayCopy(res.K, res.H, ncf*nfe*ncf*nfe*common.ne); 
//       ElementalAdditiveSchwarz(res.K, res.H, mesh.f2e, mesh.elemcon, npf, nfe, ncu, common.nf);       
// 
//       for (Int j=0; j<common.nbe; j++) {              
//         Int e1 = common.eblks[3*j]-1;
//         Int e2 = common.eblks[3*j+1];          
//         Inverse(common.cublasHandle, &res.K[ncf*common.nfe*ncf*common.nfe*e1], tmp.tempn, res.ipiv, ncf*nfe, e2-e1, backend); 
//       }          
//     }
//     else if (common.preconditioner==2) { // Block ILU0
//       //hdgBlockILU0(res.K, res.H, res, mesh, tmp, common, common.cublasHandle, backend);
//       Int nfse = common.nfse; // number of faces in each superelement
//       Int nse  = common.nse;  // number of superelements
//       Int N = nse*ncf*ncf;  
// 
//       AssembleBlockILU0(res.K, res.H, mesh.f2e, mesh.elemcon, mesh.face, mesh.row_ptr, mesh.col_ind, npf, nfe, ncu, nfse, nse);
// 
//       for (int i = 0; i < nfse; ++i) {
//           int diag_idx = common.ind_ii[i];
// 
//           // Invert all diagonal blocks at diag_idx, in-place (batched)
//           double *diag_blocks = &res.K[diag_idx * N];
//           Inverse(common.cublasHandle, diag_blocks, tmp.tempn, res.ipiv, ncf, nse, backend); 
// 
//           int nj = common.num_ji[i];
//           for (int j = 0; j < nj; ++j) {
//               int idx_ji = common.ind_ji[j + i * nj];
// 
//               // Multiply all nse blocks: block_ji = block_ji * block_diag, in-place
//               double *block_ji = &res.K[idx_ji * N];
//               PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, one, block_ji, ncf, diag_blocks, ncf, zero, tmp.tempn, ncf, nse, backend);             
//               ArrayCopy(block_ji, tmp.tempn, N);
// 
//               int nl = common.num_jl[j + i * nj];
//               for (int l = 0; l < nl; ++l) {
//                   int idx_jl = common.ind_jl[l + j * nl + i * nj * nl];
//                   int idx_il = common.ind_il[l + j * nl + i * nj * nl];
// 
//                   double *block_jl = &res.K[idx_jl * N];
//                   double *block_il = &res.K[idx_il * N];
//                   PGEMNMStridedBached(common.cublasHandle, ncf, 1, ncf, minusone, block_ji, ncf, block_il, ncf, one, block_jl, ncf, nse, backend);                             
//               }
//           }
//       }      
//     }    
}

void CDiscretization::hdgAssembleResidual(dstype *b, Int backend)
{
    int n = common.npe*common.ncu;
    int m = common.npf*common.nfe*common.ncu;
    int ne = common.ne1;
    ArraySetValue(res.Rh, zero, m*ne);
    ArraySetValue(res.Ru, zero, n*ne);

#ifdef HAVE_MPI     
    hdgAssembleResidualMPI(b, sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
#else    
    // b, K, H, F, Ru    
    ResidualHDG(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
    //uEquationHDG(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
    hdgAssembleRHS(b, res.Rh, mesh, common);      
#endif
}

// residual evaluation
void CDiscretization::evalResidual(Int backend)
{
    // compute the residual vector
    Residual(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
}

// residual evaluation
void CDiscretization::evalResidual(dstype* Ru, dstype* u, Int backend)
{ 
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1);  

    // compute the residual vector R(u)
    Residual(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);

    // copy the residual vector to Ru
    ArrayCopy(Ru, res.Ru, common.ndof1);
}

// q evaluation
void CDiscretization::evalQ(Int backend)
{
    // compute the flux q
    ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
}

void CDiscretization::evalQSer(Int backend)
{
    // compute the flux q    
    GetUhat(sol, res, app, master, mesh, tmp, common, common.cublasHandle, 0, common.nbf, backend);        
    GetQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, 0, common.nbe, 0, common.nbf, backend);        
}

void CDiscretization::evalQ(dstype* q, dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1);

    // compute the flux q
    ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);

    // get q from udg
    ArrayExtract(q, sol.udg, common.npe, common.nc, common.ne, 0, common.npe, 
            common.ncu, common.ncu+common.ncq, 0, common.ne1);
}

// matrix-vector product
void CDiscretization::evalMatVec(dstype* Jv, dstype* v, dstype* u, dstype* Ru, Int backend)
{    
    MatVec(Jv, sol, res, app, master, mesh, tmp, common, common.cublasHandle, v, u, Ru, backend); 
}

// matrix-vector product
void CDiscretization::evalMatVec(dstype* Jv, dstype* v, dstype* u, dstype* Ru, Int spatialScheme, Int backend)
{    
    if (spatialScheme == 0) {// LDG
      MatVec(Jv, sol, res, app, master, mesh, tmp, common, common.cublasHandle, v, u, Ru, backend); 
    }
    else if (spatialScheme == 1) { // HDG  
      hdgMatVec(Jv, res.H, v, res.Rh, res.Rq, res, app, mesh, common, tmp, common.cublasHandle, backend);
    }
}

void CDiscretization::updateUDG(dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1);

    if (common.ncq>0)
        // compute the flux q
        ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
}

void CDiscretization::updateU(dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1);
}

void CDiscretization::evalAVfield(dstype* avField, dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne);
    
    // compute the flux q
    if (common.ncq>0)        
        ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);

    // compute the av field
    AvfieldDriver(avField, sol.xdg, sol.udg, sol.odg, sol.wdg, mesh, master, app, sol, tmp, common, backend);    
}

void CDiscretization::evalAVfield(dstype* avField, Int backend)
{    
    
#ifdef  HAVE_MPI    
    Int bsz = common.npe*common.nc;
    Int nudg = common.npe*common.nc;
    Int n;
    
    /* copy some portion of u to buffsend */
    //for (n=0; n<common.nelemsend; n++)         
    //    ArrayCopy(&tmp.buffsend[bsz*n], &sol.udg[nudg*common.elemsend[n]], bsz, backend);           
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendudg, bsz*common.nelemsend);

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
                   EXASIM_COMM_LOCAL, &common.requests[request_counter]);
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
                   EXASIM_COMM_LOCAL, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);

    /* copy buffrecv to udg */
    //for (n=0; n<common.nelemrecv; n++) 
    //    ArrayCopy(&sol.udg[nudg*common.elemrecv[n]], &tmp.buffrecv[bsz*n], bsz, backend);        
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvudg, bsz*common.nelemrecv);
#endif
  
    // compute the av field
    AvfieldDriver(avField, sol.xdg, sol.udg, sol.odg, sol.wdg, mesh, master, app, sol, tmp, common, backend);           
}

void CDiscretization::evalOutput(dstype* output, Int backend)
{
#ifdef  HAVE_MPI    
    Int bsz = common.npe*common.nc;
    Int nudg = common.npe*common.nc;
    Int n;
    
    /* copy some portion of u to buffsend */
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendudg, bsz*common.nelemsend);

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
                   EXASIM_COMM_LOCAL, &common.requests[request_counter]);
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
                   EXASIM_COMM_LOCAL, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }

    // non-blocking receive solutions on exterior and outer elements from neighbors
    /* wait until all send and receive operations are completely done */
    MPI_Waitall(request_counter, common.requests, common.statuses);
        
    /* copy buffrecv to udg */
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvudg, bsz*common.nelemrecv);
#endif
            
    // compute the output field
    OutputDriver(output, sol.xdg, sol.udg, sol.odg, sol.wdg, mesh, master, app, sol, tmp, common, backend);    
//     void OutputDriver(dstype *f, dstype *xg, dstype *udg, dstype *odg, dstype *wdg, meshstruct &mesh, 
//         masterstruct &master, appstruct &app, solstruct &sol, tempstruct &temp, 
//         commonstruct &common, Int nge, Int e1, Int e2, Int backend)

}


void CDiscretization::evalMonitor(dstype* output,  dstype* udg, dstype* wdg, Int nc, Int backend)
{
    // compute the output field
    MonitorDriver(output, nc, sol.xdg, udg, sol.odg, wdg, mesh, master, app, sol, tmp, common, backend);    
}

void CDiscretization::DG2CG(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend)
{
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in utm
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne);
        
        // make it a CG field and store in res.Ru
        ArrayDG2CG(res.Ru, utm, mesh.cgent2dgent, mesh.rowent2elem, common.ndofucg);
        
        // convert CG field to DG field
        GetArrayAtIndex(utm, res.Ru, mesh.cgelcon, common.npe*common.ne);
        
        // insert utm into ucg
        ArrayInsert(ucg, utm, common.npe, ncucg, common.ne, 0, common.npe, i, i+1, 0, common.ne);
    }
}

void CDiscretization::DG2CG2(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend)
{
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in utm
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne);

        // make it a CG field and store in res.Ru
        ArrayDG2CG2(res.Ru, utm, mesh.colent2elem, mesh.rowent2elem, common.ndofucg, common.npe);
        
        // convert CG field to DG field
        GetArrayAtIndex(utm, res.Ru, mesh.cgelcon, common.npe*common.ne);
        
        // insert utm into ucg
        ArrayInsert(ucg, utm, common.npe, ncucg, common.ne, 0, common.npe, i, i+1, 0, common.ne);
    }
}

void CDiscretization::DG2CG3(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend)
{
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in utm
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne);
        
        // make it a CG field and store in res.Ru
        ArrayDG2CG(&ucg[i*common.ndofucg], utm, mesh.cgent2dgent, mesh.rowent2elem, common.ndofucg);
    }
}

Int CDiscretization::getFacesOnInterface(Int **faces, const Int boundarycondition)
{
    int nintfaces = getinterfacefaces(mesh.bf, common.nfe, common.ne1, boundarycondition);
    int *intfaces = nullptr; 
    TemplateMalloc(&intfaces, nintfaces, 0);

    getinterfacefaces(intfaces, mesh.bf, common.nfe, common.ne1, boundarycondition, nintfaces);

    TemplateMalloc(faces, nintfaces, common.backend);
    TemplateCopytoDevice(*faces, intfaces, nintfaces, common.backend);                           

    CPUFREE(intfaces);

    return nintfaces;
}

void CDiscretization::getDGNodesOnInterface(dstype* xdgint, const Int* faces, const Int nfaces)
{
    // npf * nfaces * ncx
    GetBoudaryNodes(xdgint, sol.xdg, faces, mesh.perm, common.nfe, 
                  common.npf, common.npe, common.ncx, common.ncx, nfaces);
}

void CDiscretization::getUDGOnInterface(dstype* udgint, const Int* faces, const Int nfaces)
{
    GetBoudaryNodes(udgint, sol.udg, faces, mesh.perm, common.nfe, 
                  common.npf, common.npe, common.nc, common.nc, nfaces);
}

void CDiscretization::getWDGOnInterface(dstype* wdgint, const Int* faces, const Int nfaces)
{
    GetBoudaryNodes(wdgint, sol.wdg, faces, mesh.perm, common.nfe, 
                  common.npf, common.npe, common.ncw, common.ncw, nfaces);
}

void CDiscretization::getODGOnInterface(dstype* odgint, const Int* faces, const Int nfaces)
{
    GetBoudaryNodes(odgint, sol.odg, faces, mesh.perm, common.nfe, 
                  common.npf, common.npe, common.nco, common.nco, nfaces);
}

void CDiscretization::getUHATOnInterface(dstype* uhint, const Int* faces, const Int nfaces)
{
    GetBoudaryNodes(uhint, sol.uh, faces, mesh.elemcon, common.nfe, 
                  common.npf, common.ncu, nfaces);
}

void CDiscretization::getNormalVectorOnInterface(dstype* nlint, dstype* xdgint, const Int nfaces)
{  
    Int nd = common.nd; 
    Int npf = common.npf; 
    Int nn = npf*nfaces; 
    Int ncx = common.ncx;    
    Int n2 = 0;    // jac
    Int n3 = nn;   // Jg
  
    if (nd==1) {
        FaceGeom1D(&tmp.tempn[n2], nlint, xdgint, nn);    
    }
    else if (nd==2){
        Node2Gauss(common.cublasHandle, &tmp.tempn[n3], xdgint, &master.shapfnt[npf*npf], npf, npf, nfaces*nd, common.backend);                
        FaceGeom2D(&tmp.tempn[n2], nlint, &tmp.tempn[n3], nn);
    }
    else if (nd==3) {
        Node2Gauss(common.cublasHandle, &tmp.tempn[n3], xdgint, &master.shapfnt[npf*npf], npf, npf, nfaces*nd, common.backend);                     
        Node2Gauss(common.cublasHandle, &tmp.tempn[n3+nn*nd], xdgint, &master.shapfnt[2*npf*npf], npf, npf, nfaces*nd, common.backend);                
        FaceGeom3D(&tmp.tempn[n2], nlint, &tmp.tempn[n3], nn);
    }
}

void CDiscretization::getFieldsAtGaussPointsOnInterface(dstype* xdggint, dstype* xdgint, const Int nfaces, const Int ncx)
{
    Node2Gauss(common.cublasHandle, xdggint, xdgint, master.shapfgt, common.ngf, common.npf, nfaces*ncx, common.backend);    
}

void CDiscretization::getInterfaceFluxesAtNodalPoints(dstype *flux, dstype* xdgint, dstype* nlint, const Int* faces, const Int nfaces)  
{    
    Int npf = common.npf;
    dstype *udgint = &tmp.tempn[0]; // reuse tempg for udgint
    dstype *odgint = &tmp.tempn[npf * nfaces * common.nc];
    dstype *wdgint = &tmp.tempn[npf * nfaces * common.nc + npf * nfaces * common.nco];
    dstype *uhint = &tmp.tempn[npf * nfaces * common.nc + npf * nfaces * common.nco + npf * nfaces * common.ncw];
    
    this->getUDGOnInterface(udgint, faces, nfaces);
    this->getODGOnInterface(odgint, faces, nfaces);
    this->getWDGOnInterface(wdgint, faces, nfaces);
    this->getUHATOnInterface(uhint, faces, nfaces);
    
    FintDriver(flux, xdgint, udgint, odgint, wdgint, uhint, nlint, mesh, 
        master, app, sol, tmp, common, nfaces*npf, 1, common.backend);        
}

void CDiscretization::getInterfaceFluxesAtGaussPoints(dstype *flux, dstype* xdggint, dstype* nlgint, const Int* faces, const Int nfaces)  
{    
    Int npf = common.npf;
    Int ngf = common.ngf;

    dstype *udgint = &tmp.tempn[0]; // reuse tempg for udgint
    dstype *odgint = &tmp.tempn[npf * nfaces * common.nc];
    dstype *wdgint = &tmp.tempn[npf * nfaces * common.nc + npf * nfaces * common.nco];
    dstype *uhint = &tmp.tempn[npf * nfaces * common.nc + npf * nfaces * common.nco + npf * nfaces * common.ncw];
    
    this->getUDGOnInterface(udgint, faces, nfaces);
    this->getODGOnInterface(odgint, faces, nfaces);
    this->getWDGOnInterface(wdgint, faces, nfaces);
    this->getUHATOnInterface(uhint, faces, nfaces);

    dstype *udggint = &tmp.tempg[0]; // reuse tempg2 for udggint
    dstype *odggint = &tmp.tempg[ngf * nfaces * common.nc];
    dstype *wdggint = &tmp.tempg[ngf * nfaces * (common.nc + common.nco)];
    dstype *uhgint = &tmp.tempg[ngf * nfaces * (common.nc + common.nco + common.ncw)];

    this->getFieldsAtGaussPointsOnInterface(udggint, udgint, nfaces, common.ncu);
    this->getFieldsAtGaussPointsOnInterface(odggint, odgint, nfaces, common.nco);
    this->getFieldsAtGaussPointsOnInterface(wdggint, wdgint, nfaces, common.ncw);
    this->getFieldsAtGaussPointsOnInterface(uhgint, uhint, nfaces, common.ncu);
    
    FintDriver(flux, xdggint, udggint, odggint, wdggint, uhgint, nlgint, mesh, 
        master, app, sol, tmp, common, nfaces*ngf, 1, common.backend);        
}

void CDiscretization::computeAverageSolutionsOnBoundary() 
{   
    if ( common.saveSolBouFreq>0 ) {
        for (Int j=0; j<common.nbf; j++) {
            Int ib = common.fblks[3*j+2];            
            if (ib == common.ibs) {     
                Int f1 = common.fblks[3*j]-1;
                Int f2 = common.fblks[3*j+1];                      
                Int npf = common.npf; // number of nodes on master face      
                Int npe = common.npe; // number of nodes on master face      
                Int nf = f2-f1;
                Int nn = npf*nf; 
                Int nc = common.nc; // number of compoments of (u, q, p)            
                Int ncu = common.ncu;
                Int ncw = common.ncw;
                GetArrayAtIndex(tmp.tempn, sol.udg, &mesh.findudg1[npf*nc*f1], nn*nc);                
                ArrayAXPBY(sol.bouudgavg, sol.bouudgavg, tmp.tempn, one, one, nn*nc);            
                ArrayAddScalar(&sol.bouudgavg[nn*nc], one, 1);
              
                if (common.spatialScheme==1)
                  GetFaceNodesHDG(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2);
                else
                  GetElemNodes(tmp.tempn, sol.uh, npf, ncu, 0, ncu, f1, f2);
                ArrayAXPBY(sol.bouuhavg, sol.bouuhavg, tmp.tempn, one, one, nn*ncu);            
                ArrayAddScalar(&sol.bouuhavg[nn*ncu], one, 1);                              

                if (ncw>0) {
                    GetFaceNodes(tmp.tempn, sol.wdg, mesh.facecon, npf, ncw, npe, ncw, f1, f2, 1);      
                    ArrayAXPBY(sol.bouwdgavg, sol.bouwdgavg, tmp.tempn, one, one, nn*ncw);            
                    ArrayAddScalar(&sol.bouwdgavg[nn*ncw], one, 1);                                  
                }
            }
        }                                        
    }
}

#endif        
