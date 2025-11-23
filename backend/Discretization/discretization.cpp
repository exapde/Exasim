/*
 * Discretization Module
 * =====================
 * This file implements the core routines for initializing and managing the discretization structures
 * used in the Exasim backend for both CPU and GPU architectures. It supports multiple spatial schemes
 * (LDG, HDG) and various preconditioners, and is designed for parallel execution with MPI and GPU acceleration.

 * Main Components:
 * ----------------
 * - crs_init: Initializes the compressed row storage (CRS) structures for superelement-based preconditioning.
 * - CDiscretization: Main class encapsulating all discretization data and operations.
 *   - Constructor: Initializes all data structures, allocates memory, and sets up geometry and solution fields.
 *   - Destructor: Releases all allocated resources and handles for both CPU and GPU.
 *   - compGeometry: Computes and stores geometric quantities for elements and faces.
 *   - compMassInverse: Computes and stores the inverse of the mass matrix.
 *   - hdgAssembleLinearSystem: Assembles the HDG linear system and applies the selected preconditioner.
 *   - hdgAssembleResidual: Assembles the HDG residual vector.
 *   - evalResidual: Evaluates the residual vector for the current solution.
 *   - evalQ: Computes the flux q for the current solution.
 *   - evalQSer: Serial evaluation of flux q for non-wave problems.
 *   - evalMatVec: Computes matrix-vector products for Jacobian-vector operations.
 *   - updateUDG/updateU: Updates the solution fields with new values.
 *   - evalAVfield: Computes artificial viscosity fields, with MPI support for distributed domains.
 *   - evalOutput: Computes output quantities, with MPI support for distributed domains.
 *   - evalMonitor: Computes monitoring quantities for the solution.
 *   - DG2CG/DG2CG2/DG2CG3: Converts DG fields to CG fields using various mapping strategies.

 * Features:
 * ---------
 * - Supports both CPU and GPU backends (CUDA/HIP).
 * - MPI parallelization for distributed memory architectures.
 * - Flexible memory management for host and device.
 * - Multiple preconditioners: Block Jacobi, Elemental Additive Schwarz, Superelement Additive Schwarz (ILU0).
 * - Handles both LDG and HDG spatial schemes.
 * - Modular design for geometry, solution, and output computations.

 * Usage:
 * ------
 * Instantiate CDiscretization with appropriate input files and parallelization parameters.
 * Use member functions to assemble systems, evaluate residuals, compute fluxes, and manage solution fields.

 * Note:
 * -----
 * This file includes several implementation files (.cpp) directly for modularity and to support template-based
 * memory management and device/host operations.
 */
#ifndef __DISCRETIZATION
#define __DISCRETIZATION

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

#include "discretization.h"
#include "ioutilities.cpp"

#ifdef HAVE_TEXT2CODE
#include "../Model/ModelDrivers.cpp"
#else
#include "../Model/KokkosDrivers.cpp"
#endif

#include "../Model/modeltemplate.hpp"

#include "connectivity.cpp"
#include "readbinaryfiles.cpp"
#include "setstructs.cpp"
#include "residual.cpp"
#include "matvec.cpp"
#include "qoicalculation.cpp"

void crs_init(commonstruct& common, meshstruct& mesh, int *elem, int nse, int nese)
{            
    common.nse = nse;
    common.nese = nese;
    
    int *row_ptr = NULL; 
    int *col_ind = NULL; 
    int *face = NULL; 
    int *f2eelem = NULL; 
    int *f2e = NULL; 
    TemplateMalloc(&f2e, 4*common.nf, 0);
    TemplateCopytoHost(f2e, mesh.f2e, 4*common.nf, common.backend); 
    
    int nfelem = crs_faceordering(&row_ptr, &col_ind, &face, &f2eelem, elem, f2e, common.nse, common.nese, common.nfe, common.nf);

    common.nfse = nfelem;
    common.nnz = row_ptr[common.nfse];      
        
    int n = 2*(common.nfe-1);
    TemplateMalloc(&common.ind_ii, nfelem, 0);
    TemplateMalloc(&common.ind_ji, nfelem*n, 0);
    TemplateMalloc(&common.ind_jl, nfelem*n*n, 0);
    TemplateMalloc(&common.ind_il, nfelem*n*n, 0);
    TemplateMalloc(&common.num_ji, nfelem, 0);
    TemplateMalloc(&common.num_jl, nfelem*n, 0);
    TemplateMalloc(&common.Lind_ji, nfelem*n*2, 0);
    TemplateMalloc(&common.Uind_ji, nfelem*n*2, 0);
    TemplateMalloc(&common.Lnum_ji, nfelem*2, 0);
    TemplateMalloc(&common.Unum_ji, nfelem*3, 0);
    for (int i=0; i<nfelem; i++) common.ind_ii[i] = -1;
    for (int i=0; i<nfelem*n; i++) common.ind_ji[i] = -1;
    for (int i=0; i<nfelem*n*n; i++) common.ind_jl[i] = -1;
    for (int i=0; i<nfelem*n*n; i++) common.ind_il[i] = -1;
    for (int i=0; i<nfelem; i++) common.num_ji[i] = 0;
    for (int i=0; i<nfelem*n; i++) common.num_jl[i] = 0;
    for (int i=0; i<nfelem*n*2; i++) common.Lind_ji[i] = -1;
    for (int i=0; i<nfelem*n*2; i++) common.Uind_ji[i] = -1;
    for (int i=0; i<nfelem*2; i++) common.Lnum_ji[i] = 0;
    for (int i=0; i<nfelem*3; i++) common.Unum_ji[i] = 0;    
    
    crs_indexingilu0(common.ind_ii, common.ind_ji, common.ind_jl, common.ind_il, common.num_ji, common.num_jl, 
            common.Lind_ji, common.Uind_ji, common.Lnum_ji, common.Unum_ji, row_ptr, col_ind, common.nfe, nfelem);

//     print2iarray(f2eelem, common.nfe, nfelem);
//     print2iarray(row_ptr, 1, nfelem+1);
//     print2iarray(col_ind, 1, row_ptr[nfelem]);    
//     print2iarray(common.ind_ii, 1, nfelem);
//     print2iarray(common.ind_ji, n, nfelem);
//     print2iarray(common.ind_jl, n*n, nfelem);
//     print2iarray(common.ind_il, n*n, nfelem);
//     print2iarray(common.num_ji, 1, nfelem);
//     print2iarray(common.num_jl, n, nfelem);
//     print2iarray(common.Lind_ji, n*2, nfelem);
//     print2iarray(common.Uind_ji, n*2, nfelem);
//     print2iarray(common.Lnum_ji, 2, nfelem);
//     print2iarray(common.Unum_ji, 3, nfelem);

    TemplateMalloc(&mesh.row_ptr, nfelem+1, common.backend);
    TemplateMalloc(&mesh.col_ind, row_ptr[nfelem],common.backend);
    TemplateMalloc(&mesh.face, nse*nfelem, common.backend);
    TemplateCopytoDevice(mesh.row_ptr, row_ptr, nfelem+1, common.backend);                       
    TemplateCopytoDevice(mesh.col_ind, col_ind, row_ptr[nfelem], common.backend);    
    TemplateCopytoDevice(mesh.face, face, nse*nfelem, common.backend);      
    
//     writearray2file(common.fileout + "elem.bin", elem, nse*nese, 0);
//     writearray2file(common.fileout + "f2e.bin", f2e, 4*common.nf, 0);
//     
//     writearray2file(common.fileout + "ind_ii.bin", common.ind_ii, nfelem, 0);
//     writearray2file(common.fileout + "ind_ji.bin", common.ind_ji, n*nfelem, 0);
//     writearray2file(common.fileout + "ind_jl.bin", common.ind_jl, n*n*nfelem, 0);
//     writearray2file(common.fileout + "ind_il.bin", common.ind_il, n*n*nfelem, 0);
//     writearray2file(common.fileout + "num_ji.bin", common.num_ji, nfelem, 0);
//     writearray2file(common.fileout + "num_jl.bin", common.num_jl, n*nfelem, 0);
//     writearray2file(common.fileout + "Lind_ji.bin", common.Lind_ji, 2*n*nfelem, 0);
//     writearray2file(common.fileout + "Uind_ji.bin", common.Uind_ji, 2*n*nfelem, 0);
//     writearray2file(common.fileout + "Lnum_ji.bin", common.Lnum_ji, 2*nfelem, 0);
//     writearray2file(common.fileout + "Unum_ji.bin", common.Unum_ji, 3*nfelem, 0);
//     
//     writearray2file(common.fileout + "row_ptr.bin", mesh.row_ptr, nfelem+1, common.backend);
//     writearray2file(common.fileout + "col_ind.bin", mesh.col_ind, row_ptr[nfelem], common.backend);
//     writearray2file(common.fileout + "face.bin", mesh.face, nse*nfelem, common.backend);
    
    CPUFREE(row_ptr);
    CPUFREE(col_ind);
    CPUFREE(face);
    CPUFREE(f2eelem);
    CPUFREE(f2e);
}
      
// Both CPU and GPU constructor
CDiscretization::CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, 
        Int fileoffset, Int omprank, Int backend) 
{
    common.backend = backend;

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
      
      if (common.coupledboundarycondition>0) {
        common.nintfaces = getinterfacefaces(mesh.bf, common.eblks, nbe, nfe, common.coupledboundarycondition);
        int *intfaces = nullptr; // store interface faces
        TemplateMalloc(&intfaces, common.nintfaces, 0);
        getinterfacefaces(intfaces, mesh.bf, common.eblks, nbe, nfe, common.coupledboundarycondition, common.nintfaces);
        TemplateMalloc(&mesh.intfaces, common.nintfaces, common.backend);
        TemplateCopytoDevice(mesh.intfaces, intfaces, common.nintfaces, common.backend);                       
        mesh.szintfaces = common.nintfaces;
        
        CPUFREE(intfaces);
        
        TemplateMalloc(&sol.xdgint, ncx*npf*common.nintfaces, common.backend);
        GetBoudaryNodes(sol.xdgint, sol.xdg, mesh.intfaces, mesh.perm, nfe, npf, npe, ncx, ncx, common.nintfaces);
        sol.szxdgint = ncx*npf*common.nintfaces;                
        //print2darray(sol.xdgint, npf*common.nintfaces, ncx);
      }

      if (common.mpiRank==0) 
        printf("Maximum number of boundary conditions = %d \n", maxbc);        

      // print2iarray(mesh.bf, nfe, ne);
      // print2iarray(common.eblks, 3, nbe);

      int *boufaces = nullptr; // store boundary faces
      TemplateMalloc(&common.nboufaces, 1 + maxbc*nbe, 0);
      TemplateMalloc(&boufaces, nboufaces, 0);
      getboundaryfaces(common.nboufaces, boufaces, mesh.bf, common.eblks, nbe, nfe, maxbc, nboufaces);
      TemplateMalloc(&mesh.boufaces, nboufaces, common.backend);
      TemplateCopytoDevice(mesh.boufaces, boufaces, nboufaces, common.backend);                       
      mesh.szboufaces = nboufaces;

      CPUFREE(boufaces);
      //CPUFREE(mesh.bf);            
                          
      if ((common.preconditioner==2) && (common.szcartgridpart > 0)) {              
        if (common.cartgridpart[0]==2) {          
          int *elem = NULL;                
          int nse  = gridpartition2d(&elem, common.cartgridpart[1], common.cartgridpart[2], common.cartgridpart[3], common.cartgridpart[4], common.cartgridpart[5]);       
          int nese = common.cartgridpart[3]*common.cartgridpart[4];    
          crs_init(common, mesh, elem, nse, nese);
          CPUFREE(elem);
        }
        else if (common.cartgridpart[0]==3) {
          int *elem = NULL;   
          int nse  = gridpartition3d(&elem, common.cartgridpart[1], common.cartgridpart[2], common.cartgridpart[3], common.cartgridpart[4], common.cartgridpart[5], common.cartgridpart[6], common.cartgridpart[7]);       
          int nese = common.cartgridpart[4]*common.cartgridpart[5]*common.cartgridpart[6];      
          crs_init(common, mesh, elem, nse, nese);
          CPUFREE(elem);
        }                               
      }
      
      res.szH = npf*nfe*ncu*npf*nfe*ncu*common.ne; // HDG elemental matrices     
      res.szK = (npe*ncu*npe*ncu + npe*ncu*npe*ncq + npf*nfe*ncu*npe*ncq + npf*nfe*ncu*npe*ncu)*neb;          
      if (common.preconditioner==0)      // Block Jacobition preconditioner
        res.szP = ncu*npf*ncu*npf*nf;
      else if (common.preconditioner==1) // Elemental additive Schwarz preconditioner
        res.szP = npf*nfe*ncu*npf*nfe*ncu*common.ne;        
      else if (common.preconditioner==2) // Superelement additive Schwarz preconditioner
        res.szP = npf*ncu*npf*ncu*common.nse*common.nnz;        
      res.szV = ncu*npf*nf*(common.gmresRestart+1); // Krylov vectors in GMRES
      res.szK = max(res.szK, res.szP + res.szV);              
      res.szF = npe*ncu*npf*nfe*ncu*common.ne;      
      res.szipiv = max(max(npf*nfe,npe)*ncu*neb, ncu*npf*common.nfb);
            
      TemplateMalloc(&res.H, res.szH, backend);
      TemplateMalloc(&res.K, res.szK, backend);      
      TemplateMalloc(&res.F, res.szF, backend);
      TemplateMalloc(&res.ipiv, res.szipiv, backend); // fix big here     
            
      // B, D, G, K share the same memmory block 
      // It is also used for storing both the preconditioner matrix and sys.v
      res.D = &res.K[npf*nfe*ncu*npe*ncu*neb];
      res.B = &res.K[npf*nfe*ncu*npe*ncu*neb + npe*ncu*npe*ncu*neb];
      res.G = &res.K[npf*nfe*ncu*npe*ncu*neb + npe*ncu*npe*ncu*neb + npe*ncu*npe*ncq*neb];        
      
      if (common.coupledinterface>0) {
        res.szRi = npf*ncu12*common.ncie;
        res.szKi = npf*ncu12*npe*ncu*common.ncie;
        res.szHi = npf*ncu12*npf*nfe*ncu*common.ncie;
        TemplateMalloc(&res.Ri, res.szRi, backend);
        TemplateMalloc(&res.Ki, res.szKi, backend);
        TemplateMalloc(&res.Hi, res.szHi, backend);
      }
      
      if (common.mpiRank==0) 
        printf("Memory allocation ...\n");        

// #ifdef HAVE_CUDA
//       int n = npe*ncu;
//       int batchSize = neb;
//       TemplateMalloc(&res.ipiv, n * batchSize * sizeof(Int), backend);
//       TemplateMalloc(&res.info,  batchSize * sizeof(Int), backend);     

//       dstype **Dp_h = (dstype **)malloc(batchSize*sizeof(dstype *));
//       cudaMalloc(&res.Dptr, batchSize*sizeof(dstype *));
//       Dp_h[0] = res.D;
//       for (Int i = 1; i < batchSize; i++)
//         Dp_h[i] = Dp_h[i-1]+(n*n);
//       cudaMemcpy(res.Dptr,Dp_h,batchSize*sizeof(dstype *),cudaMemcpyHostToDevice);
      
//       dstype **Dinvp_h = (dstype **)malloc(batchSize*sizeof(dstype *));      
//       cudaMalloc(&res.Dinvptr,batchSize*sizeof(dstype *));
//       Dinvp_h[0] = tmp.tempn;
//       for (Int i = 1; i < batchSize; i++)
//         Dinvp_h[i] = Dinvp_h[i-1] + (n*n);
//       cudaMemcpy(res.Dinvptr, Dinvp_h, batchSize*sizeof(dstype *),cudaMemcpyHostToDevice);          

//       free(Dp_h);
//       free(Dinvp_h);
// #endif

      // compute uhat by getting u on faces
        // std::cout <<"app.read_uh in discretization.cpp is : " << common.read_uh<<endl;
      if (!common.read_uh){
          if (common.mpiRank==0) 
              printf("===============================Constructing uh==========================\n");
          GetFaceNodes(sol.uh, sol.udg, mesh.f2e, mesh.perm, npf, ncu, npe, nc, nf);
      }
      else {
          if (common.mpiRank==0) 
              printf("================================Reading uh==============================\n");
          // print2darray(sol.uh, npf*ncu*10, 2);
      }

      if (common.mpiRank==0) 
        printf("Finish GetFaceNodes ... \n");        

      // print3darray(sol.udg, common.npe, common.nc, common.ne);
      // print2darray(sol.uh, common.ncu*common.npf, common.nf);
      // GetFaceNodes(tmp.tempg, sol.xdg, mesh.f2e, mesh.perm, common.npf, common.nd, common.npe, common.nd, common.nf);
      // print3darray(tmp.tempg, common.nd, common.npf, common.nf);

      if (common.ncq > 0) {        
        if (common.coupledinterface>0) {
          res.szGi = npf*ncu12*npe*ncq*common.ncie;          
          TemplateMalloc(&res.Gi, res.szGi, backend);
        }
        
        // compute M^{-1} * C and store it in res.C
        // compute M^{-1} * E and store it in res.E
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

#endif        
