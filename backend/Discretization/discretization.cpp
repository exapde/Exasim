#ifndef __DISCRETIZATION
#define __DISCRETIZATION

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

#include "discretization.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "../../build/model/KokkosDrivers.cpp"
#include "readbinaryfiles.cpp"
#include "setstructs.cpp"
#include "residual.cpp"
#include "matvec.cpp"

// Both CPU and GPU constructor
CDiscretization::CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, 
        Int ompthreads, Int omprank, Int backend) 
{
    common.backend = backend;

    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
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
                mpiprocs, mpirank, ompthreads, omprank);                    
                
        // copy data from cpu memory to gpu memory
        gpuInit(sol, res, app, master, mesh, tmp, common, 
            hsol, hres, happ, hmaster, hmesh, htmp, hcommon);                
                
        if (common.spatialScheme > 0)  { // HDG        
          TemplateMalloc(&mesh.bf, hcommon.nfe*hcommon.ne, 0);
          for (int i=0; i<hcommon.nfe*hcommon.ne; i++) mesh.bf[i] = hmesh.bf[i];   
        }
        
        if (common.mpiRank==0) printf("free CPU memory \n");
          
        // release CPU memory
        happ.freememory(hcommon.cpuMemory);        
        hmaster.freememory(hcommon.cpuMemory);        
        hmesh.freememory(hcommon.cpuMemory);        
        hsol.freememory(hcommon.cpuMemory);        
        htmp.freememory(hcommon.cpuMemory);        
        hres.freememory(hcommon.cpuMemory);        
        hcommon.freememory();             
#endif        
    }
    else  {// CPU
        cpuInit(sol, res, app, master, mesh, tmp, common, filein, fileout, 
                mpiprocs, mpirank, ompthreads, omprank);    
    }

    // compute the geometry quantities
    if (common.mpiRank==0) printf("start compGeometry... \n");
    compGeometry(backend);        
    if (common.mpiRank==0) printf("finish compGeometry... \n");        

    // compute the inverse of the mass matrix
    if (common.spatialScheme == 0) compMassInverse(backend);    
    
    // moved from InitSolution to here
    if ((common.ncq>0) & (common.wave==0) & (common.spatialScheme == 0) ) evalQSer(backend); 
    
    if (common.spatialScheme > 0)  { // HDG
      Int neb = common.neb; // maximum number of elements per block
      Int npe = common.npe; // number of nodes on master element
      Int npf = common.npf; // number of nodes on master face
      Int nfe = common.nfe; // number of faces on master element
      Int ne = common.ne; // number of elements in this subdomain
      Int nf = common.nf; // number of faces in this subdomain
      Int nc = common.nc; // number of compoments of (u, q)
      Int ncu = common.ncu; // number of compoments of (u)
      Int ncq = common.ncq; // number of compoments of (q)      
      Int nbe = common.nbe; // number of blocks for elements

      if (common.mpiRank==0) 
        printf("Init HDG Discretization ... \n");        

      int nboufaces = 0; // number of boundary faces
      int maxbc = 0; // maximum number of boundary conditions
      for (int i=0; i<nfe*ne; i++) {
        if (mesh.bf[i] > 0) nboufaces++;
        maxbc = max(maxbc, mesh.bf[i]);
      }
      common.maxnbc = maxbc;      

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

      if (common.mpiRank==0) 
        printf("Number of boundary faces = %d \n", nboufaces);        

      // print1iarray(common.nboufaces, 1 + maxbc*nbe);
      // print1iarray(mesh.boufaces, nboufaces);      
      CPUFREE(boufaces);
      CPUFREE(mesh.bf);

      TemplateMalloc(&res.D, npe*ncu*npe*ncu*neb, backend);
      TemplateMalloc(&res.F, npe*ncu*npf*nfe*ncu*common.ne, backend);
      TemplateMalloc(&res.K, max(npf*nfe*ncu*npe*ncu*neb, ncu*npf*ncu*npf*nf), backend);
      TemplateMalloc(&res.H, npf*nfe*ncu*npf*nfe*ncu*common.ne, backend);
      TemplateMalloc(&res.ipiv, npe * ncu * neb * sizeof(Int), backend);      

      res.szD = npe*ncu*npe*ncu*neb;
      res.szF = npe*ncu*npf*nfe*ncu*common.ne;
      res.szK = max(npf*nfe*ncu*npe*ncu*neb, ncu*npf*ncu*npf*nf);
      res.szH = npf*nfe*ncu*npf*nfe*ncu*common.ne;
      res.szipiv = npe * ncu * neb;

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
      GetFaceNodes(sol.uh, sol.udg, mesh.f2e, mesh.perm, npf, ncu, npe, nc, nf);

      if (common.mpiRank==0) 
        printf("Finish GetFaceNodes ... \n");        

      // print3darray(sol.udg, common.npe, common.nc, common.ne);
      // print2darray(sol.uh, common.ncu*common.npf, common.nf);
      // GetFaceNodes(tmp.tempg, sol.xdg, mesh.f2e, mesh.perm, common.npf, common.nd, common.npe, common.nd, common.nf);
      // print3darray(tmp.tempg, common.nd, common.npf, common.nf);

      if (common.ncq > 0) {        
        TemplateMalloc(&res.B, npe*ncu*npe*ncq*neb, backend);
        TemplateMalloc(&res.G, npf*nfe*ncu*npe*ncq*neb, backend);
        res.szB = npe*ncu*npe*ncq*neb;
        res.szG = npf*nfe*ncu*npe*ncq*neb;

        // compute M^{-1} * C and store it in res.C
        // compute M^{-1} * E and store it in res.E
        qEquation(sol, res, app, master, mesh, tmp, common, backend);      

        if (common.mpiRank==0) 
          printf("Finish qEquation ... \n");        

        // compute the flux q = -nabla u and store it in sol.udg
        if (common.wave == 0) hdgGetQ(sol.udg, sol.uh, sol, res, mesh, tmp, common, backend);
        //print3darray(sol.udg, npe, nc, common.ne);        
        //error("stop");

        if (common.mpiRank==0) 
          printf("Finish hdgGetQ ... \n");        
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
    app.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: app memory is freed successfully.\n");
    master.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: master memory is freed successfully.\n");
    mesh.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: mesh memory is freeed successfully.\n");
    sol.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: sol memory is freed successfully.\n");
    tmp.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: tmp memory is freed successfully.\n");
    res.freememory(common.cpuMemory);
    if (common.mpiRank==0) printf("CDiscretization destructor: res memory is freed successfully.\n");
    common.freememory();
    if (common.mpiRank==0) printf("CDiscretization destructor: common memory is freed successfully.\n");

#ifdef HAVE_CUDA    
    if (common.cpuMemory==0) {
        CHECK(cudaEventDestroy(common.eventHandle));
        CHECK_CUBLAS(cublasDestroy(common.cublasHandle));
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
    // fix bug here: tmp.tempn is not enough memory to store ncu*npf*ncu*npf*nf 
    hdgBlockJacobi(res.K, res.H, tmp.tempn, res.ipiv, mesh, common, common.cublasHandle, backend);      
#endif

    // dstype nrmUDG = PNORM(common.cublasHandle, common.npe*common.nc*common.ne1, sol.udg, backend);  
    // dstype nrmUH = PNORM(common.cublasHandle, common.npf*common.ncu*common.nf, sol.uh, backend);
    // cout<<"ne1 = "<<common.ne1<<", nrmUDG="<<nrmUDG<<", nrmUH="<<nrmUH<<endl;
    // dstype nrmRh = PNORM(common.cublasHandle, m*common.ne1, res.Rh, backend);    
    // dstype nrmH = PNORM(common.cublasHandle, m*m*common.ne1, res.H, backend);    
    // cout<<"ne1 = "<<common.ne1<<", nrmRh="<<nrmRh<<", nrmH="<<nrmH<<endl;
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
      hdgMatVec(Jv, res.H, v, res.Rh, res.Rq, mesh, common, tmp, common.cublasHandle, backend);
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

