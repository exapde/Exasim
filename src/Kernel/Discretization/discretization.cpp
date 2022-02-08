#ifndef __DISCRETIZATION
#define __DISCRETIZATION

#ifdef HAVE_CUDA
#include "gpuDeviceInfo.cpp"
#endif

#include "discretization.h"
#include "errormsg.cpp"
#include "ioutilities.cpp"
#include "readbinaryfiles.cpp"
#include "setstructs.cpp"
#include "quadrature.cpp"
#include "residual.cpp"
#include "matvec.cpp"

// CPU constructor
CDiscretization::CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, 
        Int ompthreads, Int omprank) 
{
    cpuInit(sol, res, app, master, mesh, tmp, common, filein, fileout, 
            mpiprocs, mpirank, ompthreads, omprank);    
}

// GPU constructor
#ifdef HAVE_CUDA
CDiscretization::CDiscretization(CDiscretization& hdata)
{    
    gpuInit(sol, res, app, master, mesh, tmp, common, 
            hdata.sol, hdata.res, hdata.app, hdata.master, hdata.mesh, hdata.tmp, hdata.common);          
}
#endif

// Both CPU and GPU constructor
CDiscretization::CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, 
        Int ompthreads, Int omprank, Int backend) 
{
        //common.ncarray = new Int[nomodels]; 
        //sol.udgarray = new dstype*[nomodels]; // array of pointers pointing to udg
    
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
        
        // allocate data for structs in CPU memory
        cpuInit(hsol, hres, happ, hmaster, hmesh, htmp, hcommon, filein, fileout, 
                mpiprocs, mpirank, ompthreads, omprank);                    
                
        // copy data from cpu memory to gpu memory
        gpuInit(sol, res, app, master, mesh, tmp, common, 
            hsol, hres, happ, hmaster, hmesh, htmp, hcommon);                
                
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
    else // CPU
        cpuInit(sol, res, app, master, mesh, tmp, common, filein, fileout, 
                mpiprocs, mpirank, ompthreads, omprank);            
}

// destructor 
CDiscretization::~CDiscretization()
{    
    app.freememory(common.cpuMemory);
    master.freememory(common.cpuMemory);
    mesh.freememory(common.cpuMemory);
    sol.freememory(common.cpuMemory);
    tmp.freememory(common.cpuMemory);
    res.freememory(common.cpuMemory);
    common.freememory();
           
    // free ncarray and udgarray
    
    
#ifdef HAVE_CUDA    
    if (common.cpuMemory==0) {
        CHECK(cudaEventDestroy(common.eventHandle));
        CHECK_CUBLAS(cublasDestroy(common.cublasHandle));
    }
#endif    
}

// Compute and store the geometry
void CDiscretization::compGeometry(Int backend) {
    ElemGeom(sol, master, mesh, tmp, common, common.cublasHandle, backend);   
    FaceGeom(sol, master, mesh, tmp, common, common.cublasHandle, backend);   
}

// Compute and store the inverse of the mass matrix
void CDiscretization::compMassInverse(Int backend) {
    ComputeMinv(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);    
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
   
    //ArraySetValue(sol.udg, zero, common.npe*common.nc*common.ne, backend);
 
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1, backend);  

    // compute the residual vector R(u)
    Residual(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
// error("stops");

    // copy the residual vector to Ru
    ArrayCopy(Ru, res.Ru, common.ndof1, backend);        
    // printf("Ru: %f\n", Ru[0]);
    // error("stop");
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
            0, common.ncu, 0, common.ne1, backend);  

    // compute the flux q
    ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);

    // get q from udg
    ArrayExtract(q, sol.udg, common.npe, common.nc, common.ne, 0, common.npe, 
            common.ncu, common.ncu+common.ncq, 0, common.ne1, backend);  
}

// matrix-vector product
void CDiscretization::evalMatVec(dstype* Jv, dstype* v, dstype* u, dstype* Ru, Int backend)
{    
    MatVec(Jv, sol, res, app, master, mesh, tmp, common, common.cublasHandle, v, u, Ru, backend); 
}

void CDiscretization::updateUDG(dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1, backend);  

    if (common.ncq>0)
        // compute the flux q
        ComputeQ(sol, res, app, master, mesh, tmp, common, common.cublasHandle, backend);
}

void CDiscretization::updateU(dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne1, backend);  
}

void CDiscretization::evalAVfield(dstype* avField, dstype* u, Int backend)
{
    // insert u into udg
    ArrayInsert(sol.udg, u, common.npe, common.nc, common.ne, 0, common.npe, 
            0, common.ncu, 0, common.ne, backend);  
    
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
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendudg, bsz*common.nelemsend, backend); 

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
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvudg, bsz*common.nelemrecv, backend);
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
    GetArrayAtIndex(tmp.buffsend, sol.udg, mesh.elemsendudg, bsz*common.nelemsend, backend); 

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
    PutArrayAtIndex(sol.udg, tmp.buffrecv, mesh.elemrecvudg, bsz*common.nelemrecv, backend);
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
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne, backend);         
        
        // make it a CG field and store in res.Ru
        ArrayDG2CG(res.Ru, utm, mesh.cgent2dgent, mesh.rowent2elem, common.ndofucg, backend);        
        
        // convert CG field to DG field
        GetArrayAtIndex(utm, res.Ru, mesh.cgelcon, common.npe*common.ne, backend);        
        
        // insert utm into ucg
        ArrayInsert(ucg, utm, common.npe, ncucg, common.ne, 0, common.npe, i, i+1, 0, common.ne, backend);          
    }
}

void CDiscretization::DG2CG2(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend)
{
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in utm
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne, backend);         

        // make it a CG field and store in res.Ru
        ArrayDG2CG2(res.Ru, utm, mesh.colent2elem, mesh.rowent2elem, common.ndofucg, common.npe, backend);        
        
        // convert CG field to DG field
        GetArrayAtIndex(utm, res.Ru, mesh.cgelcon, common.npe*common.ne, backend);        
        
        // insert utm into ucg
        ArrayInsert(ucg, utm, common.npe, ncucg, common.ne, 0, common.npe, i, i+1, 0, common.ne, backend);          
    }
}

void CDiscretization::DG2CG3(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend)
{
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in utm
        ArrayExtract(utm, udg, common.npe, ncudg, common.ne, 0, common.npe, i, i+1, 0, common.ne, backend);         
        
        // make it a CG field and store in res.Ru
        ArrayDG2CG(&ucg[i*common.ndofucg], utm, mesh.cgent2dgent, mesh.rowent2elem, common.ndofucg, backend);                
    }
}

#endif        
