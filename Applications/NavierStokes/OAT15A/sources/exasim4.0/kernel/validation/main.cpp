#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
//#include <stdlib.h>
//#include <chrono>
#include <sys/unistd.h>


#ifdef _OPENMP
#define HAVE_OPENMP
#else
#define HAVE_ONETHREAD
#endif

#ifdef _CUDA
#define HAVE_CUDA
#endif

#ifdef _MPI
#define HAVE_MPI
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_CUDA    
#include <cuda_runtime.h>
#include <cuda.h>
#include "cublas_v2.h"
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef TIMING
#include <chrono>
#endif

using namespace std;

#include "../common/common.h"     // declaration of variables and structs
#include "../common/core.h"       // interface to core library
#include "../common/pblas.h"      // wrappers for blas libaries and MPI     
#include "../common/application.h"// interface to application library

#include "../discretization/discretization.cpp" // discretization class
#include "../preconditioner/preconditioner.cpp" // preconditioner class
#include "../solver/solver.cpp"                 // solver class
#include "../solution/solution.cpp"             // solution class

int main(int argc, char** argv) 
{   
    if( argc == 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }                
    
    string filein  = string(argv[1]); // input files
    string fileout  = string(argv[2]); // output files           
    Int mpiprocs, mpirank, ncores, nthreads, backend;    
    int compnodes = 1;
    
#ifdef HAVE_MPI    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes    
    MPI_Comm_size(MPI_COMM_WORLD, &mpiprocs);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    int shmrank;
    MPI_Comm_rank(shmcomm, &shmrank);
    printf("Using %dth local process \n",shmrank);
    
    // Determine the number of compute nodes
//     char names[mpiprocs][MPI_MAX_PROCESSOR_NAME];
//     int len;
//     MPI_Get_processor_name( names[mpirank], &len );
//     MPI_Allgather( MPI_IN_PLACE, 0, 0, names[0], MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD );
//     int indexes[mpiprocs];    
//     for ( int i = 0; i < mpiprocs; i++ ) {
//         int found = 0;
//         for ( int j = 0; j < compnodes; j++ ) {
//             if ( strncmp( names[i], names[indexes[j]], MPI_MAX_PROCESSOR_NAME ) == 0 ) {
//                 found = 1;
//                 break;
//             }
//         }
//         if ( found == 0 ) {
//             indexes[compnodes++] = i;
//         }
//     }    
#else        
    // get number of MPI processes and MPI ranks
    mpiprocs = 1;
    mpirank = 0;
#endif                
        
#ifdef HAVE_OPENMP    
    // set OpenMP threads
    ncores = omp_get_num_procs();
    nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);            
#else
    ncores = 1; 
    nthreads = 1;    
#endif
    
#ifdef HAVE_OPENMP
    backend=1; // multi-thread
#else
    backend=0; // single thread
#endif
#ifdef HAVE_CUDA  // CUDA          
    backend=2;
#endif
        
    if (mpirank==0) { 
        printf("Solving problem %s with backend = %d!!!\n\n",filein.c_str(),backend);                
        printf("Using %d MPI processes \n",mpiprocs);
        printf("Using %d compute nodes \n",compnodes);
        printf("Using %d CPU cores per MPI process \n",ncores);
        printf("Using %d threads per MPI process \n",nthreads);    
    }
    
#ifdef HAVE_CUDA            
        //printf("Using the following GPU per MPI process \n"); 
        //gpuDeviceInfo(mpirank);#ifdef HAVE_CUDA            
        int procspernode = mpiprocs/compnodes;
        int localrank = mpirank%procspernode;
        int device;
        char hostname[128];
        printf("Using the following GPU per MPI process \n");
        gpuDeviceInfo(localrank);
        cudaGetDevice( &device );
        gethostname( hostname, sizeof(hostname) );
        int cpu = sched_getcpu();
        printf("MPI rank %d using GPU %d and CPU %d on host %s\n", mpirank, device, cpu, hostname );
        size_t available, total;
        cudaMemGetInfo(&available, &total);
        cout<<"Available GPU Memory: "<<available<<" Total GPU Memory: "<<total<<endl;
#endif                           
    
    Int ngpus = 0;
    Int gpuid = 0;        
    if (backend==2)
        ngpus = 1;
        
    CSolution Sol(filein, fileout, mpiprocs, mpirank, ngpus, gpuid, backend);       
    
//     Sol.disc.compMassInverse(backend);
//     Sol.disc.evalQ(backend);
//     writearray2file("UDG_np0.bin", Sol.disc.sol.udg, Sol.disc.common.ndofudg, backend);
    
    if (Sol.disc.common.runmode==0) {
        Sol.SolveProblem(backend);   
    }
    else if (Sol.disc.common.runmode==1){
        // compute the geometry quantities
        Sol.disc.compGeometry(backend);

        // compute the inverse of the mass matrix
        Sol.disc.compMassInverse(backend);

        // compute q        
        Sol.disc.evalQSer(backend);
        
//         if (mpirank==0)
//             printf("Comute Mass Inverse...\n");
//         Sol.disc.compMassInverse(backend);
//         if (mpirank==0)
//             printf("Comute q...\n");
//         Sol.disc.evalQ(backend);
//         if (mpirank==0)
//             printf("Save solution...\n");
        
        Sol.disc.common.saveSolOpt = 1;
        SaveSolutions(Sol.disc.sol, Sol.solv.sys, Sol.disc.common, backend);            
    }
    
    //Sol.SolveProblem(backend);        
        
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
