#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
//#include <stdlib.h>
#include <chrono>


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
    
#ifdef HAVE_MPI    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes    
    MPI_Comm_size(MPI_COMM_WORLD, &mpiprocs);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
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
        printf("Using %d CPU cores per MPI process \n",ncores);
        printf("Using %d threads per MPI process \n",nthreads);         
#ifdef HAVE_CUDA            
        printf("Using the following GPU per MPI process \n"); 
        gpuDeviceInfo();
#endif                
    }                           
    
    Int ngpus = 0;
    Int gpuid = 0;        
    if (backend==2)
        ngpus = 1;
        
    CSolution Sol(filein, fileout, mpiprocs, mpirank, ngpus, gpuid, backend);           
    Sol.disc.compMassInverse(backend);
     
    auto begin = std::chrono::high_resolution_clock::now();   
    
    for (int i = 0; i<100; i++)
        Sol.disc.evalQ(backend);
    
    auto end = std::chrono::high_resolution_clock::now();
    double tim = std::chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
    printf("Computing time: %g miliseconds\n", tim);
    
//     Sol.disc.evalQ(backend);
//     writearray2file("UDG_np0.bin", Sol.disc.sol.udg, Sol.disc.common.ndofudg, backend);
    
//    Sol.SolveProblem(backend);        
        
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
