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

#include "../Common/common.h"     // declaration of variables and structs
#include "../Common/core.h"       // interface to core library
#include "../Common/pblas.h"      // wrappers for blas libaries and MPI     
#include "../Common/application.h"// interface to application library

#include "../Discretization/discretization.cpp" // discretization class
#include "../Preconditioning/preconditioner.cpp" // preconditioner class
#include "../Solver/solver.cpp"                 // solver class
#include "../Solution/solution.cpp"             // solution class


void opuStgHomoTurb2D10(dstype *up, dstype *xdg, dstype *stgdata, dstype *uc, dstype t, int M, int N)
{                                 
    // loop over each grid point
    for (int m=0; m<M; m++) {    
        
        dstype un = 0.0;
        dstype vn = 0.0;

        // loop over wavenumbers
        for (int n=0; n<N; n++) {
            // wavenumber at mode n
            dstype kn  = stgdata[0*N+n];
            dstype amp = stgdata[1*N+n];

            // random numbers
            dstype phi = stgdata[2*N+n]; 
            dstype dx = stgdata[3*N+n];
            dstype dy = stgdata[4*N+n]; 
            dstype sigmax  = stgdata[6*N+n];
            dstype sigmay  = stgdata[7*N+n]; 
            dstype omega = stgdata[9*N+n];
            
            // angle
            dstype an = kn*(dx*(xdg[0*M+m]-uc[0]*t) + dy*(xdg[1*M+m]-uc[1]*t)) + phi + omega*t;
            // Fourier mode
            dstype bn = amp*cos(an);

            // fluctuating velocity field
            un = un + sigmax*bn;
            vn = vn + sigmay*bn;
        }

        up[0*M+m] = un;
        up[1*M+m] = vn;          
    }        
}


int main(int argc, char** argv) 
{   
    if( argc >= 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }                
    
    string filein  = string(argv[1]); // input files
    string fileout  = string(argv[2]); // output files           
    Int restart, mpiprocs, mpirank, shmrank, ncores, nthreads, backend;    
    //int compnodes = 50;
 
    restart = 0;
    if (argc>=4) {
        string mystr = string(argv[3]);
        restart = stoi(mystr);
    }             
    
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
    
    MPI_Comm_rank(shmcomm, &shmrank);
#else        
    // get number of MPI processes and MPI ranks
    mpiprocs = 1;
    mpirank = 0;
    shmrank = 0;
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
          
    if (backend==2) {
        if (mpirank==0) 
            printf("Using %d processors to solve the problem on GPU platform...\n", mpiprocs);
    }
    else {
        if (mpirank==0) 
            printf("Using %d processors to solve the problem on CPU platform...\n", mpiprocs);
    }
    
#ifdef HAVE_CUDA            
    int device;    
    cudaSetDevice(shmrank); 
    //gpuDeviceInfo(shmrank);
    cudaGetDevice( &device );
//     char hostname[128];
//     gethostname( hostname, sizeof(hostname) );
//     int cpu = sched_getcpu();    
//     printf("MPI rank %d using GPU %d and CPU %d on host %s\n", mpirank, device, cpu, hostname );
    size_t available, total;
    cudaMemGetInfo(&available, &total);
    cout<<"Available GPU Memory: "<<available<<" Total GPU Memory: "<<total<<endl;
#endif                           
    
    Int ngpus = 0;
    Int gpuid = 0;        
    if (backend==2)
        ngpus = 1;
       
    CSolution Solobj(filein, fileout, mpiprocs, mpirank, ngpus, gpuid, backend);       
    
    if (restart>0) {
        Solobj.disc.common.timestepOffset = restart;
        Solobj.disc.common.time = restart*Solobj.disc.common.dt[0];
        if (Solobj.disc.common.mpiRank==0)
            cout<<"Restart the simulation at time step: "<<restart<<endl;        
    }    
    
    if (Solobj.disc.common.runmode==0) {
        if (restart>0) {
            Solobj.disc.common.currentstep = -1;
            Solobj.ReadSolutions(backend);  
        }
        Solobj.SolveProblem(backend);           
    }
    else if (Solobj.disc.common.runmode==1){
        Solobj.disc.compGeometry(backend);
        Solobj.disc.compMassInverse(backend);
        Solobj.disc.common.currentstep = -1;
        Solobj.ReadSolutions(backend);   
        if (Solobj.disc.common.ncq>0)
            Solobj.disc.evalQ(backend);
        Solobj.disc.common.saveSolOpt = 1;
        Solobj.SaveSolutions(backend);      
        Solobj.SaveOutputCG(backend);            
    }
    else if (Solobj.disc.common.runmode==2){
        Solobj.disc.compGeometry(backend);
        Solobj.disc.compMassInverse(backend);
        for (Int istep=0; istep<Solobj.disc.common.tsteps; istep++)            
        {                                            
            if (((istep+1) % Solobj.disc.common.saveSolFreq) == 0)  
            {
                Solobj.disc.common.currentstep = istep;
                Solobj.ReadSolutions(backend);                    
                if (Solobj.disc.common.ncq>0)
                    Solobj.disc.evalQ(backend);
                Solobj.SaveOutputCG(backend);            
            }
        }
    }    
    else if (Solobj.disc.common.runmode==3){
        Solobj.disc.compGeometry(backend);
        Solobj.disc.compMassInverse(backend);                        
        for (Int istep=0; istep<Solobj.disc.common.tsteps; istep++)            
        {                                            
            if (((istep+1) % Solobj.disc.common.saveSolFreq) == 0)  
            {
                Solobj.disc.common.currentstep = istep;
                Solobj.ReadSolutions(backend);                    
                Solobj.SaveOutputCG(backend);            
            }
        }
    }
    else if (Solobj.disc.common.runmode==4) {
        Solobj.disc.common.currentstep = -1;
        Solobj.ReadSolutions(backend);
        Solobj.SolveProblem(backend);   
    }
    else if (Solobj.disc.common.runmode==5){
        Solobj.disc.compGeometry(backend);
        Solobj.disc.compMassInverse(backend);
        Solobj.disc.evalQSer(backend);
        Solobj.disc.common.saveSolOpt = 1;
        string filename = Solobj.disc.common.fileout + "_np" + NumberToString(Solobj.disc.common.mpiRank) + ".bin";                    
        writearray2file(filename, Solobj.disc.sol.udg, Solobj.disc.common.ndofudg1, backend);   
    }        
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
