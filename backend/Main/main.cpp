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

#ifdef _ENZYME
#define HAVE_ENZYME
#endif

#ifdef _MUTATIONPP
#define HAVE_MPP
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

#ifdef HAVE_MPP
#include <mutation++.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef TIMING
#include <chrono>
#endif

#include <Kokkos_Core.hpp>

using namespace std;

#include "../Common/common.h"        // declaration of variables and structs
#include "../Common/cpuimpl.h"       // interface to cpu functions
#include "../Common/kokkosimpl.h"    // interface to kokkos functions
#include "../Common/pblas.h"         // wrappers for blas libaries and MPI     

#include "../Discretization/discretization.cpp" // discretization class
#include "../Preconditioning/preconditioner.cpp" // preconditioner class
#include "../Solver/solver.cpp"                 // solver class
#include "../Solution/solution.cpp"             // solution class

int main(int argc, char** argv) 
{   
  Kokkos::initialize(argc, argv);
  {        
    if( argc >= 4 ) {
    }
    else {      
      printf("Usage: ./cppfile nummodels InputFile OutputFile\n");
      return 1;
    }                
    
    Int nummodels, restart, mpiprocs, mpirank, shmrank, ncores, nthreads, backend;    
    string mystr = string(argv[1]);
    nummodels = stoi(mystr);  // number of pde models
    string filein[nummodels]; 
    string fileout[nummodels];
    
    // two-physics and two-domain problems  
    int mpiprocs0 = 0;
    if (nummodels>100) {
      mpiprocs0 = nummodels - 100;
      nummodels = 2;
    }
    
    for (int i=0; i<nummodels; i++) {
        filein[i]  = string(argv[2*i+2]); // input files
        fileout[i]  = string(argv[2*i+3]); // output files        
        //cout<<filein[i]<<endl;
        //cout<<fileout[i]<<endl;
    }
    
    restart = 0;
    if (argc>=(2*nummodels+3)) {
        mystr = string(argv[2*nummodels+2]);
        restart = stoi(mystr);
    }             
    
    // reset nummodels
    if (mpiprocs0 > 0) nummodels = 1;
    
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
      
    if ((mpiprocs0 > 0) && (mpiprocs<= mpiprocs0)) {
      printf("For two-domain problem, total number of MPI processors (%d) must be greater than # MPI processors on the first domain (%d)\n", mpiprocs, mpiprocs0);
      return 1;
    }
      
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
    
    Int fileoffset = 0;
    Int gpuid = 0;            
      
    // initialize PDE models
    CSolution** pdemodel = new CSolution*[nummodels];     
    // initialize file streams
    ofstream* out = new ofstream[nummodels];    
                
    for (int i=0; i<nummodels; i++) {        
        if (mpiprocs0==0) {
          pdemodel[i] = new CSolution(filein[i], fileout[i], mpiprocs, mpirank, fileoffset, gpuid, backend);                 
        }
        else if (mpiprocs0 > 0) {
          if (mpirank < mpiprocs0) { 
            pdemodel[i] = new CSolution(filein[0], fileout[0], mpiprocs, mpirank, fileoffset, gpuid, backend);       
          }
          else {
            fileoffset = mpiprocs0;
            pdemodel[i] = new CSolution(filein[1], fileout[1], mpiprocs, mpirank, fileoffset, gpuid, backend);       
          }
        }
        pdemodel[i]->disc.common.nomodels = nummodels;
        
        // can move these to the constructor 
        pdemodel[i]->disc.common.ncarray = new Int[nummodels]; 
        pdemodel[i]->disc.sol.udgarray = new dstype*[nummodels]; // array of pointers pointing to udg
        
        if (pdemodel[i]->disc.common.timestepOffset>0)
            restart = pdemodel[i]->disc.common.timestepOffset;    
        
        if (restart>0) {
            pdemodel[i]->disc.common.timestepOffset = restart;
            pdemodel[i]->disc.common.time = restart*pdemodel[i]->disc.common.dt[0];            
        }        
        
        if (pdemodel[i]->disc.common.mpiRank==0 && pdemodel[i]->disc.common.saveResNorm==1) {
            string filename = pdemodel[i]->disc.common.fileout + "_residualnorms" + NumberToString(i) + ".bin";
            out[i].open(filename.c_str(), ios::out | ios::binary);            
            if (!out[i]) error("Unable to open file " + filename);        
        }            
    }            
        
    // set ncarray and udgarray for each PDE model
    for (int i=0; i<nummodels; i++) 
        for (int j=0; j<nummodels; j++) {
            pdemodel[i]->disc.common.ncarray[j] = pdemodel[j]->disc.common.nc;            
            pdemodel[i]->disc.sol.udgarray[j] = &pdemodel[j]->disc.sol.udg[0]; // model[i], model[j]
        }    
   
//     if (pdemodel[0]->disc.common.subproblem==1) { // subproblem
//         for (int i=0; i<nummodels; i++) {
//             if (restart>0) {
//                 pdemodel[i]->disc.common.currentstep = -1;
//                 pdemodel[i]->ReadSolutions(backend);  
//             }            
//         }      
//         pdemodel[0]->SolveProblem(pdemodel[1], backend);               
//     }
//     else 
    
    
    if ((pdemodel[0]->disc.common.tdep==1) && (pdemodel[0]->disc.common.runmode==0)) {
                        
        // initialize 
        for (int i=0; i<nummodels; i++) {
            if (restart>0) {
                pdemodel[i]->disc.common.currentstep = -1;
                pdemodel[i]->ReadSolutions(backend);  
            }            
            pdemodel[i]->InitSolution(backend);             
        }
        
        // initial time
        dstype time = pdemodel[0]->disc.common.time;           
        
        // time stepping with DIRK schemes
        for (Int istep=0; istep<pdemodel[0]->disc.common.tsteps; istep++)            
        {            
            for (int i=0; i<nummodels; i++) {
                // current timestep        
                pdemodel[i]->disc.common.currentstep = istep;

                // store previous solutions to calculate the source term        
                PreviousSolutions(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.common, backend);
            }
                                
            // compute the solution at the next step
            for (Int j=0; j<pdemodel[0]->disc.common.tstages; j++) {     
                
                if (pdemodel[0]->disc.common.mpiRank==0)
                    printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",istep+1,j+1,time + pdemodel[0]->disc.common.dt[istep]*pdemodel[0]->disc.common.DIRKcoeff_t[j]);                                
                
                for (int i=0; i<nummodels; i++) {
                    // current timestage
                    pdemodel[i]->disc.common.currentstage = j;

                    // current time
                    pdemodel[i]->disc.common.time = time + pdemodel[i]->disc.common.dt[istep]*pdemodel[i]->disc.common.DIRKcoeff_t[j];
        
                    // update source term             
                    UpdateSource(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.app, pdemodel[i]->disc.res, pdemodel[i]->disc.common, backend);

                    // solve the problem 
                    pdemodel[i]->SteadyProblem(out[i], backend);                             

                    // update solution 
                    UpdateSolution(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.app, pdemodel[i]->disc.res, pdemodel[i]->disc.tmp, pdemodel[i]->disc.common, backend);
                    //UpdateSolution(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.common, backend);                     
                }
            }
        
            for (int i=0; i<nummodels; i++) {
                //compute time-average solution
                if (pdemodel[i]->disc.common.compudgavg == 1) {
                    ArrayAXPBY(pdemodel[i]->disc.sol.udgavg, pdemodel[i]->disc.sol.udgavg, pdemodel[i]->disc.sol.udg, one, one, pdemodel[i]->disc.common.ndofudg1);            
                    ArrayAddScalar(&pdemodel[i]->disc.sol.udgavg[pdemodel[i]->disc.common.ndofudg1], one, 1);
                }

                // save solutions into binary files                
                pdemodel[i]->SaveSolutions(backend); 
                pdemodel[i]->SaveSolutionsOnBoundary(backend); 
                if (pdemodel[i]->disc.common.nce>0)
                    pdemodel[i]->SaveOutputCG(backend);                
            }
            
            // update time
            time = time + pdemodel[0]->disc.common.dt[istep];                    
        }                   
    }
    else {        
        for (int i=0; i<nummodels; i++) {                                
            if (pdemodel[i]->disc.common.runmode==0) {
                if (restart>0) {
                    pdemodel[i]->disc.common.currentstep = -1;
                    pdemodel[i]->ReadSolutions(backend);  
                }
                pdemodel[i]->SolveProblem(out[i], backend);           
            }
            else if (pdemodel[i]->disc.common.runmode==1){
                pdemodel[i]->disc.common.currentstep = -1;
                pdemodel[i]->ReadSolutions(backend);   
                if (pdemodel[i]->disc.common.ncq>0)
                    pdemodel[i]->disc.evalQ(backend);
                pdemodel[i]->disc.common.saveSolOpt = 1;
                pdemodel[i]->SaveSolutions(backend);      
                pdemodel[i]->SaveOutputCG(backend);            
            }
            else if (pdemodel[i]->disc.common.runmode==2){
                for (Int istep=0; istep<pdemodel[i]->disc.common.tsteps; istep++)            
                {                                            
                    if (((istep+1) % pdemodel[i]->disc.common.saveSolFreq) == 0)  
                    {
                        pdemodel[i]->disc.common.currentstep = istep;
                        pdemodel[i]->ReadSolutions(backend);                    
                        if (pdemodel[i]->disc.common.ncq>0)
                            pdemodel[i]->disc.evalQ(backend);
                        pdemodel[i]->SaveOutputCG(backend);            
                    }
                }
            }    
            else if (pdemodel[i]->disc.common.runmode==3){
                for (Int istep=0; istep<pdemodel[i]->disc.common.tsteps; istep++)            
                {                                            
                    if (((istep+1) % pdemodel[i]->disc.common.saveSolFreq) == 0)  
                    {
                        pdemodel[i]->disc.common.currentstep = istep;
                        pdemodel[i]->ReadSolutions(backend);                    
                        pdemodel[i]->SaveOutputCG(backend);            
                    }
                }
            }
            else if (pdemodel[i]->disc.common.runmode==4) {
                pdemodel[i]->disc.common.currentstep = -1;
                pdemodel[i]->ReadSolutions(backend);
                pdemodel[i]->SolveProblem(out[i], backend);   
            }
            else if (pdemodel[i]->disc.common.runmode==5){
                pdemodel[i]->disc.evalQSer(backend);
                pdemodel[i]->disc.common.saveSolOpt = 1;
                string filename = pdemodel[i]->disc.common.fileout + "_np" + NumberToString(pdemodel[i]->disc.common.mpiRank) + ".bin";                    
                //string filename = pdemodel[i]->disc.common.fileout + "_np" + NumberToString(pdemodel[i]->disc.common.mpiRank) + ".bin";
                writearray2file(filename, pdemodel[i]->disc.sol.udg, pdemodel[i]->disc.common.ndofudg1, backend);   
            }                
        }
    }
    
    for (int i=0; i<nummodels; i++) {                
        if (pdemodel[i]->disc.common.mpiRank==0 && pdemodel[i]->disc.common.saveResNorm==1)             
            out[i].close();                                
    }
    
    // After using pdemodel, delete each CSolution object and then the array of pointers
    for (int i = 0; i < nummodels; i++) {
        delete[] pdemodel[i]->disc.common.ncarray; // Delete the ncarray
        delete[] pdemodel[i]->disc.sol.udgarray; // Delete the udgarray
        delete pdemodel[i]; // Delete each CSolution object
    }

    delete[] pdemodel; // Delete the array of pointers
    delete[] out; // Delete the array of ofstream objects

#ifdef HAVE_MPI
  MPI_Finalize();
#endif         
    
  }
  Kokkos::finalize();  
  return 0;             
     
}
