#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string.h>
//#include <stdlib.h>
//#include <chrono>
//#include <sys/unistd.h>

#ifdef _OPENMP
#define HAVE_OPENMP
#else
#define HAVE_ONETHREAD
#endif

#ifdef _CUDA
#define HAVE_GPU
#define HAVE_CUDA
#endif

#ifdef _HIP
#define HAVE_GPU
#define HAVE_HIP
#endif

#ifdef _TEXT2CODE
#define HAVE_TEXT2CODE
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

#ifdef HAVE_HIP
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <hipblas/hipblas.h>
#include <rocblas/rocblas.h>
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
  
    Int nummodels, restart, mpiprocs, mpirank, shmrank, backend;    
    
#ifdef HAVE_MPI    
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

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
  
#ifdef HAVE_TEXT2CODE
    if (mpirank==0) printf("TEXT2CODE ENABLED...\n");    
#endif
    
// #ifdef HAVE_OPENMP    
//     // set OpenMP threads
//     ncores = omp_get_num_procs();
//     nthreads = omp_get_max_threads();
//     omp_set_num_threads(nthreads);            
// #else
//     ncores = 1; 
//     nthreads = 1;    
// #endif
    
#ifdef HAVE_OPENMP
    backend=1; // multi-thread
#else
    backend=0; // single thread
#endif
#ifdef HAVE_CUDA  // CUDA          
    backend=2;
#endif
#ifdef HAVE_HIP  // HIP          
    backend=3;
#endif
            
    if (backend==2) {
        if (mpirank==0) 
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on CUDA platform...\n", mpiprocs, mpirank, shmrank);
    }
    else if (backend==3) {
        if (mpirank==0) 
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on HIP platform...\n", mpiprocs, mpirank, shmrank);
    }
    else {
        if (mpirank==0) 
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on CPU platform...\n", mpiprocs, mpirank, shmrank);
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
    
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    
    std::cout << "MPI Rank: " << mpirank << ", Device Count: " << deviceCount << ", CUDA Device: " << device 
              << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB" 
              << std::endl;
#endif                           
    
#ifdef HAVE_HIP
    int device;
    // Set the device based on the rank
    CHECK(hipSetDevice(shmrank));

    // Get the current device ID
    CHECK(hipGetDevice(&device));

    // Query memory information for the current device
    size_t available, total;
    CHECK(hipMemGetInfo(&available, &total));

    int deviceCount = 0;
    hipGetDeviceCount(&deviceCount);
    
    // Print memory information
    std::cout << "MPI Rank: " << mpirank << ", Device Count: " << deviceCount << ", HIP Device: " << device 
              << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB" 
              << std::endl;
#endif
    
  Kokkos::initialize(argc, argv);
  {        
    if( argc >= 4 ) {
    }
    else {      
      printf("Usage: ./cppfile nummodels InputFile OutputFile\n");
      return 1;
    }                
    
    //Int nummodels, restart, mpiprocs, mpirank, shmrank, ncores, nthreads, backend;    
    string mystr = string(argv[1]);
    nummodels = stoi(mystr);  // number of pde models
    string filein[10]; 
    string fileout[10];
    
    // two-physics and two-domain problems  
    int mpiprocs0 = 0;
    if (nummodels>100) {
      mpiprocs0 = nummodels - 100;
      nummodels = 2;
    }
    
    if ((mpiprocs0 > 0) && (mpiprocs<= mpiprocs0)) {
      printf("For two-domain problem, total number of MPI processors (%d) must be greater than # MPI processors on the first domain (%d)\n", mpiprocs, mpiprocs0);
      return 1;
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
    
    if (pdemodel[0]->disc.common.AVdistfunction==1) {
      avdistfunc(pdemodel, out, nummodels, backend);
    }    
    else if ((pdemodel[0]->disc.common.tdep==1) && (pdemodel[0]->disc.common.runmode==0)) {
                        
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
    
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    ///////////////////// Start Pseudo-time step //////////////////////
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    else if ((pdemodel[0]->disc.common.tdep==1) && (pdemodel[0]->disc.common.runmode==10)) {
                        
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
        dstype monitor_diff, monitor_scale, delta_monitor;
        int NLiters = pdemodel[0]->disc.common.nonlinearSolverMaxIter;

        // time stepping with DIRK schemes
        for (Int istep=0; istep<pdemodel[0]->disc.common.tsteps; istep++)            
        {            
            for (int i=0; i<nummodels; i++) {
                pdemodel[i]->disc.common.nonlinearSolverMaxIter = 1;

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
                    // auto begin = chrono::high_resolution_clock::now();   
                    
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
                    
                    // TODO: input wprev
                    pdemodel[i]->disc.evalMonitor(pdemodel[i]->disc.tmp.tempn,  pdemodel[i]->disc.sol.udg, pdemodel[i]->disc.sol.wdg, pdemodel[i]->disc.common.nc, backend);
                    pdemodel[i]->disc.evalMonitor(pdemodel[i]->disc.tmp.tempg,  pdemodel[i]->solv.sys.udgprev, pdemodel[i]->solv.sys.wprev, pdemodel[i]->disc.common.ncu, backend);
                    
                    ArrayAXPBY(pdemodel[i]->disc.tmp.tempn, pdemodel[i]->disc.tmp.tempn, pdemodel[i]->disc.tmp.tempg, 1.0, -1.0, pdemodel[i]->disc.common.npe*pdemodel[i]->disc.common.ncm*pdemodel[i]->disc.common.ne);            
                    
                    monitor_diff  = PNORM(pdemodel[i]->disc.common.cublasHandle,  pdemodel[i]->disc.common.npe*pdemodel[i]->disc.common.ncm*pdemodel[i]->disc.common.ne, pdemodel[i]->disc.tmp.tempn, backend);
                    monitor_scale = PNORM(pdemodel[i]->disc.common.cublasHandle, pdemodel[i]->disc.common.npe*pdemodel[i]->disc.common.ncm*pdemodel[i]->disc.common.ne,  pdemodel[i]->disc.tmp.tempg, backend);

                    delta_monitor = monitor_diff / monitor_scale;
                    std::cout << "delta_monitor: " << delta_monitor << std::endl;

                    if ((delta_monitor > 1.0 || pdemodel[i]->solv.sys.alpha < 0.1))
                    {
                        std::cout << "Linesearch failed or large change in solution: reducing timestep" << std::endl;
                        
                        // Revert time step
                        pdemodel[i]->disc.common.time = pdemodel[i]->disc.common.time - pdemodel[i]->disc.common.dt[istep]*pdemodel[i]->disc.common.DIRKcoeff_t[j];

                        // Copy udg old to udg
                        ArrayCopy(pdemodel[i]->disc.sol.udg, pdemodel[i]->solv.sys.udgprev, pdemodel[i]->disc.common.npe*pdemodel[i]->disc.common.ncu*pdemodel[i]->disc.common.ne2);

                        // Compute a new UH
                        GetFaceNodes(pdemodel[i]->disc.sol.uh, pdemodel[i]->disc.sol.udg, pdemodel[i]->disc.mesh.f2e, pdemodel[i]->disc.mesh.perm, pdemodel[i]->disc.common.npf, pdemodel[i]->disc.common.ncu, pdemodel[i]->disc.common.npe, pdemodel[i]->disc.common.nc, pdemodel[i]->disc.common.nf);

                        // Recompute gradient from udg old
                        hdgGetQ(pdemodel[i]->disc.sol.udg, pdemodel[i]->disc.sol.uh, pdemodel[i]->disc.sol, pdemodel[i]->disc.res, pdemodel[i]->disc.mesh, pdemodel[i]->disc.tmp, pdemodel[i]->disc.common, backend);
    
                        // decrease timestep by 10
                        pdemodel[i]->disc.common.dt[istep+1] = pdemodel[i]->disc.common.dt[istep] / 10;
                        std::cout << "next time step: " << pdemodel[i]->disc.common.dt[istep+1] << std::endl;
                        if (pdemodel[i]->disc.common.dt[istep+1] < 1e-8){
                            std::cout << "WARNING: PTC stalled" << std::endl;
                            istep = pdemodel[0]->disc.common.tsteps+1;
                        }
                    }
                    else if (delta_monitor < 0.1)
                    {
                        if (pdemodel[i]->disc.common.linearSolverIter < 200){
                            // increase timestep by 2
                            pdemodel[i]->disc.common.dt[istep+1] = pdemodel[i]->disc.common.dt[istep]*2;
                        }
                        else{
                            std::cout << "Too many GMRES iterations, not increasing timestep" << std::endl;
                            pdemodel[i]->disc.common.dt[istep+1] = pdemodel[i]->disc.common.dt[istep]*1;
                        }
                        

                        if (delta_monitor < 1e-4 && pdemodel[i]->disc.common.dt[istep] > 1e-3)
                        {
                            std::cout << "Steady solve..." << std::endl;
                            pdemodel[i]->disc.common.tdep=0;
                            pdemodel[i]->disc.common.nonlinearSolverMaxIter = NLiters;
                            pdemodel[i]->SolveProblem(out[i], backend);
                            istep = pdemodel[0]->disc.common.tsteps+1;
                            // Solve steady and terminate
                        }    
                    }
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
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    ///////////////////// End Pseudo-time step //////////////////////
    ///////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
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
  }
  Kokkos::finalize();  
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif         
  
  return 0;             
     
}
