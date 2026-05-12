#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <filesystem>

// Use C++ header for C string utilities
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>

#ifdef _TEXT2CODE
#define HAVE_TEXT2CODE
#endif

#if defined(HAVE_TEXT2CODE) || defined(HAVE_BUILTINMODEL)
#define HAVE_SHARED_MODEL_LIB
#endif

#ifdef _MPI
#define HAVE_MPI
#endif

#ifdef HAVE_MPI
#include <mpi.h>
MPI_Comm EXASIM_COMM_WORLD = MPI_COMM_NULL;
MPI_Comm EXASIM_COMM_LOCAL = MPI_COMM_NULL;
#endif

#include <Kokkos_Core.hpp>

// g++ -O2 -std=c++17 -I../../kokkos/buildserial/include postprocess.cpp ../../kokkos/buildserial/lib/libkokkoscore.a -o postprocess -lblas -llapack
// mpicxx -O2 -std=c++17 -D _MPI -I../../kokkos/buildserial/include postprocess.cpp ../../kokkos/buildserial/lib/libkokkoscore.a -o postprocess -lblas -llapack

// g++ -O2 -std=c++17 -I/Users/cuongnguyen/Documents/GitHub/Exasim/kokkos/buildserial/include postprocess.cpp /Users/cuongnguyen/Documents/GitHub/Exasim/kokkos/buildserial/lib/libkokkoscore.a -o postprocess -lblas -llapack
// mpicxx -O2 -std=c++17 -D _MPI -I/Users/cuongnguyen/Documents/GitHub/Exasim/kokkos/buildserial/include postprocess.cpp /Users/cuongnguyen/Documents/GitHub/Exasim/kokkos/buildserial/lib/libkokkoscore.a -o postprocess -lblas -llapack

using namespace std;

#include "../Common/common.h"        // declaration of variables and structs
#include "../Common/cpuimpl.h"       // interface to cpu functions
#include "../Common/kokkosimpl.h"    // interface to kokkos functions
#include "../Common/pblas.h"         // wrappers for blas libaries and MPI     

#include "../Discretization/postdiscretization.cpp" // discretization class
#include "../Preconditioning/postpreconditioner.cpp" // preconditioner class
#include "../Solver/postsolver.cpp"                 // solver class
#include "../Visualization/visualization.cpp"  //  visualization class
#include "../Solution/postsolution.cpp"             // solution class

int main(int argc, char** argv) 
{   
  
    Int nummodels=1, mpiprocs=1, mpirank=0, shmrank=0, backend=0, localprocs=1, localrank=0;    
    
#ifdef HAVE_MPI    
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // set EXASIM_COMM_WORLD and EXASIM_COMM_LOCAL
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;

    // Get the number of processes    
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);

    // Get the rank of the process
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
   
    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    
    MPI_Comm_rank(shmcomm, &shmrank);
#endif                
                      
  Kokkos::initialize(argc, argv);
  {        

    string filein[10]; 
    string fileout[10];
    string exasimpath = "";  
    int mpiprocs0 = 0;
    int restart = 0;
    int postmode = 0;
    int nsca = 0;
    int nvec = 0;
    int nten = 0;
    int nsurf = 0;
    int nvqoi = 0;
    int builtinmodelID=0;
    
    if (argc < 3) {
      printf("Usage: ./postprocess nummodels InputFile(s) OutputFile(s) [timestep] [postmode]\n");
      return 1;
    }
    
    string mystr = string(argv[1]);
    try {
        nummodels = stoi(mystr);  // number of pde models
    } catch (...) {
        if (mpirank==0) std::cerr << "Invalid nummodels: " << mystr << std::endl;
        return 1;
    }        
      
    // two-physics and two-domain problems      
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
    }
        
    if (argc >= (2*nummodels + 3)) {
        mystr = string(argv[2*nummodels+2]);
        try { restart = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid timestep value: " << mystr << std::endl;
            return 1;
        }
    }

    if (argc >= (2*nummodels + 4)) {
        mystr = string(argv[2*nummodels+3]);
        try { postmode = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid postmode value: " << mystr << std::endl;
            return 1;
        }
    }    
    if (argc >= (2*nummodels + 5)) {
        mystr = string(argv[2*nummodels+4]);
        try { nsca = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid nsca value: " << mystr << std::endl;
            return 1;
        }
    }
    if (argc >= (2*nummodels + 6)) {
        mystr = string(argv[2*nummodels+5]);
        try { nvec = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid nvec value: " << mystr << std::endl;
            return 1;
        }
    }
    if (argc >= (2*nummodels + 7)) {
        mystr = string(argv[2*nummodels+6]);
        try { nten = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid nten value: " << mystr << std::endl;
            return 1;
        }
    }
    if (argc >= (2*nummodels + 8)) {
        mystr = string(argv[2*nummodels+7]);
        try { nsurf = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid nsurf value: " << mystr << std::endl;
            return 1;
        }
    }
    if (argc >= (2*nummodels + 9)) {
        mystr = string(argv[2*nummodels+8]);
        try { nvqoi = stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid nvqoi value: " << mystr << std::endl;
            return 1;
        }
    }
    
    std::filesystem::path cwd = std::filesystem::current_path();
    exasimpath = trimToSubstringAtLastOccurence(cwd, "Exasim");            
    if (exasimpath == "") exasimpath = trimToSubstringAtLastOccurence(fileout[0], "Exasim");                     
    if (mpirank==0) std::cout << "exasimpath = "<<exasimpath<<std::endl;
      
    // reset nummodels
    if (mpiprocs0 > 0) {
      nummodels = 1;
#ifdef HAVE_MPI          
      int color = (mpirank < mpiprocs0) ? 0 : 1;
      // split the communicator based on the color
      MPI_Comm_split(MPI_COMM_WORLD, color, mpirank, &EXASIM_COMM_LOCAL);

      // Get the number of processes    
      MPI_Comm_size(EXASIM_COMM_LOCAL, &localprocs);

      // Get the rank of the process
      MPI_Comm_rank(EXASIM_COMM_LOCAL, &localrank);   
#endif      
    }
                        
    Int fileoffset = 0;
    Int gpuid = 0;            
      
    // initialize PDE models
    CSolution** pdemodel = new CSolution*[nummodels];     
            
    for (int i=0; i<nummodels; i++) {        
        if (mpiprocs0==0) {
          pdemodel[i] = new CSolution(filein[i], fileout[i], exasimpath, 
                        mpiprocs, mpirank, fileoffset, gpuid, backend, 
                        builtinmodelID, nsca, nvec, nten, nsurf, nvqoi);                 
        }
        else if (mpiprocs0 > 0) {
          if (mpirank < mpiprocs0) { 
            pdemodel[i] = new CSolution(filein[0], fileout[0], exasimpath, 
                              mpiprocs, mpirank, fileoffset, gpuid, backend, 
                              builtinmodelID, nsca, nvec, nten, nsurf, nvqoi);       
          }
          else {
            fileoffset = mpiprocs0;
            cout<<filein[1]<<endl;
            cout<<fileout[1]<<endl;
            pdemodel[i] = new CSolution(filein[1], fileout[1], exasimpath, 
                              mpiprocs, mpirank, fileoffset, gpuid, backend, 
                              builtinmodelID, nsca, nvec, nten, nsurf, nvqoi);       
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
    }            
               
    for (int i=0; i<nummodels; i++) {                                
        if (postmode == 0) {
            pdemodel[i]->disc.common.currentstep = -1;
            pdemodel[i]->ReadSolutions(backend);   
            pdemodel[i]->SaveQoI(backend);
            if (mpirank==0) cout<<"save paraview = "<<pdemodel[i]->vis.savemode<<endl;
            if (pdemodel[i]->vis.savemode > 0) pdemodel[i]->SaveParaview(backend); 
            pdemodel[i]->SaveOutputCG(backend);            
        }
    }
        
    // After using pdemodel, delete each CSolution object and then the array of pointers
    for (int i = 0; i < nummodels; i++) {
        delete[] pdemodel[i]->disc.common.ncarray; // Delete the ncarray
        delete[] pdemodel[i]->disc.sol.udgarray; // Delete the udgarray
        delete pdemodel[i]; // Delete each CSolution object
    }

    delete[] pdemodel; // Delete the array of pointers    
  }
  Kokkos::finalize();  
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif         
  
  return 0;                  
}
