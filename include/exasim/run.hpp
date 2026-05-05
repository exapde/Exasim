// SPDX-License-Identifier: see LICENSE
//
// <exasim/run.hpp> — `exasim::run<M>(argc, argv)` façade.
//
// The MPI / Kokkos / argv / DIRK / multi-domain / pseudo-time loop
// from `backend/Main/main.cpp` factored into a templated function so
// that user main.cpp's collapse to:
//
//     #include <exasim/run.hpp>
//     #include "my_model.hpp"
//     int main(int argc, char** argv) {
//         return exasim::run<MyModel>(argc, argv);
//     }
//
// Calling `exasim::run<exasim::detail::AbiAdapter>(argc, argv)`
// reproduces the legacy `cput2cEXASIM` / `cpumpit2cEXASIM` /
// `cpuEXASIM` binaries' main bodies. So `backend/Main/main.cpp` is
// itself just a thin caller of `exasim::run<>` (M defaults to
// AbiAdapter), preserving every legacy code path.
//
// Side effects (MPI_Init/Finalize, Kokkos::initialize/finalize) live
// inside this function — calling it twice in one process would be a
// bug. That matches the previous main.cpp behavior.

#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <filesystem>

#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <cstring>

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

#if defined(HAVE_TEXT2CODE) || defined(HAVE_BUILTINMODEL)
#define HAVE_SHARED_MODEL_LIB
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

#ifdef HAVE_SHARED_MODEL_LIB
#include <array>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#endif

#include <Kokkos_Core.hpp>

// PR #73 review NB2: backend headers used to rely on a global
// `using namespace std;` here, leaking the entire std namespace
// into every consumer TU. That directive is gone — backend code
// now qualifies `std::*` directly.

#include <backend/Common/common.h>
#include <backend/Common/cpuimpl.h>
#include <backend/Common/kokkosimpl.h>
#include <backend/Common/pblas.h>

#include <backend/Discretization/discretization.hpp>
#include <backend/Preconditioning/preconditioner.hpp>
#include <backend/Solver/solver.hpp>
#include <backend/Visualization/visualization.hpp>
#include <backend/PointLocator/pointlocator.hpp>
#include <backend/Solution/solution.hpp>

#ifdef HAVE_SHARED_MODEL_LIB
#include <backend/Preprocessing/preprocessing.hpp>
#endif

#include "detail/abi_adapter.hpp"

namespace exasim {
namespace detail {

// RAII guards for MPI and Kokkos init/finalize. exasim::run<M>() has
// many early-return error paths after MPI_Init / Kokkos::initialize;
// without RAII, each one had to remember to call the matching
// finalize, and several didn't — leaving MPI/Kokkos half-shut-down
// on the way out. With these guards, every return point unwinds the
// stack and the destructors fire in reverse construction order, so
// Kokkos::finalize precedes MPI_Finalize automatically.
//
// Each guard is "ownership-aware": if MPI/Kokkos is already
// initialized when the guard is constructed, the guard does NOT
// call finalize on destruction. This makes run<M>() safe to call
// from a host that already manages MPI/Kokkos lifecycle (e.g. a
// long-lived facade that calls run<M>() multiple times).

#ifdef HAVE_MPI
struct MpiScope {
    bool _owned = false;
    MpiScope(int* argc, char*** argv) {
        int already = 0;
        MPI_Initialized(&already);
        if (!already) {
            MPI_Init(argc, argv);
            _owned = true;
        }
    }
    ~MpiScope() {
        if (_owned) {
            int finalized = 0;
            MPI_Finalized(&finalized);
            if (!finalized) MPI_Finalize();
        }
    }
    MpiScope(const MpiScope&) = delete;
    MpiScope& operator=(const MpiScope&) = delete;
};
#endif

struct KokkosScope {
    bool _owned = false;
    KokkosScope(int& argc, char** argv) {
        if (!Kokkos::is_initialized() && !Kokkos::is_finalized()) {
            Kokkos::initialize(argc, argv);
            _owned = true;
        }
    }
    ~KokkosScope() {
        if (_owned && Kokkos::is_initialized() && !Kokkos::is_finalized()) {
            Kokkos::finalize();
        }
    }
    KokkosScope(const KokkosScope&) = delete;
    KokkosScope& operator=(const KokkosScope&) = delete;
};

} // namespace detail

template <class M = exasim::detail::AbiAdapter>
inline int run(int argc, char** argv) {

    Int nummodels=1, mpiprocs=1, mpirank=0, shmrank=0, backend=0, localprocs=1, localrank=0;
    (void)localprocs; (void)localrank;

#ifdef HAVE_MPI
    detail::MpiScope _mpi_scope(&argc, &argv);

    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;

    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);

    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);

    MPI_Comm_rank(shmcomm, &shmrank);
#endif

#ifdef HAVE_OPENMP
    backend=1;
#else
    backend=0;
#endif
#ifdef HAVE_CUDA
    backend=2;
#endif
#ifdef HAVE_HIP
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
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) {
        if (mpirank==0) std::cerr << "No CUDA devices available" << std::endl;
        return 1;
    }
    cudaSetDevice(shmrank % deviceCount);
    cudaGetDevice( &device );
    size_t available, total;
    cudaMemGetInfo(&available, &total);

    cudaGetDeviceCount(&deviceCount);

    char busid[64] = {0};
    CHECK( cudaDeviceGetPCIBusId(busid, sizeof(busid), device) );

    std::cout << "MPI Rank: " << mpirank << ", Device Count: " << deviceCount << ", CUDA Device: " << device
              << ", PCI: " << busid << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB"
              << std::endl;
#endif

#ifdef HAVE_HIP
    int device;
    int deviceCount = 0;
    hipGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) {
        if (mpirank==0) std::cerr << "No HIP devices available" << std::endl;
        return 1;
    }
    CHECK(hipSetDevice(shmrank % deviceCount));

    CHECK(hipGetDevice(&device));

    size_t available, total;
    CHECK(hipMemGetInfo(&available, &total));

    hipGetDeviceCount(&deviceCount);

    char busid[64] = {0};
    CHECK( hipDeviceGetPCIBusId(busid, sizeof(busid), device) );

    std::cout << "MPI Rank: " << mpirank << ", Device Count: " << deviceCount << ", HIP Device: " << device
              << ", PCI: " << busid << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB"
              << std::endl;
#endif

  detail::KokkosScope _kokkos_scope(argc, argv);
  {

    std::string filein[10];
    std::string fileout[10];
    std::string exasimpath = "";
    int mpiprocs0 = 0;
    int restart = 0;
    int builtinmodelID=0;

#ifdef HAVE_SHARED_MODEL_LIB
    if (argc < 2) {
        if (mpirank==0) std::cerr << "Usage: ./Exasim <pdeapp.txt>\n";
        return 1;
    }

    if (std::filesystem::exists(argv[1])) {
        if (mpirank==0) std::cout << "Running Exasim with this text file ("<< argv[1] << ") ... \n\n";
    } else {
        error("Error: Input file does not exist.\n");
    }

    InputParams params = parseInputFile(argv[1], mpirank);
    PDE pde = initializePDE(params, mpirank);
    nummodels = 1;
    filein[0] = pde.datainpath + "/";
    fileout[0] = make_path(pde.dataoutpath, "out");
    exasimpath = pde.exasimpath;
    builtinmodelID = pde.builtinmodelID;
    if (mpirank==0) std::cout << "builtinmodelID = "<< builtinmodelID << "\n";

    if (pde.gendatain == 0) {
      CPreprocessing preproc(argv[1], mpirank, mpiprocs);
      if (mpiprocs == 1)
        preproc.SerialPreprocessing();
      else {
#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
        preproc.ParallelPreprocessing(EXASIM_COMM_LOCAL);
#endif
      }
    }

#else
    if (argc < 3) {
      printf("Usage: ./cppfile nummodels InputFile(s) OutputFile(s) [restart]\n");
      return 1;
    }

    std::string mystr = std::string(argv[1]);
    try {
        nummodels = std::stoi(mystr);
    } catch (...) {
        if (mpirank==0) std::cerr << "Invalid nummodels: " << mystr << std::endl;
        return 1;
    }

    if (nummodels>100) {
      mpiprocs0 = nummodels - 100;
      nummodels = 2;
    }

    if ((mpiprocs0 > 0) && (mpiprocs<= mpiprocs0)) {
      printf("For two-domain problem, total number of MPI processors (%d) must be greater than # MPI processors on the first domain (%d)\n", mpiprocs, mpiprocs0);
      return 1;
    }

    for (int i=0; i<nummodels; i++) {
        filein[i]  = std::string(argv[2*i+2]);
        fileout[i]  = std::string(argv[2*i+3]);
    }

    if (argc >= (2*nummodels + 3)) {
        mystr = std::string(argv[2*nummodels+2]);
        try { restart = std::stoi(mystr); }
        catch (...) {
            if (mpirank==0) std::cerr << "Invalid restart value: " << mystr << std::endl;
            return 1;
        }
    }

    std::filesystem::path cwd = std::filesystem::current_path();
    exasimpath = trimToSubstringAtLastOccurence(cwd, "Exasim");
    if (exasimpath == "") exasimpath = trimToSubstringAtLastOccurence(fileout[0], "Exasim");
#endif

    if (mpirank==0) std::cout << "exasimpath = "<<exasimpath<<std::endl;

    if (mpiprocs0 > 0) {
      nummodels = 1;
#ifdef HAVE_MPI
      int color = (mpirank < mpiprocs0) ? 0 : 1;
      MPI_Comm_split(MPI_COMM_WORLD, color, mpirank, &EXASIM_COMM_LOCAL);

      MPI_Comm_size(EXASIM_COMM_LOCAL, &localprocs);
      MPI_Comm_rank(EXASIM_COMM_LOCAL, &localrank);
#endif
    }

    Int fileoffset = 0;
    Int gpuid = 0;

    // initialize PDE models — templated on M (defaults to AbiAdapter
    // for the legacy main.cpp wrapper; user code instantiates with
    // their own Model type).
    CSolution<M>** pdemodel = new CSolution<M>*[nummodels];
    std::ofstream* out = new std::ofstream[nummodels];

    for (int i=0; i<nummodels; i++) {
        if (mpiprocs0==0) {
          pdemodel[i] = new CSolution<M>(filein[i], fileout[i], exasimpath, mpiprocs, mpirank, fileoffset, gpuid, backend, builtinmodelID);
        }
        else if (mpiprocs0 > 0) {
          if (mpirank < mpiprocs0) {
            pdemodel[i] = new CSolution<M>(filein[0], fileout[0], exasimpath, mpiprocs, mpirank, fileoffset, gpuid, backend, builtinmodelID);
          }
          else {
            fileoffset = mpiprocs0;
            std::cout<<filein[1]<<std::endl;
            std::cout<<fileout[1]<<std::endl;
            pdemodel[i] = new CSolution<M>(filein[1], fileout[1], exasimpath, mpiprocs, mpirank, fileoffset, gpuid, backend, builtinmodelID);
          }
        }

        pdemodel[i]->disc.common.nomodels = nummodels;

        pdemodel[i]->disc.common.ncarray = new Int[nummodels];
        pdemodel[i]->disc.sol.udgarray = new dstype*[nummodels];

        if (pdemodel[i]->disc.common.timestepOffset>0)
            restart = pdemodel[i]->disc.common.timestepOffset;

        if (restart>0) {
            pdemodel[i]->disc.common.timestepOffset = restart;
            pdemodel[i]->disc.common.time = restart*pdemodel[i]->disc.common.dt[0];
        }

        if (pdemodel[i]->disc.common.mpiRank==0 && pdemodel[i]->disc.common.saveResNorm==1) {
            std::string filename = pdemodel[i]->disc.common.fileout + "_residualnorms" + NumberToString(i) + ".bin";
            out[i].open(filename.c_str(), std::ios::out | std::ios::binary);
            if (!out[i]) error("Unable to open file " + filename);
        }
    }

    for (int i=0; i<nummodels; i++)
        for (int j=0; j<nummodels; j++) {
            pdemodel[i]->disc.common.ncarray[j] = pdemodel[j]->disc.common.nc;
            pdemodel[i]->disc.sol.udgarray[j] = &pdemodel[j]->disc.sol.udg[0];
        }


    if (pdemodel[0]->disc.common.AVdistfunction==1) {
      avdistfunc(pdemodel, out, nummodels, backend);
    }
    else if ((pdemodel[0]->disc.common.tdep==1) && (pdemodel[0]->disc.common.runmode==0)) {

        for (int i=0; i<nummodels; i++) {
            if (restart>0) {
                pdemodel[i]->disc.common.currentstep = -1;
                pdemodel[i]->ReadSolutions(backend);
            }
            pdemodel[i]->InitSolution(backend);
        }

        dstype time = pdemodel[0]->disc.common.time;

        for (Int istep=0; istep<pdemodel[0]->disc.common.tsteps; istep++)
        {
            for (int i=0; i<nummodels; i++) {
                pdemodel[i]->disc.common.currentstep = istep;
                PreviousSolutions(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.common, backend);
            }

            for (Int j=0; j<pdemodel[0]->disc.common.tstages; j++) {

                if (pdemodel[0]->disc.common.mpiRank==0)
                    printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",istep+1,j+1,time + pdemodel[0]->disc.common.dt[istep]*pdemodel[0]->disc.common.DIRKcoeff_t[j]);

                for (int i=0; i<nummodels; i++) {
                    pdemodel[i]->disc.common.currentstage = j;
                    pdemodel[i]->disc.common.time = time + pdemodel[i]->disc.common.dt[istep]*pdemodel[i]->disc.common.DIRKcoeff_t[j];

                    UpdateSource(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.app, pdemodel[i]->disc.res, pdemodel[i]->disc.common, backend);

                    pdemodel[i]->SteadyProblem(out[i], backend);

                    UpdateSolution<M>(pdemodel[i]->disc.sol, pdemodel[i]->solv.sys, pdemodel[i]->disc.app, pdemodel[i]->disc.res, pdemodel[i]->disc.tmp, pdemodel[i]->disc.common, backend);
                }
            }

            for (int i=0; i<nummodels; i++) {
                if (pdemodel[i]->disc.common.compudgavg == 1) {
                    ArrayAXPBY(pdemodel[i]->disc.sol.udgavg, pdemodel[i]->disc.sol.udgavg, pdemodel[i]->disc.sol.udg, one, one, pdemodel[i]->disc.common.ndofudg1);
                    ArrayAddScalar(&pdemodel[i]->disc.sol.udgavg[pdemodel[i]->disc.common.ndofudg1], one, 1);
                }

                pdemodel[i]->disc.computeAverageSolutionsOnBoundary();

                pdemodel[i]->SaveSolutions(backend);
                pdemodel[i]->SaveQoI(backend);
                if (pdemodel[i]->vis.savemode > 0) pdemodel[i]->SaveParaview(backend);
                pdemodel[i]->SaveSolutionsOnBoundary(backend);
                if (pdemodel[i]->disc.common.nce>0)
                    pdemodel[i]->SaveOutputCG(backend);
            }

            time = time + pdemodel[0]->disc.common.dt[istep];
        }
    }
    // Pseudo-time stepping (single model only).
    else if ((pdemodel[0]->disc.common.tdep==1) && (pdemodel[0]->disc.common.runmode==10 || pdemodel[0]->disc.common.runmode==11)) {
        if (restart>0) {
            pdemodel[0]->disc.common.currentstep = -1;
            pdemodel[0]->ReadSolutions(backend);
        }
        pdemodel[0]->InitSolution(backend);
        pdemodel[0]->SteadyProblem_PTC(out[0], backend);
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
                pdemodel[i]->SaveQoI(backend);
                if (pdemodel[i]->vis.savemode > 0) pdemodel[i]->SaveParaview(backend);
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
                std::string filename = pdemodel[i]->disc.common.fileout + "_np" + NumberToString(pdemodel[i]->disc.common.mpiRank) + ".bin";
                writearray2file(filename, pdemodel[i]->disc.sol.udg, pdemodel[i]->disc.common.ndofudg1, backend);
            }
        }
    }

    for (int i=0; i<nummodels; i++) {
        if (pdemodel[i]->disc.common.mpiRank==0 && pdemodel[i]->disc.common.saveResNorm==1)
            out[i].close();
    }

    for (int i = 0; i < nummodels; i++) {
        delete[] pdemodel[i]->disc.common.ncarray;
        delete[] pdemodel[i]->disc.sol.udgarray;
        delete pdemodel[i];
    }

    delete[] pdemodel;
    delete[] out;
  }
  // _kokkos_scope and (if HAVE_MPI) _mpi_scope finalize here via RAII
  // when the function returns; explicit Kokkos::finalize / MPI_Finalize
  // are gone so every early-return path is balanced for free.
  return 0;
}

} // namespace exasim
