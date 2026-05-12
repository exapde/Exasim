#pragma once

#include <fstream>
#include <memory>
#include <string>
#include <vector>

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

#ifdef _BUILTINMODEL
#define HAVE_BUILTINMODEL
#endif

#ifdef _KOKKOSKERNEL
#define HAVE_KOKKOSKERNEL
#endif

#if defined(HAVE_TEXT2CODE) || defined(HAVE_BUILTINMODEL) || defined(HAVE_KOKKOSKERNEL)
#define HAVE_SHARED_MODEL_LIB
#endif

#ifdef _UDERDEFINED
#define HAVE_UDERDEFINED
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

#include "driver_abi.hpp"
#include "modeldefaults.hpp"

#if defined(HAVE_MPI) || defined(_MPI)
#include <mpi.h>
extern MPI_Comm EXASIM_COMM_WORLD;
extern MPI_Comm EXASIM_COMM_LOCAL;
#else
using MPI_Comm = int;
static constexpr MPI_Comm MPI_COMM_NULL = 0;
#endif

// Forward declaration keeps this public header light. The implementation file
// should include the full CSolution definition before defining the methods.
class CSolution;

struct ExasimPoint {
    dstype x = 0.0;
    dstype y = 0.0;
    dstype z = 0.0;
};

class ExasimSolver
{
public:
    ExasimSolver();
    ~ExasimSolver();

    ExasimSolver(const ExasimSolver&) = delete;
    ExasimSolver& operator=(const ExasimSolver&) = delete;
    ExasimSolver(ExasimSolver&&) = delete;
    ExasimSolver& operator=(ExasimSolver&&) = delete;

    int InitializeEnvironment(int argc, char** argv, MPI_Comm comm);
    int Initialize(int argc, char** argv, MPI_Comm comm);    
    int ParseInputs(int argc, char** argv);
    int SetModelDefinition(const int modelnumber, const int builtinmodelID,
                           const ExasimDriverABI& abi);
    int InitializeModels();
    int BuildModels();
    int CoupleModels();
    int OpenOutputStreams();
    void CloseOutputStreams();
    void DestroyModels();
    void FinalizeEnvironment();    
    int Finalize();

    int IntializeMeshInterface(const int modelnumber,
                               const int ncuext, const int ncuint,
                               const int ibc, const int comperm,
                               const int offset, MPI_Comm comm);
    int SaveState(const int modelnumber);
    int RestoreState(const int modelnumber);
    int ClearSavedState(const int modelnumber);
    std::vector<ExasimPoint> getInterfacePoints() const;
    void getInterfaceFluxes(std::vector<dstype>& send_flux) const;
    void setInterfaceFluxes(const std::vector<dstype>& recv_flux);

    int InitializeSolution();
    int InitializeSolution(const int modelnumber);
    int RunAVDistanceFunction();
    int RunTimeDependent();
    int RunTimeDependent(const int modelnumber);
    int RunTimeDependent(int start, int steps);
    int RunTimeDependent(const int modelnumber, int start, int steps);
    int RunPseudoTime();
    int RunPseudoTime(const int modelnumber);
    int RunSteady();
    int RunSteady(const int modelnumber);
    int RunSolveProblemOrPostprocess();
    int Solve();
    int Solve(const int modelnumber);

    bool IsInitialized() const;
    int NumModels() const;
    int NumModelDefinitions() const;
    int BuiltinModelID(const int modelnumber) const;
    int Backend() const;
    CSolution* Model(int i);

    int _ncuext; // number of components of external fields
    int _ncuint; // number of components of internal fields
    int _ibc; // boundary condition index for the interface
    int _comperm; // interface flux map
    int _offset; // component offset in odg
    int _rank = 0;
    int _size = 1;
    MPI_Comm _comm = MPI_COMM_NULL;
    int *faces= nullptr;        // faces on the interface
    int ncx;          // number of compoments of (xdg)
    int nfaces;       // number of faces on the interface
    int npf;          // number of nodes on master face           
    int ngf;          // number of gauss poInts on master face
    dstype *xdgint= nullptr;   // XDG nodes on the interface
    dstype *nlint= nullptr;    // normal vectors on the interface
    dstype *xdggint= nullptr;  // XDG gauss nodes on the interface
    dstype *nlgint= nullptr;   // normal vectors at Gauss nodes on the interface

private:
    std::vector<std::unique_ptr<CSolution>> models_;
    std::vector<ExasimDriverABI> model_abis_;
    std::vector<std::ofstream> residual_outputs_;
    std::vector<std::string> filein_;
    std::vector<std::string> fileout_;
    std::string exasimpath_;

    MPI_Comm world_comm_ = MPI_COMM_NULL;
    MPI_Comm local_comm_ = MPI_COMM_NULL;

    int nummodels_ = 0;
    int backend_ = 0;
    int restart_ = 0;
    std::vector<int> builtinmodelID_;
    int interface_modelnumber_ = -1;

    int mpiprocs_ = 1;
    int mpirank_ = 0;
    int shmrank_ = 0;
    int localprocs_ = 1;
    int localrank_ = 0;
    int mpiprocs0_ = 0;

    bool initialized_ = false;
    bool owns_mpi_ = false;
    bool owns_kokkos_ = false;    
};
