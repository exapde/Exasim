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

#ifdef _ENZYME
#define HAVE_ENZYME
#endif

#ifdef _MUTATIONPP
#define HAVE_MPP
#endif

#ifdef _MPI
#define HAVE_MPI
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
MPI_Comm EXASIM_COMM_WORLD = MPI_COMM_NULL;
MPI_Comm EXASIM_COMM_LOCAL = MPI_COMM_NULL;
#endif

#ifdef HAVE_SHARED_MODEL_LIB
#include <array>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#endif

#ifdef TIMING
#include <chrono>
#endif

#include <Kokkos_Core.hpp>

using namespace std;

#include "../Common/common.h"
#include "../Common/cpuimpl.h"
#include "../Common/kokkosimpl.h"
#include "../Common/pblas.h"

#include "../Model/ModelDispatch/model_drivers_abi.cpp"
#include "../Discretization/discretization.cpp"
#include "../Preconditioning/preconditioner.cpp"
#include "../Solver/solver.cpp"
#include "../Visualization/visualization.cpp"
#include "../PointLocator/pointlocator.cpp"
#include "../Solution/solution.cpp"

#ifdef HAVE_SHARED_MODEL_LIB
#include "../Preprocessing/preprocessing.cpp"
#endif

#include "ExasimSolver.hpp"

namespace {

static vector<CSolution*> ModelPointers(vector<unique_ptr<CSolution>>& models)
{
    vector<CSolution*> ptrs(models.size(), nullptr);
    for (size_t i = 0; i < models.size(); ++i)
        ptrs[i] = models[i].get();
    return ptrs;
}

static bool IsValidModelABI(const ExasimDriverABI& abi)
{
    return abi.abi_version == kExasimDriverABIVersion &&
           abi.struct_size == sizeof(ExasimDriverABI) &&
           abi.KokkosFlux &&
           abi.KokkosSource &&
           abi.KokkosSourcew &&
           abi.KokkosTdfunc &&
           abi.KokkosAvfield &&
           abi.KokkosEoS &&
           abi.KokkosEoSdu &&
           abi.KokkosEoSdw &&
           abi.KokkosFbou &&
           abi.KokkosUbou &&
           abi.KokkosFhat &&
           abi.KokkosUhat &&
           abi.KokkosStab &&
           abi.KokkosOutput &&
           abi.KokkosMonitor &&
           abi.KokkosVisScalars &&
           abi.KokkosVisVectors &&
           abi.KokkosVisTensors &&
           abi.KokkosQoIvolume &&
           abi.KokkosQoIboundary &&
           abi.KokkosInitu &&
           abi.KokkosInitq &&
           abi.KokkosInitudg &&
           abi.KokkosInitwdg &&
           abi.KokkosInitodg &&
           abi.cpuInitu &&
           abi.cpuInitq &&
           abi.cpuInitudg &&
           abi.cpuInitwdg &&
           abi.cpuInitodg &&
           abi.HdgFlux &&
           abi.HdgSource &&
           abi.HdgSourcew &&
           abi.HdgSourcewonly &&
           abi.HdgEoS &&
           abi.HdgFbou &&
           abi.HdgFbouonly &&
           abi.HdgFint &&
           abi.HdgFintonly &&
           abi.HdgFext &&
           abi.HdgFextonly;
}

} // namespace

ExasimSolver::ExasimSolver() = default;

ExasimSolver::~ExasimSolver()
{
    Finalize();
}

int ExasimSolver::Initialize(int argc, char** argv, MPI_Comm comm)
{
    if (initialized_)
        return 0;

    int err = InitializeEnvironment(argc, argv, comm);
    if (err) return err;

    err = ParseInputs(argc, argv);
    if (err) {
        FinalizeEnvironment();
        return err;
    }

    return InitializeModels();
}

int ExasimSolver::SetModelDefinition(const int modelnumber,
                                     const int builtinmodelID,
                                     const ExasimDriverABI& abi)
{
    if (modelnumber < 0)
        return 1;

    if (initialized_ || !models_.empty())
        return 1;

    if ((nummodels_ > 0) && (modelnumber >= NumModelDefinitions()))
        return 1;

    const size_t index = static_cast<size_t>(modelnumber);
    if (builtinmodelID_.size() <= index)
        builtinmodelID_.resize(index + 1, 0);
    if (model_abis_.size() <= index)
        model_abis_.resize(index + 1);

    builtinmodelID_[index] = builtinmodelID;
    model_abis_[index] = abi;
    return 0;
}

int ExasimSolver::InitializeModels()
{
    if (initialized_)
        return 0;

    int err = BuildModels();
    if (err) {
        Finalize();
        return err;
    }

    err = CoupleModels();
    if (err) {
        Finalize();
        return err;
    }

    err = OpenOutputStreams();
    if (err) {
        Finalize();
        return err;
    }

    initialized_ = true;
    return 0;
}

int ExasimSolver::Solve()
{
    if (!initialized_ || models_.empty())
        return 1;

    int status = 0;
    commonstruct& common = models_[0]->disc.common;

    if (common.AVdistfunction == 1) {
        status = RunAVDistanceFunction();
    }
    else if ((common.tdep == 1) && (common.runmode == 0)) {
        status = RunTimeDependent();
    }
    else if ((common.tdep == 1) && (common.runmode == 10 || common.runmode == 11)) {
        status = RunPseudoTime();
    }
    else if ((common.tdep == 0) && (common.runmode == 0)) {
        status = RunSteady();
    }
    else {
        status = RunSolveProblemOrPostprocess();
    }

    return status;
}

int ExasimSolver::Solve(const int i)
{
    if (!initialized_ || models_.empty())
        return 1;

    int status = 0;
    commonstruct& common = models_[i]->disc.common;

    if ((common.tdep == 1) && (common.runmode == 0)) {
        status = RunTimeDependent(i);
    }
    else if ((common.tdep == 1) && (common.runmode == 10 || common.runmode == 11)) {
        status = RunPseudoTime(i);
    }
    else if ((common.tdep == 0) && (common.runmode == 0)) {
        status = RunSteady(i);
    }
    else {
        status = -1;
    }

    return status;
}

int ExasimSolver::Finalize()
{
    CloseOutputStreams();
    DestroyModels();
    FinalizeEnvironment();
    initialized_ = false;
    return 0;
}

bool ExasimSolver::IsInitialized() const
{
    return initialized_;
}

int ExasimSolver::NumModels() const
{
    return nummodels_;
}

int ExasimSolver::NumModelDefinitions() const
{
    return (mpiprocs0_ > 0) ? 2 : nummodels_;
}

int ExasimSolver::BuiltinModelID(const int modelnumber) const
{
    if (modelnumber < 0)
        return 0;

    const size_t index = static_cast<size_t>(modelnumber);
    if (index >= builtinmodelID_.size())
        return 0;

    return builtinmodelID_[index];
}

int ExasimSolver::Backend() const
{
    return backend_;
}

CSolution* ExasimSolver::Model(int i)
{
    if (i < 0 || i >= static_cast<int>(models_.size()))
        return nullptr;
    return models_[i].get();
}

int ExasimSolver::InitializeEnvironment(int argc, char** argv, MPI_Comm comm)
{
    (void) argc;
    (void) argv;

#ifdef HAVE_MPI
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if (!mpiInitialized) {
        MPI_Init(&argc, &argv);
        owns_mpi_ = true;
    }

    world_comm_ = (comm == MPI_COMM_NULL) ? MPI_COMM_WORLD : comm;
    local_comm_ = world_comm_;
    EXASIM_COMM_WORLD = world_comm_;
    EXASIM_COMM_LOCAL = local_comm_;

    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs_);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank_);

    MPI_Comm shmcomm = MPI_COMM_NULL;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank_);
    MPI_Comm_free(&shmcomm);
#else
    (void) comm;
    mpiprocs_ = 1;
    mpirank_ = 0;
    shmrank_ = 0;
#endif

#ifdef HAVE_OPENMP
    backend_ = 1;
#else
    backend_ = 0;
#endif
#ifdef HAVE_CUDA
    backend_ = 2;
#endif
#ifdef HAVE_HIP
    backend_ = 3;
#endif

    if (backend_ == 2) {
        if (mpirank_ == 0)
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on CUDA platform...\n",
                   mpiprocs_, mpirank_, shmrank_);
    }
    else if (backend_ == 3) {
        if (mpirank_ == 0)
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on HIP platform...\n",
                   mpiprocs_, mpirank_, shmrank_);
    }
    else {
        if (mpirank_ == 0)
            printf("Using %d processors with %d rank and %d shmrank to solve the problem on CPU platform...\n",
                   mpiprocs_, mpirank_, shmrank_);
    }

#ifdef HAVE_CUDA
    int device;
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) {
        if (mpirank_ == 0) std::cerr << "No CUDA devices available" << std::endl;
        return 1;
    }
    cudaSetDevice(shmrank_ % deviceCount);
    cudaGetDevice(&device);

    size_t available, total;
    cudaMemGetInfo(&available, &total);
    cudaGetDeviceCount(&deviceCount);

    char busid[64] = {0};
    CHECK(cudaDeviceGetPCIBusId(busid, sizeof(busid), device));

    std::cout << "MPI Rank: " << mpirank_ << ", Device Count: " << deviceCount
              << ", CUDA Device: " << device << ", PCI: " << busid
              << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB"
              << std::endl;
#endif

#ifdef HAVE_HIP
    int device;
    int deviceCount = 0;
    hipGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) {
        if (mpirank_ == 0) std::cerr << "No HIP devices available" << std::endl;
        return 1;
    }
    CHECK(hipSetDevice(shmrank_ % deviceCount));
    CHECK(hipGetDevice(&device));

    size_t available, total;
    CHECK(hipMemGetInfo(&available, &total));
    hipGetDeviceCount(&deviceCount);

    char busid[64] = {0};
    CHECK(hipDeviceGetPCIBusId(busid, sizeof(busid), device));

    std::cout << "MPI Rank: " << mpirank_ << ", Device Count: " << deviceCount
              << ", HIP Device: " << device << ", PCI: " << busid
              << ", Available Memory: " << available / (1024.0 * 1024.0) << " MB"
              << ", Total Memory: " << total / (1024.0 * 1024.0) << " MB"
              << std::endl;
#endif

    if (!Kokkos::is_initialized()) {
        Kokkos::initialize(argc, argv);
        owns_kokkos_ = true;
    }

    return 0;
}

int ExasimSolver::ParseInputs(int argc, char** argv)
{
    filein_.clear();
    fileout_.clear();
    exasimpath_.clear();
    nummodels_ = 1;
    mpiprocs0_ = 0;
    restart_ = 0;
    const bool preserveModelDefinitions =
        !builtinmodelID_.empty() || !model_abis_.empty();

#ifdef HAVE_SHARED_MODEL_LIB
    if (argc < 2) {
        if (mpirank_ == 0) std::cerr << "Usage: ./Exasim <pdeapp.txt>\n";
        return 1;
    }

    if (std::filesystem::exists(argv[1])) {
        if (mpirank_ == 0)
            std::cout << "Running Exasim with this text file (" << argv[1] << ") ... \n\n";
    } else {
        error("Error: Input file does not exist.\n");
    }

    InputParams params = parseInputFile(argv[1], mpirank_);
    PDE pde = initializePDE(params, mpirank_);
    nummodels_ = 1;
    filein_.push_back(pde.datainpath + "/");
    fileout_.push_back(make_path(pde.dataoutpath, "out"));
    exasimpath_ = pde.exasimpath;
    if (!preserveModelDefinitions)
        builtinmodelID_.assign(1, pde.builtinmodelID);
    if (mpirank_ == 0)
        std::cout << "builtinmodelID = " << BuiltinModelID(0) << "\n";

    if (pde.gendatain == 0) {
        CPreprocessing preproc(argv[1], mpirank_, mpiprocs_);
        if (mpiprocs_ == 1)
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

    string mystr = string(argv[1]);
    try {
        nummodels_ = stoi(mystr);
    } catch (...) {
        if (mpirank_ == 0) std::cerr << "Invalid nummodels: " << mystr << std::endl;
        return 1;
    }

    if (nummodels_ > 100) {
        mpiprocs0_ = nummodels_ - 100;
        nummodels_ = 2;
    }

    if ((mpiprocs0_ > 0) && (mpiprocs_ <= mpiprocs0_)) {
        printf("For two-domain problem, total number of MPI processors (%d) must be greater than # MPI processors on the first domain (%d)\n",
               mpiprocs_, mpiprocs0_);
        return 1;
    }

    for (int i = 0; i < nummodels_; i++) {
        filein_.push_back(string(argv[2 * i + 2]));
        fileout_.push_back(string(argv[2 * i + 3]));
    }

    if (argc >= (2 * nummodels_ + 3)) {
        mystr = string(argv[2 * nummodels_ + 2]);
        try {
            restart_ = stoi(mystr);
        } catch (...) {
            if (mpirank_ == 0) std::cerr << "Invalid restart value: " << mystr << std::endl;
            return 1;
        }
    }

    std::filesystem::path cwd = std::filesystem::current_path();
    exasimpath_ = trimToSubstringAtLastOccurence(cwd, "Exasim");
    if (exasimpath_ == "")
        exasimpath_ = trimToSubstringAtLastOccurence(fileout_[0], "Exasim");
#endif

    if (mpirank_ == 0) std::cout << "exasimpath = " << exasimpath_ << std::endl;

    if (mpiprocs0_ > 0) {
        nummodels_ = 1;
#ifdef HAVE_MPI
        int color = (mpirank_ < mpiprocs0_) ? 0 : 1;
        MPI_Comm_split(EXASIM_COMM_WORLD, color, mpirank_, &local_comm_);
        EXASIM_COMM_LOCAL = local_comm_;
        MPI_Comm_size(EXASIM_COMM_LOCAL, &localprocs_);
        MPI_Comm_rank(EXASIM_COMM_LOCAL, &localrank_);
#endif
    }

    if (!preserveModelDefinitions && builtinmodelID_.empty())
        builtinmodelID_.assign(NumModelDefinitions(), 0);

    if (!builtinmodelID_.empty() &&
        builtinmodelID_.size() != static_cast<size_t>(NumModelDefinitions())) {
        if (mpirank_ == 0)
            std::cerr << "Number of builtin model IDs does not match the number of model definitions." << std::endl;
        return 1;
    }

    return 0;
}

int ExasimSolver::BuildModels()
{
    const int numModelDefinitions = NumModelDefinitions();

    if (builtinmodelID_.size() != static_cast<size_t>(numModelDefinitions)) {
        if (mpirank_ == 0)
            std::cerr << "ExasimSolver requires one builtin model ID per model definition." << std::endl;
        return 1;
    }

    if (model_abis_.size() != static_cast<size_t>(numModelDefinitions)) {
        if (mpirank_ == 0)
            std::cerr << "ExasimSolver requires one ExasimDriverABI per model definition." << std::endl;
        return 1;
    }

    for (int i = 0; i < numModelDefinitions; i++) {
        if (!IsValidModelABI(model_abis_[i])) {
            if (mpirank_ == 0)
                std::cerr << "Invalid ExasimDriverABI for model definition " << i << "." << std::endl;
            return 1;
        }
    }

    models_.clear();
    residual_outputs_.clear();
    models_.reserve(nummodels_);
    residual_outputs_.resize(nummodels_);

    Int fileoffset = 0;
    Int gpuid = 0;

    for (int i = 0; i < nummodels_; i++) {
        const int modelDefinition = (mpiprocs0_ > 0 && mpirank_ >= mpiprocs0_) ? 1 : i;

        if (mpiprocs0_ == 0) {
            models_.push_back(std::make_unique<CSolution>(
                filein_[i], fileout_[i], exasimpath_, mpiprocs_, mpirank_,
                fileoffset, gpuid, backend_, builtinmodelID_[modelDefinition],
                model_abis_[modelDefinition]));
        }
        else if (mpiprocs0_ > 0) {
            if (mpirank_ < mpiprocs0_) {
                models_.push_back(std::make_unique<CSolution>(
                    filein_[0], fileout_[0], exasimpath_, mpiprocs_, mpirank_,
                    fileoffset, gpuid, backend_, builtinmodelID_[modelDefinition],
                    model_abis_[modelDefinition]));
            }
            else {
                fileoffset = mpiprocs0_;
                cout << filein_[1] << endl;
                cout << fileout_[1] << endl;
                models_.push_back(std::make_unique<CSolution>(
                    filein_[1], fileout_[1], exasimpath_, mpiprocs_, mpirank_,
                    fileoffset, gpuid, backend_, builtinmodelID_[modelDefinition],
                    model_abis_[modelDefinition]));
            }
        }

        models_[i]->disc.common.nomodels = nummodels_;
        models_[i]->disc.common.ncarray = new Int[nummodels_];
        models_[i]->disc.sol.udgarray = new dstype*[nummodels_];

        if (models_[i]->disc.common.timestepOffset > 0)
            restart_ = models_[i]->disc.common.timestepOffset;

        if (restart_ > 0) {
            models_[i]->disc.common.timestepOffset = restart_;
            models_[i]->disc.common.time = restart_ * models_[i]->disc.common.dt[0];
        }
    }

    return 0;
}

int ExasimSolver::CoupleModels()
{
    for (int i = 0; i < nummodels_; i++) {
        for (int j = 0; j < nummodels_; j++) {
            models_[i]->disc.common.ncarray[j] = models_[j]->disc.common.nc;
            models_[i]->disc.sol.udgarray[j] = &models_[j]->disc.sol.udg[0];
        }
    }

    return 0;
}

int ExasimSolver::OpenOutputStreams()
{
    for (int i = 0; i < nummodels_; i++) {
        if (models_[i]->disc.common.mpiRank == 0 &&
            models_[i]->disc.common.saveResNorm == 1) {
            string filename = models_[i]->disc.common.fileout +
                              "_residualnorms" + NumberToString(i) + ".bin";
            residual_outputs_[i].open(filename.c_str(), ios::out | ios::binary);
            if (!residual_outputs_[i])
                error("Unable to open file " + filename);
        }
    }

    return 0;
}

int ExasimSolver::IntializeMeshInterface(const int modelnumber,
                                         const int ncuext, const int ncuint,
                                         const int ibc, const int comperm,
                                         const int offset, MPI_Comm comm)
{
    if (modelnumber < 0 ||
        modelnumber >= static_cast<int>(models_.size()) ||
        !models_[modelnumber])
        return 1;

    CSolution& model = *models_[modelnumber];
    const int interfaceBackend = model.disc.common.backend;

    if (faces) {
        TemplateFree(faces, interfaceBackend);
        faces = nullptr;
    }
    if (xdgint) {
        TemplateFree(xdgint, interfaceBackend);
        xdgint = nullptr;
    }
    if (nlint) {
        TemplateFree(nlint, interfaceBackend);
        nlint = nullptr;
    }
    if (xdggint) {
        TemplateFree(xdggint, interfaceBackend);
        xdggint = nullptr;
    }
    if (nlgint) {
        TemplateFree(nlgint, interfaceBackend);
        nlgint = nullptr;
    }
    if (model.disc.sol.uext) {
        TemplateFree(model.disc.sol.uext, interfaceBackend);
        model.disc.sol.uext = nullptr;
        model.disc.sol.szuext = 0;
    }

    _ncuext = ncuext;
    _ncuint = ncuint;
    _ibc = ibc;
    _comperm = comperm;
    _offset = offset;
    _comm = comm;
    interface_modelnumber_ = modelnumber;

#ifdef HAVE_MPI
    MPI_Comm_rank(_comm, &_rank);
    MPI_Comm_size(_comm, &_size);
#else
    (void) comm;
    _rank = 0;
    _size = 1;
#endif

    backend_ = interfaceBackend;
    ncx = model.disc.common.ncx;
    npf = model.disc.common.npf;
    ngf = model.disc.common.ngf;
    nfaces = model.disc.getFacesOnInterface(&faces, ibc + 1);

    model.disc.common.ncuext = ncuext;
    TemplateMalloc(&model.disc.sol.uext, ngf * nfaces * ncuext,
                   interfaceBackend);
    model.disc.sol.szuext = ngf * nfaces * ncuext;

    TemplateMalloc(&xdgint, npf * nfaces * model.disc.common.ncx,
                   interfaceBackend);
    TemplateMalloc(&nlint, npf * nfaces * model.disc.common.nd,
                   interfaceBackend);
    TemplateMalloc(&xdggint, ngf * nfaces * model.disc.common.ncx,
                   interfaceBackend);
    TemplateMalloc(&nlgint, ngf * nfaces * model.disc.common.nd,
                   interfaceBackend);

    model.disc.getDGNodesOnInterface(xdgint, faces, nfaces);
    model.disc.getNormalVectorOnInterface(nlint, xdgint, nfaces);
    model.disc.getFieldsAtGaussPointsOnInterface(xdggint, xdgint, nfaces,
                                                 model.disc.common.ncx);
    model.disc.getFieldsAtGaussPointsOnInterface(nlgint, nlint, nfaces,
                                                 model.disc.common.nd);

    return 0;
}

int ExasimSolver::SaveState(const int modelnumber)
{
    if (modelnumber < 0 ||
        modelnumber >= static_cast<int>(models_.size()) ||
        !models_[modelnumber])
        return 1;

    models_[modelnumber]->SaveState();
    return 0;
}

int ExasimSolver::RestoreState(const int modelnumber)
{
    if (modelnumber < 0 ||
        modelnumber >= static_cast<int>(models_.size()) ||
        !models_[modelnumber])
        return 1;

    models_[modelnumber]->RestoreState();
    return 0;
}

int ExasimSolver::ClearSavedState(const int modelnumber)
{
    if (modelnumber < 0 ||
        modelnumber >= static_cast<int>(models_.size()) ||
        !models_[modelnumber])
        return 1;

    models_[modelnumber]->ClearSavedState();
    return 0;
}

std::vector<ExasimPoint> ExasimSolver::getInterfacePoints() const
{
    std::vector<ExasimPoint> interfacePoints;

    const int ngb = nfaces * ngf;
    if (ngb <= 0 || xdggint == nullptr)
        return interfacePoints;

    std::vector<double> xdggHost(ncx * ngb);
    TemplateCopytoHost(xdggHost.data(), xdggint, ncx * ngb, backend_);

    interfacePoints.reserve(ngb);
    for (int i = 0; i < ngb; ++i) {
        ExasimPoint p{};
        p.x = xdggHost[0 * ngb + i];
        p.y = (ncx > 1) ? xdggHost[1 * ngb + i] : 0.0;
        p.z = (ncx > 2) ? xdggHost[2 * ngb + i] : 0.0;
        interfacePoints.push_back(p);
    }

    return interfacePoints;
}

void ExasimSolver::getInterfaceFluxes(std::vector<double>& send_flux) const
{    
    send_flux.resize(ngf * nfaces * _ncuint);
    CSolution& model = *models_[interface_modelnumber_];
    model.disc.getInterfaceFluxesAtGaussPoints(send_flux.data(), xdggint, nlgint, faces, nfaces);
}

void ExasimSolver::setInterfaceFluxes(const std::vector<double>& recv_flux)
{
    CSolution& model = *models_[interface_modelnumber_];
    ArrayCopy(model.disc.sol.uext, recv_flux.data(), ngf * nfaces * model.disc.common.ncuext);
    model.disc.common.FextCall = _ibc + 1; // flag to indicate that external fields have been set
}

int ExasimSolver::RunAVDistanceFunction()
{
    vector<CSolution*> ptrs = ModelPointers(models_);
    avdistfunc(ptrs.data(), residual_outputs_.data(), nummodels_, backend_);
    return 0;
}

int ExasimSolver::InitializeSolution()
{
    for (int i = 0; i < nummodels_; i++) {
        if (restart_ > 0) {
            models_[i]->disc.common.currentstep = -1;
            models_[i]->ReadSolutions(backend_);
        }
        models_[i]->InitSolution(backend_);
    }

    return 0;
}

int ExasimSolver::InitializeSolution(const int i)
{
    if (restart_ > 0) {
        models_[i]->disc.common.currentstep = -1;
        models_[i]->ReadSolutions(backend_);
    }
    models_[i]->InitSolution(backend_);

    return 0;
}

int ExasimSolver::RunTimeDependent()
{
    int status = InitializeSolution();
    if (status) return status;

    return RunTimeDependent(0, models_[0]->disc.common.tsteps);
}

int ExasimSolver::RunTimeDependent(const int modelnumber)
{
    int status = InitializeSolution(modelnumber);
    if (status) return status;

    return RunTimeDependent(modelnumber, 0, models_[0]->disc.common.tsteps);
}

int ExasimSolver::RunTimeDependent(Int start, Int steps)
{
    dstype time = models_[0]->disc.common.time;
    for (Int istep = start; istep < start+steps; istep++) {
        for (int i = 0; i < nummodels_; i++) {
            models_[i]->disc.common.currentstep = istep;
            PreviousSolutions(models_[i]->disc.sol, models_[i]->solv.sys,
                              models_[i]->disc.common, backend_);
        }

        for (Int j = 0; j < models_[0]->disc.common.tstages; j++) {
            if (models_[0]->disc.common.mpiRank == 0)
                printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",
                       istep + 1, j + 1,
                       time + models_[0]->disc.common.dt[istep] *
                       models_[0]->disc.common.DIRKcoeff_t[j]);

            for (int i = 0; i < nummodels_; i++) {
                models_[i]->disc.common.currentstage = j;
                models_[i]->disc.common.time =
                    time + models_[i]->disc.common.dt[istep] *
                    models_[i]->disc.common.DIRKcoeff_t[j];

                UpdateSource(models_[i]->disc.sol, models_[i]->solv.sys,
                             models_[i]->disc.app, models_[i]->disc.driver_abi, models_[i]->disc.res,
                             models_[i]->disc.common, backend_);

                models_[i]->SteadyProblem(residual_outputs_[i], backend_);

                UpdateSolution(models_[i]->disc.sol, models_[i]->solv.sys,
                               models_[i]->disc.app, models_[i]->disc.driver_abi, models_[i]->disc.res,
                               models_[i]->disc.tmp, models_[i]->disc.common,
                               backend_);
            }
        }

        for (int i = 0; i < nummodels_; i++) {
            if (models_[i]->disc.common.compudgavg == 1) {
                ArrayAXPBY(models_[i]->disc.sol.udgavg,
                           models_[i]->disc.sol.udgavg,
                           models_[i]->disc.sol.udg, one, one,
                           models_[i]->disc.common.ndofudg1);
                ArrayAddScalar(&models_[i]->disc.sol.udgavg[
                    models_[i]->disc.common.ndofudg1], one, 1);
            }

            models_[i]->disc.computeAverageSolutionsOnBoundary();
            models_[i]->SaveSolutions(backend_);
            models_[i]->SaveQoI(backend_);
            if (models_[i]->vis.savemode > 0)
                models_[i]->SaveParaview(backend_);
            models_[i]->SaveSolutionsOnBoundary(backend_);
            if (models_[i]->disc.common.nce > 0)
                models_[i]->SaveOutputCG(backend_);
        }

        time = time + models_[0]->disc.common.dt[istep];
    }

    return 0;
}

int ExasimSolver::RunTimeDependent(const int i, Int start, Int steps)
{
    dstype time = models_[i]->disc.common.time;
    for (Int istep = start; istep < start+steps; istep++) {
        models_[i]->disc.common.currentstep = istep;
        PreviousSolutions(models_[i]->disc.sol, models_[i]->solv.sys,
                          models_[i]->disc.common, backend_);
        
        for (Int j = 0; j < models_[0]->disc.common.tstages; j++) {
            if (models_[0]->disc.common.mpiRank == 0)
                printf("\nTimestep :  %d,  Timestage :  %d,   Time : %g\n",
                       istep + 1, j + 1,
                       time + models_[0]->disc.common.dt[istep] *
                       models_[0]->disc.common.DIRKcoeff_t[j]);
            
            models_[i]->disc.common.currentstage = j;
            models_[i]->disc.common.time =
                time + models_[i]->disc.common.dt[istep] *
                models_[i]->disc.common.DIRKcoeff_t[j];

            UpdateSource(models_[i]->disc.sol, models_[i]->solv.sys,
                             models_[i]->disc.app, models_[i]->disc.driver_abi, models_[i]->disc.res,
                         models_[i]->disc.common, backend_);

            models_[i]->SteadyProblem(residual_outputs_[i], backend_);

            UpdateSolution(models_[i]->disc.sol, models_[i]->solv.sys,
                               models_[i]->disc.app, models_[i]->disc.driver_abi, models_[i]->disc.res,
                           models_[i]->disc.tmp, models_[i]->disc.common,
                           backend_);
            
        }
        
        if (models_[i]->disc.common.compudgavg == 1) {
            ArrayAXPBY(models_[i]->disc.sol.udgavg,
                       models_[i]->disc.sol.udgavg,
                       models_[i]->disc.sol.udg, one, one,
                       models_[i]->disc.common.ndofudg1);
            ArrayAddScalar(&models_[i]->disc.sol.udgavg[
                models_[i]->disc.common.ndofudg1], one, 1);
        }

        models_[i]->disc.computeAverageSolutionsOnBoundary();
        models_[i]->SaveSolutions(backend_);
        models_[i]->SaveQoI(backend_);
        if (models_[i]->vis.savemode > 0) models_[i]->SaveParaview(backend_);
        models_[i]->SaveSolutionsOnBoundary(backend_);
        if (models_[i]->disc.common.nce > 0) models_[i]->SaveOutputCG(backend_);
        
        time = time + models_[0]->disc.common.dt[istep];
    }

    return 0;
}

int ExasimSolver::RunPseudoTime()
{
    int status = InitializeSolution();
    if (status) return status;

    for (int i = 0; i < nummodels_; i++) {
        models_[i]->SteadyProblem_PTC(residual_outputs_[i], backend_);
    }
    return 0;
}

int ExasimSolver::RunPseudoTime(const int modelnumber)
{
    int status = InitializeSolution(modelnumber);
    if (status) return status;    
    models_[modelnumber]->SteadyProblem_PTC(residual_outputs_[modelnumber], backend_);

    return 0;
}

int ExasimSolver::RunSteady()
{
    int status = InitializeSolution();
    if (status) return status;

    for (int i = 0; i < nummodels_; i++) {
        models_[i]->SteadyProblem(residual_outputs_[i], backend_);

        models_[i]->SaveSolutions(backend_);    
        models_[i]->SaveQoI(backend_);
        if (models_[i]->vis.savemode > 0) models_[i]->SaveParaview(backend_); 
        models_[i]->SaveSolutionsOnBoundary(backend_);         
        if (models_[i]->disc.common.nce>0) models_[i]->SaveOutputCG(backend_);            
    }

    return 0;
}

int ExasimSolver::RunSteady(const int modelnumber)
{
    int status = InitializeSolution(modelnumber);
    if (status) return status;
    models_[modelnumber]->SteadyProblem(residual_outputs_[modelnumber], backend_);

    models_[modelnumber]->SaveSolutions(backend_);    
    models_[modelnumber]->SaveQoI(backend_);
    if (models_[modelnumber]->vis.savemode > 0) models_[modelnumber]->SaveParaview(backend_); 
    models_[modelnumber]->SaveSolutionsOnBoundary(backend_);         
    if (models_[modelnumber]->disc.common.nce>0) models_[modelnumber]->SaveOutputCG(backend_);            

    return 0;
}

int ExasimSolver::RunSolveProblemOrPostprocess()
{
    for (int i = 0; i < nummodels_; i++) {
        if (models_[i]->disc.common.runmode == 0) {
            if (restart_ > 0) {
                models_[i]->disc.common.currentstep = -1;
                models_[i]->ReadSolutions(backend_);
            }
            models_[i]->SolveProblem(residual_outputs_[i], backend_);
        }
        else if (models_[i]->disc.common.runmode == 1) {
            models_[i]->disc.common.currentstep = -1;
            models_[i]->ReadSolutions(backend_);
            if (models_[i]->disc.common.ncq > 0)
                models_[i]->disc.evalQ(backend_);
            models_[i]->disc.common.saveSolOpt = 1;
            models_[i]->SaveSolutions(backend_);
            models_[i]->SaveQoI(backend_);
            if (models_[i]->vis.savemode > 0)
                models_[i]->SaveParaview(backend_);
            models_[i]->SaveOutputCG(backend_);
        }
        else if (models_[i]->disc.common.runmode == 2) {
            for (Int istep = 0; istep < models_[i]->disc.common.tsteps; istep++) {
                if (((istep + 1) % models_[i]->disc.common.saveSolFreq) == 0) {
                    models_[i]->disc.common.currentstep = istep;
                    models_[i]->ReadSolutions(backend_);
                    if (models_[i]->disc.common.ncq > 0)
                        models_[i]->disc.evalQ(backend_);
                    models_[i]->SaveOutputCG(backend_);
                }
            }
        }
        else if (models_[i]->disc.common.runmode == 3) {
            for (Int istep = 0; istep < models_[i]->disc.common.tsteps; istep++) {
                if (((istep + 1) % models_[i]->disc.common.saveSolFreq) == 0) {
                    models_[i]->disc.common.currentstep = istep;
                    models_[i]->ReadSolutions(backend_);
                    models_[i]->SaveOutputCG(backend_);
                }
            }
        }
        else if (models_[i]->disc.common.runmode == 4) {
            models_[i]->disc.common.currentstep = -1;
            models_[i]->ReadSolutions(backend_);
            models_[i]->SolveProblem(residual_outputs_[i], backend_);
        }
        else if (models_[i]->disc.common.runmode == 5) {
            models_[i]->disc.evalQSer(backend_);
            models_[i]->disc.common.saveSolOpt = 1;
            string filename = models_[i]->disc.common.fileout + "_np" +
                              NumberToString(models_[i]->disc.common.mpiRank) + ".bin";
            writearray2file(filename, models_[i]->disc.sol.udg,
                            models_[i]->disc.common.ndofudg1, backend_);
        }
    }

    return 0;
}

void ExasimSolver::CloseOutputStreams()
{
    for (int i = 0; i < static_cast<int>(models_.size()) &&
                    i < static_cast<int>(residual_outputs_.size()); i++) {
        if (models_[i] &&
            models_[i]->disc.common.mpiRank == 0 &&
            models_[i]->disc.common.saveResNorm == 1 &&
            residual_outputs_[i].is_open()) {
            residual_outputs_[i].close();
        }
    }
}

void ExasimSolver::DestroyModels()
{
    if (!models_.empty() && models_[0]) {
        const int interfaceBackend = models_[0]->disc.common.backend;
        if (faces) {
            TemplateFree(faces, interfaceBackend);
            faces = nullptr;
        }
        if (xdgint) {
            TemplateFree(xdgint, interfaceBackend);
            xdgint = nullptr;
        }
        if (nlint) {
            TemplateFree(nlint, interfaceBackend);
            nlint = nullptr;
        }
        if (xdggint) {
            TemplateFree(xdggint, interfaceBackend);
            xdggint = nullptr;
        }
        if (nlgint) {
            TemplateFree(nlgint, interfaceBackend);
            nlgint = nullptr;
        }
    }

    for (auto& model : models_) {
        if (model) {
            delete[] model->disc.common.ncarray;
            model->disc.common.ncarray = nullptr;
            delete[] model->disc.sol.udgarray;
            model->disc.sol.udgarray = nullptr;
        }
    }

    models_.clear();
    model_abis_.clear();
    builtinmodelID_.clear();
    residual_outputs_.clear();
    filein_.clear();
    fileout_.clear();
    exasimpath_.clear();
    nummodels_ = 0;
    restart_ = 0;
    interface_modelnumber_ = -1;
}

void ExasimSolver::FinalizeEnvironment()
{
    if (owns_kokkos_ && Kokkos::is_initialized()) {
        Kokkos::finalize();
        owns_kokkos_ = false;
    }

#ifdef HAVE_MPI
    int mpiFinalized = 0;
    MPI_Finalized(&mpiFinalized);
    if (!mpiFinalized) {
        if (local_comm_ != MPI_COMM_NULL && local_comm_ != world_comm_) {
            MPI_Comm_free(&local_comm_);
            local_comm_ = world_comm_;
            EXASIM_COMM_LOCAL = local_comm_;
        }

        if (owns_mpi_) {
            MPI_Finalize();
            owns_mpi_ = false;
        }
    }

    world_comm_ = MPI_COMM_NULL;
    local_comm_ = MPI_COMM_NULL;
    EXASIM_COMM_WORLD = MPI_COMM_NULL;
    EXASIM_COMM_LOCAL = MPI_COMM_NULL;
#endif
}
