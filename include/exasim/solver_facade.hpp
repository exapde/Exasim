// SPDX-License-Identifier: see LICENSE
//
// <exasim/solver_facade.hpp> — high-level fluent wrapper.
//
// Bundles HOT.7.0–7.4 into a single class:
//
//   ExasimSolver<MyModel> solver;
//   solver.set_mesh(p, t, np, ne, nve);
//   solver.add_boundary(/*tag=*/1, [](const double* x){ return std::abs(x[1])<1e-8; });
//   solver.set_physics_params({1.0});
//   solver.set_polynomial_order(3);
//   solver.solve();
//   const dstype* udg = solver.udg();
//
// All the underlying machinery — `CPreprocessing(PDE, InputParams,
// ParsedSpec, ...)`, `meshFromArrays`, `take()`, `CSolution<M>(
// Preprocessed&&, ...)`, `SolveProblem`, `host_udg()` — is reachable
// directly for users who want finer control. This façade just sets
// HDG-friendly defaults and saves typing.
//
// The Model template parameter `M` is the same `exasim::ModelDefaults<>`-
// derived struct the templated kernels see. Compile-time constants
// (`M::nd`, `M::ncu`, `M::ncw`, …) are propagated into the PDE
// struct automatically.
//
// Lifetime: ExasimSolver owns the underlying `CSolution<M>` (held by
// `unique_ptr`). The accessors (`udg()`, `uhat()`) return non-owning
// pointers into `disc.sol.*` and are valid until the ExasimSolver
// destructs.

#pragma once

#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <exasim/run.hpp>                         // pulls common.h, all backend headers, `using namespace std;`
#include <exasim/model.hpp>                       // exasim::ModelDefaults, is_model_v
#include "../../backend/Preprocessing/structs.hpp"
#include "../../backend/Preprocessing/buildstructs.hpp"

namespace exasim {

template <class M>
class ExasimSolver {
    static_assert(is_model_v<M>,
                  "ExasimSolver<M>: M must satisfy the Model contract.");

public:
    ExasimSolver()
    {
        // HDG-friendly defaults. The user's PDE may override any of
        // these via pde() / params() before solve().
        pde_.discretization = "hdg";
        pde_.platform       = "cpu";
        pde_.gendatain      = 1;
        pde_.builtinmodelID = 1;
        pde_.saveOutputs    = 0;        // façade users default to in-memory
        pde_.porder         = 1;
        pde_.pgauss         = 2;
        pde_.torder         = 1;
        pde_.nstage         = 1;
        pde_.tdep           = 0;
        pde_.NewtonIter     = 20;
        pde_.NewtonTol      = 1e-6;
        pde_.GMRESiter      = 200;
        pde_.GMRESrestart   = 50;
        pde_.GMREStol       = 1e-8;
        pde_.tau            = {1.0};
        pde_.dt             = {0.0};
        pde_.neb            = 4096;
        pde_.nfb            = 8192;
        pde_.ibs            = 1;

        // Compile-time dimensions from the Model.
        pde_.nd  = M::nd;
        pde_.ncu = M::ncu;
        pde_.ncw = M::ncw;
        pde_.nc  = M::ncu * (1 + M::nd);   // ModelD layout
        pde_.physicsparam.assign(/*size=*/std::max(1, M::nparam), 0.0);

        // exasimpath default: env var `EXASIM_DIR` if set, else ".".
        // setcommonstruct expects this to point at the source/install
        // tree containing `backend/Preprocessing/{master,gauss}nodes.bin`.
        if (const char* d = std::getenv("EXASIM_DIR")) pde_.exasimpath = d;
        else                                            pde_.exasimpath = ".";
        pde_.datapath    = ".";
        pde_.datainpath  = "./datain";
        pde_.dataoutpath = "./dataout";
        pde_.modelfile   = "";              // unused; Model is hand-supplied
        pde_.meshfile    = "";              // unused; mesh comes from set_mesh()
    }

    // --- Required setters --------------------------------------------------

    // Vertex coords + element connectivity. `p` is `nd × np` doubles
    // (column-major); `t` is `nve × ne` 0-based ints (column-major).
    // Currently the façade infers `nd` from `M::nd` and asserts on
    // mismatch.
    void set_mesh(const double* p, const int* t, int np, int ne, int nve)
    {
        mesh_p_.assign(p, p + (size_t)M::nd * np);
        mesh_t_.assign(t, t + (size_t)nve * ne);
        mesh_np_   = np;
        mesh_ne_   = ne;
        mesh_nve_  = nve;
        distributed_ = false;
    }

    // HOT.7.8 — distributed mesh: each rank supplies its slice. The
    // connectivity uses **global** node indices in [0, np_global);
    // global element IDs are implicit (rank r owns
    // [elmdist[r], elmdist[r+1])). ParMETIS will repartition for
    // load balance, so any reasonable initial split works (rows,
    // contiguous blocks, etc.).
    void set_mesh_distributed(const double* p_local, const int* t_local,
                              int np_local, int ne_local, int nve,
                              int np_global, int ne_global)
    {
        mesh_p_.assign(p_local, p_local + (size_t)M::nd * np_local);
        mesh_t_.assign(t_local, t_local + (size_t)nve   * ne_local);
        mesh_np_   = np_local;
        mesh_ne_   = ne_local;
        mesh_nve_  = nve;
        mesh_np_global_ = np_global;
        mesh_ne_global_ = ne_global;
        distributed_    = true;
    }

    // Adds a boundary tag with a typed predicate. Tags are 1-based by
    // legacy convention — the order they're added becomes the
    // `int ib` discriminator the Model's fbou_hdg / fbou methods see.
    void add_boundary(int tag, BoundaryPred pred)
    {
        params_.boundaryConditions.push_back(tag);
        params_.boundaryPreds.push_back(std::move(pred));
        params_.curvedBoundaries.push_back(0);
        params_.curvedBoundaryExprs.push_back("");
    }

    // Curved boundary: same as add_boundary but supplies a level-set
    // function `f(x)` whose zero set is the curved surface (see
    // doc/tutorial_cpu_gpu_mpi.md §5 for the projector formula).
    void add_curved_boundary(int tag, BoundaryPred classify, BoundaryLevelSet level_set)
    {
        params_.boundaryConditions.push_back(tag);
        params_.boundaryPreds.push_back(std::move(classify));
        params_.curvedBoundaries.push_back(1);
        params_.curvedBoundaryExprs.push_back("");
        params_.curvedBoundaryLevelSets.push_back(std::move(level_set));
    }

    // --- Common knobs ------------------------------------------------------

    void set_polynomial_order(int porder)        { pde_.porder = porder; }
    void set_quadrature_order(int pgauss)        { pde_.pgauss = pgauss; }
    void set_physics_params(std::vector<double> p) { pde_.physicsparam = std::move(p); }
    void set_newton_tol(double tol)              { pde_.NewtonTol = tol; }
    void set_gmres_tol(double tol)               { pde_.GMREStol  = tol; }
    void set_tau(std::vector<double> t)          { pde_.tau = std::move(t); }
    void set_exasim_dir(std::string dir)         { pde_.exasimpath = std::move(dir); }

    // HOT.7.9 — load runtime configuration from a legacy pdeapp.txt
    // file. Same parser the file-driven CSolution<M>(filein,...)
    // path uses; populates pde_/params_/spec_ exactly as the legacy
    // path would. After this call you can still override individual
    // fields via pde() / params() / spec() escape hatches before
    // calling solve().
    //
    // Mesh: if neither set_mesh() nor set_mesh_distributed() has
    // been called by the time solve() runs, the facade falls through
    // to initializeMesh / initializeParMesh on `pde_.meshfile`. So
    // a complete legacy app can be driven by:
    //
    //     ExasimSolver<MyModel> solver;
    //     solver.load_pdeapp("pdeapp.txt");
    //     solver.solve(mpiprocs, mpirank);
    //
    // (no set_mesh* needed — grid.bin is read off disk).
    void load_pdeapp(const std::string& filename, int mpirank = 0)
    {
        params_ = parseInputFile(filename, mpirank);
        pde_    = initializePDE(params_, mpirank);
        if (!pde_.modelfile.empty())
            spec_ = TextParser::parseFile(make_path(pde_.datapath, pde_.modelfile));
        if (spec_.exasimpath.empty()) spec_.exasimpath = pde_.exasimpath;
        pdeapp_loaded_ = true;
    }

    // Path layout. `dir == ""` keeps in-memory mode (saveOutputs = 0,
    // no files written). Setting a non-empty dir flips saveOutputs
    // back on so the runtime writes outudg_np<r>.bin etc. there.
    void set_output_dir(std::string dir)
    {
        if (dir.empty()) {
            pde_.saveOutputs = 0;
        } else {
            pde_.dataoutpath = std::move(dir);
            pde_.saveOutputs = 1;
        }
    }

    // --- Escape hatches ----------------------------------------------------
    PDE&         pde()    noexcept { return pde_; }
    InputParams& params() noexcept { return params_; }
    ParsedSpec&  spec()   noexcept { return spec_; }
    const PDE&         pde()    const noexcept { return pde_; }
    const InputParams& params() const noexcept { return params_; }
    const ParsedSpec&  spec()   const noexcept { return spec_; }

    // --- Run ---------------------------------------------------------------

    // Construct the in-memory pipeline and solve. Returns 0 on
    // success. After this call, udg() / uhat() / wdg() return host
    // pointers into the converged state (CPU only for now; GPU/HIP
    // wiring lands in HOT.7.6).
    int solve(int mpiprocs = 1, int mpirank = 0)
    {
        const bool have_inmem_mesh = (mesh_np_ > 0);
        const bool have_file_mesh  = !pde_.meshfile.empty();
        if (!have_inmem_mesh && !have_file_mesh)
            error("ExasimSolver::solve(): no mesh — call set_mesh / set_mesh_distributed "
                  "or load_pdeapp(...) with a meshfile.");
        if (params_.boundaryConditions.empty())
            error("ExasimSolver::solve(): no boundary tags — call add_boundary "
                  "or load_pdeapp(...) with boundaryconditions set.");

        pde_.mpiprocs = mpiprocs;

        // Backend matches the binary's compile-time defines. Same
        // convention <exasim/run.hpp> uses for the legacy main path.
        int backend = 0;
#ifdef HAVE_OPENMP
        backend = 1;
#endif
#ifdef HAVE_CUDA
        backend = 2;
#endif
#ifdef HAVE_HIP
        backend = 3;
#endif
        // When saveOutputs == 0 the binary won't open any output
        // files, so an empty `fileout` is fine and we skip creating
        // the dataout directory.
        const std::string fileout = (pde_.saveOutputs != 0)
                                        ? (pde_.dataoutpath + "/out")
                                        : std::string{};
        const int gpuid      = 0;
        const int fileoffset = 0;

        CPreprocessing preproc(pde_, params_, spec_, mpirank, mpiprocs);

        if (mpiprocs == 1) {
            if (have_inmem_mesh) {
                preproc.mesh = meshFromArrays(mesh_p_.data(), mesh_t_.data(),
                                              mesh_np_, mesh_ne_, mesh_nve_, M::nd,
                                              preproc.params, preproc.pde);
                auto pre = preproc.take();
                pre.save_outputs = (pde_.saveOutputs != 0);
                solver_ = std::make_unique<CSolution<M>>(
                    std::move(pre), fileout, pde_.exasimpath,
                    mpiprocs, mpirank, fileoffset, gpuid, backend, pde_.builtinmodelID);
            } else {
                // Legacy file path: SerialPreprocessing reads the
                // mesh from pde.meshfile, writes datain bins,
                // CSolution<M>(filein,...) reads them back.
                std::filesystem::create_directories(pde_.datainpath);
                if (pde_.saveOutputs != 0) std::filesystem::create_directories(pde_.dataoutpath);
                preproc.SerialPreprocessing();
                const std::string filein = pde_.datainpath + "/";
                solver_ = std::make_unique<CSolution<M>>(
                    filein, fileout, pde_.exasimpath,
                    mpiprocs, mpirank, fileoffset, gpuid, backend, pde_.builtinmodelID);
            }
        } else {
#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
            // Populate preproc.mesh: distributed-from-arrays if the
            // caller supplied per-rank slices, else read the user-
            // provided global mesh from pde.meshfile and let ParMETIS
            // partition. Either way takeParallel() does the rest in
            // memory (no per-rank datain bins on disk).
            if (have_inmem_mesh) {
                if (!distributed_)
                    error("ExasimSolver::solve(): mpiprocs > 1 requires set_mesh_distributed(...). "
                          "set_mesh(...) is the single-rank API.");
                preproc.mesh = meshFromArraysDistributed(
                    mesh_p_.data(), mesh_t_.data(),
                    mesh_np_, mesh_ne_, mesh_nve_, M::nd,
                    mesh_np_global_, mesh_ne_global_,
                    preproc.params, preproc.pde);
            } else {
                preproc.mesh = initializeParMesh(preproc.params, preproc.spec,
                                                 preproc.pde, EXASIM_COMM_LOCAL);
            }

            // EXASIM_FACADE_INMEMORY_MPI=0 falls back to the legacy
            // file ABI (writes per-rank datain bins, reads them back).
            // Kept as an escape hatch while the in-memory takeParallel
            // path is the new default. HOT.7.13.
            const char* inmem_env = std::getenv("EXASIM_FACADE_INMEMORY_MPI");
            const bool use_files = (inmem_env && inmem_env[0] == '0');
            if (!use_files) {
                if (pde_.saveOutputs != 0) std::filesystem::create_directories(pde_.dataoutpath);
                auto pre = preproc.takeParallel(EXASIM_COMM_LOCAL);
                pre.save_outputs = (pde_.saveOutputs != 0);
                solver_ = std::make_unique<CSolution<M>>(
                    std::move(pre), fileout, pde_.exasimpath,
                    mpiprocs, mpirank, fileoffset, gpuid, backend, pde_.builtinmodelID);
            } else {
                std::filesystem::create_directories(pde_.datainpath);
                if (pde_.saveOutputs != 0) std::filesystem::create_directories(pde_.dataoutpath);
                preproc.ParallelPreprocessing(EXASIM_COMM_LOCAL);
                const std::string filein = pde_.datainpath + "/";
                solver_ = std::make_unique<CSolution<M>>(
                    filein, fileout, pde_.exasimpath,
                    mpiprocs, mpirank, fileoffset, gpuid, backend, pde_.builtinmodelID);
            }
#else
            error("ExasimSolver::solve(): mpiprocs > 1 requires HAVE_PARMETIS + HAVE_MPI.");
#endif
        }

        // CSolution expects the multi-model bookkeeping arrays even
        // for a single-model run (see <exasim/run.hpp> ~line 358).
        solver_->disc.common.nomodels = 1;
        solver_->disc.common.ncarray  = new Int[1]{ solver_->disc.common.nc };
        solver_->disc.sol.udgarray    = new dstype*[1]{ &solver_->disc.sol.udg[0] };

        std::ofstream resnorms;
        solver_->SolveProblem(resnorms, backend);
        solved_ = true;
        return 0;
    }

    // --- Access ------------------------------------------------------------

    // Returns a host-side pointer to the converged buffer. On GPU
    // backends the underlying disc.sol.* live on device; the first
    // call after solve() copies the bytes into a member-cached
    // host vector and returns its data().
    const dstype* udg()  { return host_buffer_(udg_buf_,  solver_ ? solver_->host_udg()  : nullptr,
                                                solver_ ? solver_->host_udg_size()  : 0); }
    Int           udg_size()  const { return solver_ ? solver_->host_udg_size()  : 0; }
    const dstype* uhat() { return host_buffer_(uhat_buf_, solver_ ? solver_->host_uhat() : nullptr,
                                                solver_ ? solver_->host_uhat_size() : 0); }
    Int           uhat_size() const { return solver_ ? solver_->host_uhat_size() : 0; }
    const dstype* wdg()  { return host_buffer_(wdg_buf_,  solver_ ? solver_->host_wdg()  : nullptr,
                                                solver_ ? solver_->host_wdg_size()  : 0); }
    Int           wdg_size()  const { return solver_ ? solver_->host_wdg_size()  : 0; }

    // Direct access to the underlying CSolution<M> for advanced
    // workflows (e.g. inspecting disc.common, calling disc.evalQ,
    // re-solving with updated physics_params).
    CSolution<M>*       solver()       noexcept { return solver_.get(); }
    const CSolution<M>* solver() const noexcept { return solver_.get(); }

    ~ExasimSolver()
    {
        if (solver_) {
            delete[] solver_->disc.common.ncarray;
            delete[] solver_->disc.sol.udgarray;
        }
    }

    // Non-copyable (CSolution<M> owns malloc'd C arrays through raw
    // pointers — copying would double-free).
    ExasimSolver(const ExasimSolver&) = delete;
    ExasimSolver& operator=(const ExasimSolver&) = delete;
    ExasimSolver(ExasimSolver&&) = default;
    ExasimSolver& operator=(ExasimSolver&&) = default;

private:
    PDE         pde_;
    InputParams params_;
    ParsedSpec  spec_;

    std::vector<double> mesh_p_;
    std::vector<int>    mesh_t_;
    int  mesh_np_         = 0;
    int  mesh_ne_         = 0;
    int  mesh_nve_        = 0;
    int  mesh_np_global_  = 0;
    int  mesh_ne_global_  = 0;
    bool distributed_     = false;
    bool pdeapp_loaded_   = false;

    std::unique_ptr<CSolution<M>> solver_;
    bool solved_ = false;

    // Host-side caches for GPU runs: disc.sol.{udg,uh,wdg} live on
    // device, so udg()/uhat()/wdg() copy them to host on first call.
    // For CPU backends we return the pointer directly (no copy).
    std::vector<dstype> udg_buf_, uhat_buf_, wdg_buf_;

    // Lazy device → host copy on GPU; pass-through on CPU. Returns
    // a host-resident pointer of `n` elements (or nullptr if n==0).
    const dstype* host_buffer_(std::vector<dstype>& cache,
                               const dstype* device_or_host_ptr, Int n)
    {
        if (!solver_ || n == 0 || device_or_host_ptr == nullptr) return nullptr;
        const int backend = solver_->disc.common.backend;
        if (backend < 2) return device_or_host_ptr;        // CPU / OpenMP
#ifdef HAVE_GPU
        if ((Int)cache.size() != n) cache.resize(n);
        // TemplateCopytoHost expects non-const pointers; the device
        // memory it sources from isn't being mutated, so the cast
        // is sound.
        TemplateCopytoHost(cache.data(),
                           const_cast<dstype*>(device_or_host_ptr),
                           n, backend);
        return cache.data();
#else
        return device_or_host_ptr;
#endif
    }
};

} // namespace exasim
