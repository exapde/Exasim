// solve_square.cpp — single-file Exasim consumer.
//
// Builds a Cartesian quad mesh, populates the runtime config in C++,
// then drives CPreprocessing + CSolution<Poisson2D> directly. This is
// what <exasim/run.hpp> does internally; we replace it with code we
// can edit, instrument, or wrap in a tighter outer loop.
//
// No pdeapp.txt, no pdemodel.txt, no text2code. The only file the
// program writes on its way to the solve is the mesh (grid.txt).
//
// Tutorial: doc/tutorial_cpu_gpu_mpi.md §5.

// We pull the entire <exasim/run.hpp> preamble (preprocessor defines,
// `using namespace std;`, backend Common / Discretization / Solver /
// Preprocessing / Solution headers in the right order). We don't call
// `exasim::run<>()` — main() drives CPreprocessing + CSolution<M>
// directly. The backend headers expect this exact include order
// because they were written for the legacy `backend/Main/main.cpp`.
#include <exasim/run.hpp>
#include <exasim/model.hpp>              // exasim::ModelDefaults<>

#include <filesystem>
#include <fstream>
#include <string>

// ------------------------------------------------------------------
// 1. The PDE: Poisson 2D, hand-written. No text2code involved.
//    See doc/getting_started.md for a full annotation of every method.
// ------------------------------------------------------------------
struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;
    static constexpr int ncu    = 1;
    static constexpr int ncw    = 0;
    static constexpr int nco    = 0;
    static constexpr int nparam = 1;
    static constexpr auto disc  = exasim::Discretization::HDG;
    static constexpr int Nq     = ncu * (1 + nd);

    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double[], const double uq[],
              const double[], const double[], const double mu[],
              const double[], double) {
        f[0] = mu[0] * uq[1];
        f[1] = mu[0] * uq[2];
    }
    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double[],
                const double[], const double[], const double[],
                const double[], double) {
        constexpr double pi = 3.141592653589793;
        s[0] = 2.0 * pi * pi
             * Kokkos::sin(pi * x[0]) * Kokkos::sin(pi * x[1]);
    }
    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double[], const double[],
                     const double[], const double[], const double mu[],
                     const double[], double) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int, const double[], const double[],
                  const double[], const double[], const double uh[],
                  const double[], const double tau[],
                  const double[], const double[], double) {
        fb[0] = -tau[0] * uh[0];
    }
    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int, const double[],
                         const double[], const double[], const double[],
                         const double[], const double[], const double tau[],
                         const double[], const double[], double) {
        fb_uh[0] = -tau[0];
    }
    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double[], const double[], const double[]) {
        ui[0] = 0.0;
    }
    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int, const double x[], const double uq[],
              const double v[], const double w[], const double uh[],
              const double n[], const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0]*n[0] + f_local[1]*n[1] + tau[0]*(uq[0] - uh[0]);
    }
    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int, const double[], const double[],
              const double[], const double[], const double[],
              const double[], const double[], const double[],
              const double[], double) {
        ub[0] = 0.0;
    }
};

// ------------------------------------------------------------------
// 2. Generate an n×n quad mesh on the unit square as flat arrays.
//    Layout (matches what `meshFromArrays` expects):
//        p[d + 2*j]      = d-th coord of j-th vertex   (column-major)
//        t[v + 4*e]      = 0-based index of v-th corner of e-th element
//    Quad corners CCW; meshFromArrays infers elemtype=1.
// ------------------------------------------------------------------
struct Square {
    int nv, ne;
    std::vector<double> p;   // 2 × nv
    std::vector<int>    t;   // 4 × ne
};
static Square build_square(int n) {
    Square s;
    s.nv = (n + 1) * (n + 1);
    s.ne = n * n;
    s.p.resize(2 * s.nv);
    s.t.resize(4 * s.ne);
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            int idx = j * (n + 1) + i;
            s.p[2*idx + 0] = double(i) / n;
            s.p[2*idx + 1] = double(j) / n;
        }
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            int e = j * n + i;
            int v00 = j * (n + 1) + i;        // 0-based
            s.t[4*e + 0] = v00;
            s.t[4*e + 1] = v00 + 1;
            s.t[4*e + 2] = v00 + 1 + (n + 1);
            s.t[4*e + 3] = v00 + (n + 1);
        }
    return s;
}

// ------------------------------------------------------------------
// 3. Build the runtime config in C++. No pdeapp.txt — we populate the
//    PDE / InputParams / ParsedSpec structs directly and hand them to
//    the programmatic CPreprocessing constructor.
// ------------------------------------------------------------------
static void build_runtime_config(PDE& pde, InputParams& params,
                                 ParsedSpec& spec)
{
    pde.discretization = "hdg";
    pde.platform       = "cpu";
    pde.datapath       = ".";
    pde.datainpath     = "./datain";
    pde.dataoutpath    = "./dataout";
    // exasimpath must resolve `${exasimpath}/backend/Preprocessing/masternodes.bin`.
    // EXASIM_DIR env var lets callers point at an installed/source tree
    // without rebuilding; falls back to "." for in-tree runs.
    if (const char* d = std::getenv("EXASIM_DIR")) pde.exasimpath = d;
    else                                           pde.exasimpath = ".";
    pde.meshfile       = "";   // unused: mesh comes from arrays
    pde.modelfile      = "";
    pde.gendatain      = 1;
    pde.builtinmodelID = 1;
    pde.saveOutputs    = 0;       // HOT.7.4 — keep solution in disc.sol; no dataout/ writes
    pde.porder = 3;
    pde.pgauss = 6;
    pde.torder = 1;
    pde.nstage = 1;
    pde.tdep   = 0;
    pde.nd     = 2;       pde.ncu = 1;       pde.ncw = 0;
    pde.nc     = 3;
    pde.neb    = 4096;    pde.nfb = 8192;
    pde.ibs    = 1;
    pde.NewtonIter   = 20;     pde.NewtonTol   = 1e-6;
    pde.GMRESiter    = 200;    pde.GMRESrestart = 50;
    pde.GMREStol     = 1e-8;
    pde.tau          = {1.0};
    pde.dt           = {0.0};
    pde.physicsparam = {1.0};

    params.boundaryConditions = {1, 1, 1, 1};
    // Boundary predicates instead of tinyexpr strings (HOT.7.2). The
    // four predicates classify points on the four edges of [0,1]^2.
    // Order matches `boundaryConditions` above.
    params.boundaryPreds = {
        [](const double* x){ return std::abs(x[1])     < 1e-8; },   // y == 0
        [](const double* x){ return std::abs(x[0] - 1) < 1e-8; },   // x == 1
        [](const double* x){ return std::abs(x[1] - 1) < 1e-8; },   // y == 1
        [](const double* x){ return std::abs(x[0])     < 1e-8; },   // x == 0
    };
    params.curvedBoundaries     = {0, 0, 0, 0};   // none curved
    params.curvedBoundaryExprs  = {"", "", "", ""};
    params.tau                = {1.0};
    params.dt                 = {0.0};
    params.physicsParam       = {1.0};

    (void)spec;
}

int main(int argc, char** argv) {
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    int mpiprocs = 1, mpirank = 0;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
#else
    int mpiprocs = 1, mpirank = 0;
    (void)argc; (void)argv;
#endif

    Kokkos::initialize();
    {
        // No datain/ or dataout/ — preprocessing handoff goes through
        // `Preprocessed` (HOT.7.3); converged solution is read with
        // `sol.host_udg()` after SolveProblem (HOT.7.4).
        PDE         pde;
        InputParams params;
        ParsedSpec  spec;
        build_runtime_config(pde, params, spec);

        CPreprocessing preproc(pde, params, spec, mpirank, mpiprocs);

        // Mesh is built directly into the preprocessor — no grid.txt
        // hop. `meshFromArrays` does what `initializeMesh()` would do
        // off disk: copies p/t into the Mesh struct, infers elemtype
        // from (nd, nve), and applies the boundary / model fields
        // from `params` and `pde`.
        Square sq = build_square(/*n=*/16);
        preproc.mesh = meshFromArrays(sq.p.data(), sq.t.data(),
                                      sq.nv, sq.ne, /*nve=*/4, /*nd=*/2,
                                      preproc.params, preproc.pde);

        // Run preprocessing into an in-memory Preprocessed bundle —
        // no app.bin/master.bin/mesh.bin/sol.bin written. CSolution
        // consumes it directly via the matching constructor.
        // (HOT.7.3 — replaces SerialPreprocessing + CSolution(filein,…))
        auto pre = preproc.take();

        std::string fileout = pde.dataoutpath + "/out";
        int fileoffset = 0, gpuid = 0, backend = 0;
        CSolution<Poisson2D> sol(std::move(pre), fileout, pde.exasimpath,
                                 mpiprocs, mpirank, fileoffset,
                                 gpuid, backend, pde.builtinmodelID);

        sol.disc.common.nomodels = 1;
        sol.disc.common.ncarray  = new Int[1]{ sol.disc.common.nc };
        sol.disc.sol.udgarray    = new dstype*[1]{ &sol.disc.sol.udg[0] };

        std::ofstream resnorms;
        sol.SolveProblem(resnorms, backend);

        // HOT.7.4 — read the converged solution from memory rather
        // than from a file. host_udg() returns disc.sol.udg directly
        // (size = host_udg_size() = npe * nc * ne).
        const dstype* udg = sol.host_udg();
        Int udg_n = sol.host_udg_size();
        double maxabs = 0;
        for (Int i = 0; i < udg_n; ++i) {
            double v = std::abs((double)udg[i]);
            if (v > maxabs) maxabs = v;
        }
        if (mpirank == 0)
            std::cerr << "[solve_square] in-memory udg: " << udg_n
                      << " doubles, max|udg| = " << maxabs << "\n";

        delete[] sol.disc.common.ncarray;
        delete[] sol.disc.sol.udgarray;
    }
    Kokkos::finalize();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
