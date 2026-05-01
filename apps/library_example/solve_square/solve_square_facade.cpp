// solve_square_facade.cpp — same problem as solve_square.cpp but
// driven through the high-level `ExasimSolver<M>` façade (HOT.7.5).
//
// Compare the two side by side:
//
//   solve_square.cpp           — manual `CPreprocessing(pde, params,
//                                spec, ...) → take() → CSolution<M>`
//                                control flow, ~270 lines.
//   solve_square_facade.cpp    — `ExasimSolver<M>` wraps all of the
//                                above, ~80 lines.
//
// The façade reads compile-time `M::nd / ncu / ncw / nparam` and
// fills in HDG-friendly defaults (porder=1, pgauss=2, NewtonTol=1e-6,
// saveOutputs=0, etc.). The user supplies just the mesh, boundary
// predicates, and physics params.

#include <exasim/run.hpp>
#include <exasim/model.hpp>
#include <exasim/solver_facade.hpp>

#include <fstream>

// ------------------------------------------------------------------
// Same Poisson2D model as solve_square.cpp.
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

int main(int argc, char** argv) {
    int mpiprocs = 1, mpirank = 0, shmrank = 0;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);

    // Per-node rank — used to bind one GPU per rank.
    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
#else
    (void)argc; (void)argv;
#endif

    // Bind this rank to a GPU before Kokkos::initialize so Kokkos
    // picks the right device. Same `shmrank % deviceCount` rule
    // run.hpp uses.
#ifdef HAVE_CUDA
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) { std::cerr << "no CUDA devices\n"; return 1; }
    cudaSetDevice(shmrank % deviceCount);
#endif
#ifdef HAVE_HIP
    int deviceCount = 0;
    hipGetDeviceCount(&deviceCount);
    if (deviceCount <= 0) { std::cerr << "no HIP devices\n"; return 1; }
    hipSetDevice(shmrank % deviceCount);
#endif

    Kokkos::initialize();
    {
        // ------------------------------------------------------------
        // Build an n×n quad mesh on [0,1]². Single-rank gets the full
        // mesh; multi-rank does HOT.7.8's distributed mesh — each rank
        // owns a horizontal stripe of `n_per_rank` rows.
        //
        // Global numbering: nodes laid out left-to-right, bottom-up;
        // elements likewise. Rank r owns:
        //     rows  [r*n_per_rank, (r+1)*n_per_rank)
        // So local node IDs in the slice include rows
        //     [r*n_per_rank, (r+1)*n_per_rank]    (closed at top — shared edge with rank r+1)
        // and the slice's `t_local` references nodes by **global** ID.
        // ------------------------------------------------------------
        const int n  = 16;
        const int nv = (n + 1) * (n + 1);
        const int ne = n * n;

        exasim::ExasimSolver<Poisson2D> solver;
        if (mpiprocs == 1) {
            std::vector<double> p(2 * nv);
            std::vector<int>    t(4 * ne);
            for (int j = 0; j <= n; ++j)
                for (int i = 0; i <= n; ++i) {
                    int idx = j * (n + 1) + i;
                    p[2*idx + 0] = double(i) / n;
                    p[2*idx + 1] = double(j) / n;
                }
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i) {
                    int e = j * n + i;
                    int v00 = j * (n + 1) + i;
                    t[4*e + 0] = v00;
                    t[4*e + 1] = v00 + 1;
                    t[4*e + 2] = v00 + 1 + (n + 1);
                    t[4*e + 3] = v00 + (n + 1);
                }
            solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4);
        } else {
            // Distributed: each rank builds its row slice. We assume
            // `n` is divisible by mpiprocs; for the example mesh
            // (n=16) this works for 1, 2, 4, 8, 16 ranks.
            const int n_per_rank = n / mpiprocs;
            const int row_lo     =  mpirank      * n_per_rank;        // first owned row
            const int row_hi     = (mpirank + 1) * n_per_rank;        // one past last
            const int ne_local   = n_per_rank * n;
            // Local nodes: rows [row_lo, row_hi]   (inclusive at top → shared with rank+1)
            const int np_local   = (n + 1) * (row_hi - row_lo + 1);

            std::vector<double> p_loc(2 * np_local);
            std::vector<int>    t_loc(4 * ne_local);

            // Local node → global node ID:
            //     gid(i, j_global)  =  j_global * (n+1) + i
            // Local index  loc(i, j_local) = j_local * (n+1) + i,
            // with j_global = j_local + row_lo.
            for (int j_local = 0; j_local <= row_hi - row_lo; ++j_local) {
                int j_global = j_local + row_lo;
                for (int i = 0; i <= n; ++i) {
                    int loc = j_local * (n + 1) + i;
                    p_loc[2*loc + 0] = double(i) / n;
                    p_loc[2*loc + 1] = double(j_global) / n;
                }
            }
            for (int j_local = 0; j_local < n_per_rank; ++j_local) {
                int j_global = j_local + row_lo;
                for (int i = 0; i < n; ++i) {
                    int e = j_local * n + i;
                    // Connectivity uses GLOBAL node IDs (so neighboring
                    // ranks' connectivity matches up at the shared row).
                    int v00 = j_global * (n + 1) + i;
                    t_loc[4*e + 0] = v00;
                    t_loc[4*e + 1] = v00 + 1;
                    t_loc[4*e + 2] = v00 + 1 + (n + 1);
                    t_loc[4*e + 3] = v00 + (n + 1);
                }
            }
            solver.set_mesh_distributed(p_loc.data(), t_loc.data(),
                                        np_local, ne_local, /*nve=*/4,
                                        nv, ne);
        }
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[1])     < 1e-8; });   // y=0
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[0] - 1) < 1e-8; });   // x=1
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[1] - 1) < 1e-8; });   // y=1
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[0])     < 1e-8; });   // x=0
        solver.set_polynomial_order(3);
        solver.set_quadrature_order(6);
        solver.set_physics_params({1.0});
        solver.solve(mpiprocs, mpirank);

        // ------------------------------------------------------------
        // Read the converged state directly from memory.
        // ------------------------------------------------------------
        const dstype* udg = solver.udg();
        Int           udg_n = solver.udg_size();
        double maxabs = 0;
        for (Int i = 0; i < udg_n; ++i) {
            double v = std::abs((double)udg[i]);
            if (v > maxabs) maxabs = v;
        }
        if (mpirank == 0)
            std::cerr << "[solve_square_facade] udg: " << udg_n
                      << " doubles, max|udg| = " << maxabs << "\n";
    }
    Kokkos::finalize();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
