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
        // Build an n×n quad mesh on [0,1]² as flat arrays.
        // ------------------------------------------------------------
        const int n = 16;
        const int nv = (n + 1) * (n + 1);
        const int ne = n * n;
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

        // ------------------------------------------------------------
        // The full façade flow: mesh + boundaries + params + solve.
        // ------------------------------------------------------------
        exasim::ExasimSolver<Poisson2D> solver;
        solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4);
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
