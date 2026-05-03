// Tutorial path 06 — hand-written Model + ExasimSolver<M> with a
// distributed mesh.
//
// Each rank builds its own slice of the global Cartesian quad mesh
// on the unit square. Vertex IDs in `t_local` are global node IDs;
// ParMETIS migrates the slice and repartitions for load balance
// inside `solver.solve(mpiprocs, mpirank)`.

#include <exasim/run.hpp>
#include <exasim/solver_facade.hpp>
#include <exasim/model.hpp>

#include "my_model.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

// Splits a global range [0, total) into `size` contiguous chunks.
// Rank `rank` gets `[offset, offset + count)`. Used to give each
// rank its initial slice of nodes and elements before ParMETIS
// repartitions.
static void compute_local_range(int total, int size, int rank,
                                int& count, int& offset) {
    int base = total / size;
    int rem  = total % size;
    if (rank < rem) {
        count  = base + 1;
        offset = rank * count;
    } else {
        count  = base;
        offset = rem * (base + 1) + (rank - rem) * base;
    }
}

int main(int argc, char** argv) {
    int mpiprocs = 1, mpirank = 0, shmrank = 0;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    EXASIM_COMM_WORLD = MPI_COMM_WORLD;
    EXASIM_COMM_LOCAL = MPI_COMM_WORLD;
    MPI_Comm_size(EXASIM_COMM_WORLD, &mpiprocs);
    MPI_Comm_rank(EXASIM_COMM_WORLD, &mpirank);
    MPI_Comm shmcomm;
    MPI_Comm_split_type(EXASIM_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &shmrank);
#else
    (void)argc; (void)argv;
#endif
#ifdef HAVE_CUDA
    int nd; cudaGetDeviceCount(&nd); cudaSetDevice(shmrank % nd);
#endif
#ifdef HAVE_HIP
    int nd; hipGetDeviceCount(&nd); hipSetDevice(shmrank % nd);
#endif
    (void)shmrank;

    Kokkos::initialize();
    {
        // Global mesh: 16x16 Cartesian quads on [0,1]^2.
        const int n  = 16;
        const int nv = (n + 1) * (n + 1);
        const int ne = n * n;

        // Each rank owns a contiguous slice of the global vertex and
        // element ranges. ParMETIS will repartition based on the
        // adjacency, so this initial split can be any reasonable one.
        int np_local, ne_local, node_off, elem_off;
        compute_local_range(nv, mpiprocs, mpirank, np_local, node_off);
        compute_local_range(ne, mpiprocs, mpirank, ne_local, elem_off);

        std::vector<double> p_local(2 * np_local);
        std::vector<int>    t_local(4 * ne_local);

        for (int li = 0; li < np_local; ++li) {
            int g = node_off + li;
            int j = g / (n + 1);
            int i = g % (n + 1);
            p_local[2 * li + 0] = double(i) / n;
            p_local[2 * li + 1] = double(j) / n;
        }
        // The t_local entries are GLOBAL vertex IDs. Each rank
        // references vertices that live in its own slice and on its
        // neighbors' slices; ParMETIS will fetch and redistribute
        // during partitioning.
        for (int le = 0; le < ne_local; ++le) {
            int ge  = elem_off + le;
            int j   = ge / n;
            int i   = ge % n;
            int v00 = j * (n + 1) + i;
            t_local[4 * le + 0] = v00;
            t_local[4 * le + 1] = v00 + 1;
            t_local[4 * le + 2] = v00 + 1 + (n + 1);
            t_local[4 * le + 3] = v00 + (n + 1);
        }

        exasim::ExasimSolver<Poisson2D> solver;
        solver.set_mesh_distributed(p_local.data(), t_local.data(),
                                    np_local, ne_local, /*nve=*/4,
                                    /*np_global=*/nv, /*ne_global=*/ne);
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[1])     < 1e-8; });
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[0] - 1) < 1e-8; });
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[1] - 1) < 1e-8; });
        solver.add_boundary(/*tag=*/1,
            [](const double* x){ return std::abs(x[0])     < 1e-8; });
        solver.set_polynomial_order(3);
        solver.set_quadrature_order(6);
        solver.set_physics_params({1.0});
        solver.solve(mpiprocs, mpirank);

        const double* udg   = solver.udg();
        const Int     udg_n = solver.udg_size();

        double maxabs_local = 0.0;
        for (Int i = 0; i < udg_n; ++i) {
            double v = std::abs(udg[i]);
            if (v > maxabs_local) maxabs_local = v;
        }
        double maxabs = maxabs_local;
#ifdef HAVE_MPI
        MPI_Allreduce(&maxabs_local, &maxabs, 1, MPI_DOUBLE, MPI_MAX,
                      EXASIM_COMM_WORLD);
#endif
        if (mpirank == 0) {
            std::printf("[tutorial_06] mpiprocs=%d, max|udg| = %.5f\n",
                        mpiprocs, maxabs);
        }
    }
    Kokkos::finalize();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
