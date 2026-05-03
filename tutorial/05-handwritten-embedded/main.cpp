// Tutorial section 05 — handwritten model + ExasimSolver<M> with
// an in-memory mesh, single rank only.
//
// The mesh, boundary tagging, polynomial order, and physics
// parameters are set on the solver object. After `solver.solve()`
// the converged solution is read back via `solver.udg()`. For a
// distributed-memory variant that uses `set_mesh_distributed`, see
// section 06.

#include <exasim/run.hpp>            // pulls common preamble + kokkos + namespace std
#include <exasim/solver_facade.hpp>
#include <exasim/model.hpp>

#include "my_model.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

int main(int argc, char** argv) {
    (void)argc; (void)argv;

    Kokkos::initialize();
    {
        // Build a Cartesian quad mesh on [0,1]^2 with 16x16 elements.
        const int n  = 16;
        const int nv = (n + 1) * (n + 1);
        const int ne = n * n;

        std::vector<double> p(2 * nv);
        std::vector<int>    t(4 * ne);
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                int idx = j * (n + 1) + i;
                p[2 * idx + 0] = double(i) / n;
                p[2 * idx + 1] = double(j) / n;
            }
        }
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                int e   = j * n + i;
                int v00 = j * (n + 1) + i;
                t[4 * e + 0] = v00;
                t[4 * e + 1] = v00 + 1;
                t[4 * e + 2] = v00 + 1 + (n + 1);
                t[4 * e + 3] = v00 + (n + 1);
            }
        }

        exasim::ExasimSolver<Poisson2D> solver;
        solver.set_mesh(p.data(), t.data(), nv, ne, /*nve=*/4);
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
        solver.solve();

        const double* udg   = solver.udg();
        const Int     udg_n = solver.udg_size();

        double maxabs = 0.0;
        for (Int i = 0; i < udg_n; ++i) {
            double v = std::abs(udg[i]);
            if (v > maxabs) maxabs = v;
        }
        std::printf("[tutorial_05] udg: %lld doubles, max|udg| = %.5f\n",
                    static_cast<long long>(udg_n), maxabs);
    }
    Kokkos::finalize();
    return 0;
}
