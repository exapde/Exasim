# 06 — Handwritten model with the distributed-mesh embedded solver

The user writes the model header (`my_model.hpp`) by hand, then
each MPI rank constructs its own slice of the global mesh and
hands the slice to `solver.set_mesh_distributed(...)`. ParMETIS
migrates the slice and repartitions for load balance inside
`solver.solve(mpiprocs, mpirank)`. There is no `pdeapp.txt`, no
`grid.bin`, and no `text2code` step. The converged solution is
read back through `solver.udg()` per rank, and an `MPI_Allreduce`
collects the global maximum.

This is the embedded path for users who want a true
distributed-memory solve where each rank only ever materializes
its share of the mesh — required when the global mesh is too
large to fit on a single rank.

## Files

- `my_model.hpp` — C++ struct `Poisson2D` implementing the
  `<exasim/model.hpp>` contract.
- `main.cpp` — distributed embedded driver. Each rank builds a
  contiguous slice of a global 16×16 Cartesian quad mesh on
  `[0,1]²`, configures the solver, calls `solve(mpiprocs, mpirank)`,
  reduces the per-rank `max|udg|` across `MPI_COMM_WORLD`, and
  prints the global maximum on rank 0.
- `CMakeLists.txt` — out-of-tree build for
  `tutorial_06_handwritten_distributed`.

The runtime expects to find the master-element node tables under
`$EXASIM/backend/Preprocessing/`, so the binary needs the
environment variable `EXASIM_DIR` to be set when it runs. The
section is meaningful primarily under `mpirun -np N`; on serial
build dirs the binary still compiles and runs as a single rank
with `mpiprocs = 1`.

## Build and run

### CPU (single rank)

```bash
cd $EXASIM
cmake --build build --target tutorial_06_handwritten_distributed
EXASIM_DIR=$EXASIM $EXASIM/build/tutorial_06_handwritten_distributed
```

### GPU (single rank)

```bash
cd $EXASIM
cmake --build build_gpu --target tutorial_06_handwritten_distributed
EXASIM_DIR=$EXASIM $EXASIM/build_gpu/tutorial_06_handwritten_distributed
```

### MPI

```bash
cd $EXASIM
cmake --build build_mpi --target tutorial_06_handwritten_distributed
EXASIM_DIR=$EXASIM mpirun -np 2 \
    $EXASIM/build_mpi/tutorial_06_handwritten_distributed
```

### MPI+GPU

```bash
cd $EXASIM
cmake --build build_mpi_gpu --target tutorial_06_handwritten_distributed
EXASIM_DIR=$EXASIM mpirun -np 2 \
    $EXASIM/build_mpi_gpu/tutorial_06_handwritten_distributed
```

The expected output is `max|udg| = 3.14158` (the global maximum
absolute value of the converged mixed state for the manufactured
solution `u = sin(πx) sin(πy)`).

## File contents

### `my_model.hpp`

```cpp
#pragma once

#include <exasim/model.hpp>

struct Poisson2D : exasim::ModelDefaults<Poisson2D> {
    static constexpr int nd     = 2;
    static constexpr int ncu    = 1;
    static constexpr int ncw    = 0;
    static constexpr int nco    = 0;
    static constexpr int nparam = 1;

    static constexpr auto disc = exasim::Discretization::HDG;
    static constexpr int Nq = ncu * (1 + nd);

    KOKKOS_INLINE_FUNCTION static
    void flux(double f[], const double /*x*/[], const double uq[],
              const double /*v*/[], const double /*w*/[], const double mu[],
              const double /*uinf*/[], double /*t*/) {
        const double mu0  = mu[0];
        const double udg2 = uq[1];
        const double udg3 = uq[2];
        f[0] = mu0 * udg2;
        f[1] = mu0 * udg3;
    }

    KOKKOS_INLINE_FUNCTION static
    void source(double s[], const double x[], const double /*uq*/[],
                const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                const double /*uinf*/[], double /*t*/) {
        const double xdg1 = x[0];
        const double xdg2 = x[1];
        s[0] = Kokkos::sin(xdg1 * 3.141592653589793)
             * Kokkos::sin(xdg2 * 3.141592653589793)
             * 1.973920880217872E+1;
    }

    KOKKOS_INLINE_FUNCTION static
    void initu(double ui[], const double /*x*/[],
               const double /*uinf*/[], const double /*mu*/[]) {
        ui[0] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void flux_jac_uq(double f_uq[], const double /*x*/[],
                     const double /*uq*/[],
                     const double /*v*/[], const double /*w*/[],
                     const double mu[], const double /*uinf*/[],
                     double /*t*/) {
        for (int k = 0; k < ncu * nd * Nq; ++k) f_uq[k] = 0.0;
        f_uq[1 * (ncu * nd) + 0] = mu[0];
        f_uq[2 * (ncu * nd) + 1] = mu[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg(double fb[], int /*ib*/,
                  const double /*x*/[],  const double /*uq*/[],
                  const double /*v*/[],  const double /*w*/[],  const double uh[],
                  const double /*n*/[],  const double tau[],
                  const double /*mu*/[], const double /*uinf*/[],
                  double /*t*/) {
        fb[0] = -tau[0] * uh[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou_hdg_jac_uh(double fb_uh[], int /*ib*/,
                         const double /*x*/[],  const double /*uq*/[],
                         const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
                         const double /*n*/[],  const double tau[],
                         const double /*mu*/[], const double /*uinf*/[],
                         double /*t*/) {
        fb_uh[0] = -tau[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void fbou(double fb[], int /*ib*/,
              const double x[],  const double uq[],
              const double v[],  const double w[],  const double uh[],
              const double n[],  const double tau[],
              const double mu[], const double uinf[], double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }

    KOKKOS_INLINE_FUNCTION static
    void ubou(double ub[], int /*ib*/,
              const double /*x*/[],  const double /*uq*/[],
              const double /*v*/[],  const double /*w*/[],  const double /*uh*/[],
              const double /*n*/[],  const double /*tau*/[],
              const double /*mu*/[], const double /*uinf*/[],
              double /*t*/) {
        ub[0] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION static
    void vis_scalars(double s[], const double /*x*/[], const double uq[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[0];
        s[1] = uq[1] + uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void vis_vectors(double s[], const double /*x*/[], const double uq[],
                     const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                     const double /*uinf*/[], double /*t*/) {
        s[0] = uq[1];
        s[1] = uq[2];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_volume(double s[], const double x[], const double uq[],
                    const double /*v*/[], const double /*w*/[], const double /*mu*/[],
                    const double /*uinf*/[], double /*t*/) {
        const double t1 = 3.141592653589793;
        const double t2 = Kokkos::sin(t1 * x[0]);
        const double t3 = Kokkos::sin(t1 * x[1]);
        const double uexact = t2 * t3;
        s[0] = (uq[0] - uexact) * (uq[0] - uexact);
        s[1] = uq[0];
    }

    KOKKOS_INLINE_FUNCTION static
    void qoi_boundary(double fb[], int /*ib*/,
                      const double x[],  const double uq[],
                      const double v[],  const double w[],  const double uh[],
                      const double n[],  const double tau[],
                      const double mu[], const double uinf[],
                      double t) {
        double f_local[ncu * nd];
        flux(f_local, x, uq, v, w, mu, uinf, t);
        fb[0] = f_local[0] * n[0] + f_local[1] * n[1]
              + tau[0] * (uq[0] - uh[0]);
    }
};
```

`Poisson2D` inherits `exasim::ModelDefaults<Poisson2D>` (CRTP),
which supplies zero-fill defaults for every method on the model
contract. The HDG branch consumes `flux`, `source`,
`flux_jac_uq`, `fbou_hdg`, and `fbou_hdg_jac_uh`. The numeric
constant `1.973920880217872E+1` in `source` is `2π²` to enough
digits to match `text2code`'s SymEngine output exactly.

### `main.cpp`

```cpp
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
```

`compute_local_range(total, size, rank, count, offset)` is a tiny
helper that splits `[0, total)` into `size` contiguous chunks; the
first `total % size` ranks get one extra element so the slices
exactly cover the global range.

The `HAVE_MPI` block initializes MPI, stores the world
communicator on `EXASIM_COMM_WORLD` and `EXASIM_COMM_LOCAL`, and
splits a node-local communicator with `MPI_COMM_TYPE_SHARED` so
each rank knows its position inside its node. The `HAVE_CUDA` and
`HAVE_HIP` blocks bind one device per rank using that local rank
index.

The mesh-build loop uses `compute_local_range` twice: once for the
node range (each rank stores coordinates for `np_local` global
nodes starting at `node_off`) and once for the element range (each
rank stores incidence for `ne_local` global elements starting at
`elem_off`). The vertex IDs in `t_local` are global IDs, not local
indices into `p_local`. ParMETIS will fetch and redistribute
vertices during partitioning.

`solver.set_mesh_distributed(p_local, t_local, np_local, ne_local,
nve, np_global, ne_global)` hands each rank's slice to the solver.
The trailing `np_global` and `ne_global` parameters tell the
solver the size of the global mesh so it can size its global ID
maps. Boundary tagging, polynomial order, quadrature order, and
physics parameters are configured the same way as the in-memory
case in section 05.

`solver.solve(mpiprocs, mpirank)` runs the parallel Newton/GMRES
loop. Internally ParMETIS repartitions the mesh based on the
adjacency, exchanges halo elements, and the assembly and solver
loops run on each rank's locally owned elements.

After `solve()`, each rank's `solver.udg()` returns its locally
owned slice of the converged mixed state. The driver scans the
slice for the local maximum absolute value and uses
`MPI_Allreduce(MPI_MAX)` to collect the global maximum.

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.16)

set(_target tutorial_06_handwritten_distributed)

add_executable(${_target} main.cpp)
target_compile_features(${_target} PRIVATE cxx_std_17)
target_compile_definitions(${_target} PRIVATE _TEXT2CODE)
target_include_directories(${_target} PRIVATE
    ${EXASIM_DIR}/include
    "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
target_link_libraries(${_target} PRIVATE
    Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
link_metis(${_target})
set_target_properties(${_target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${Model_LIB_DIR}")
```

Identical to section 05's `CMakeLists.txt` modulo target name. The
target compiles `main.cpp` against the Exasim headers and links it
against Kokkos, BLAS/LAPACK, ParMETIS (via `link_metis`), and the
`libpdemodelserial.{so,dylib}` placeholder under `backend/Model/`.
