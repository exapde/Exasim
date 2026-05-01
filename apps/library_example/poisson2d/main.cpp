// apps/library_example/poisson2d/main.cpp
//
// Hand-written user-facing example. The math is in `my_model.hpp`
// (the `Poisson2D` struct). Everything else — MPI/Kokkos setup,
// argv parsing, DIRK time stepping, multi-domain coupling, file I/O —
// is in `<exasim/run.hpp>`'s `exasim::run<M>(argc, argv)`.

#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<Poisson2D>(argc, argv);
}
