// Tutorial path 04 — hand-written Model + legacy CLI driver.
//
// `exasim::run<Poisson2D>(argc, argv)` reads `pdeapp.txt` from
// `argv[1]`, runs preprocessing (mesh, master element, DMD), solves
// the discretized PDE, and writes the converged solution to
// `dataout/`. The same code path is what `cput2cEXASIM` and the
// other in-tree examples use.

#include <exasim/run.hpp>
#include "my_model.hpp"

int main(int argc, char** argv) {
    return exasim::run<Poisson2D>(argc, argv);
}
