// Out-of-tree install consumer test.
//
// Includes the public Model contract header and instantiates one of
// its concept checks against a downstream-defined `Poisson2D`. The
// test is the *build* — find_package(Exasim) must resolve, the
// public header chain (B1) must resolve from <prefix>/include/...,
// and the result must compile.
//
// We deliberately do NOT pull in <exasim/run.hpp>. That entry point
// drags in `backend/Model/KokkosDrivers.cpp`, which in turn pulls
// per-model codegen kernels (`KokkosFlux.cpp`, ...) that are not
// installable as-is — they're text2code output, generated per
// model by the consumer. A real downstream consumer of run.hpp
// runs text2code first; that's a runtime workflow concern, not an
// install-layout concern. B-Test1 only validates the install
// layout, so it stays narrow.

#include <exasim/model.hpp>
#include "my_model.hpp"

static_assert(exasim::is_model_v<Poisson2D>,
              "Poisson2D should satisfy the Exasim Model contract");

int main() { return 0; }
