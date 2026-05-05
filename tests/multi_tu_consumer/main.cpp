// Multi-TU consumer test driver. Mirrors tutorial 04's shape (hand-
// written Model + exasim::run<Poisson2D>) but split across three
// translation units instead of one. The link step is the test: any
// non-inline free function defined in <exasim/run.hpp>'s transitive
// headers will be present in unit_a.o and unit_b.o (and main.o), and
// the linker will flag the duplicate.

#include <exasim/run.hpp>
#include "my_model.hpp"

int report_a();
int report_b();

int main(int argc, char** argv) {
    // Touch the helper TUs so they aren't dead-stripped before
    // linkage finishes; they are what surface the ODR violations.
    if (report_a() + report_b() != 55) return 2;
    return exasim::run<Poisson2D>(argc, argv);
}
