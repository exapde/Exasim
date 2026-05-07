// Multi-TU consumer test, translation unit A.
//
// Includes the public Exasim header and defines a callable. The test
// is the link step of unit_a.o + unit_b.o + main.o: any non-inline
// free function defined in <exasim/run.hpp>'s transitive header set
// will appear in both unit_a.o and unit_b.o and the linker will
// report a duplicate-symbol error.
//
// The compute path is never executed; we just need both TUs to pull
// the public headers' definitions through the same way a real
// downstream multi-file consumer would.

#include <exasim/run.hpp>

int report_a() { return 42; }
