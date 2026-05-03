# 01 — Hand-written model + legacy CLI

Authoring: hand-written C++ struct `Poisson2D` in `my_model.hpp`.
Driving: `exasim::run<Poisson2D>(argc, argv)` reading `pdeapp.txt`.

Source: `apps/library_example/poisson2d/`. The struct is in
`my_model.hpp`; the three-line `main.cpp` is:

```cpp
#include <exasim/run.hpp>
#include "my_model.hpp"
int main(int argc, char** argv) {
    return exasim::run<Poisson2D>(argc, argv);
}
```

## Run

```bash
cmake --build build_cpu --target poisson2d_template
cd apps/library_example/poisson2d
mkdir -p datain dataout
$EXASIM/build_cpu/poisson2d_template ./pdeapp.txt
ls dataout/        # outudg_np0.bin, outuhat_np0.bin, outqoi.txt, ...
```

For MPI:

```bash
cd apps/library_example/poisson2d
mpirun -np 2 $EXASIM/build_mpi/poisson2d_template_mpi ./pdeapp.txt
```

## Read the solution

The runtime writes binaries that the validation harness reads.
For programmatic access from the same binary, switch to the facade
([`../02-handwritten-facade/`](../02-handwritten-facade/README.md)).
