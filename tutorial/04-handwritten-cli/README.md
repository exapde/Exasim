# 04 — Hand-written Model with the legacy CLI

This section solves Poisson 2D on the unit square using a hand-
written `Poisson2D` C++ struct and Exasim's legacy command-line
interface.

## Files

- `my_model.hpp` defines the `Poisson2D` struct that satisfies the
  contract in `<exasim/model.hpp>`. The struct provides the flux,
  source, boundary, and initial-condition methods that the templated
  FEM internals call.
- `main.cpp` instantiates the runtime by calling
  `exasim::run<Poisson2D>(argc, argv)`. The function reads
  `pdeapp.txt`, runs preprocessing, solves, and writes the
  converged state to `dataout/`.
- `pdeapp.txt` configures the runtime: discretization order,
  solver tolerances, boundary conditions, mesh path, physics
  parameters, and so on.
- `grid.bin` is the binary mesh file consumed by the runtime.
- `CMakeLists.txt` builds the binary as `tutorial_04_handwritten_cli`.

## Build and run

```bash
cd $EXASIM
cmake --build build --target tutorial_04_handwritten_cli
cd $EXASIM/tutorial/04-handwritten-cli
mkdir -p datain dataout
$EXASIM/build/tutorial_04_handwritten_cli ./pdeapp.txt
```

The runtime writes `dataout/outudg_np0.bin` and friends. The
test driver `tutorial/run-all.sh` calls these steps and verifies
that the output files exist.
