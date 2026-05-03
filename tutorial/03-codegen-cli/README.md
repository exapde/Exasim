# 03 — Codegen Model with the legacy CLI

This section solves Poisson 2D on the unit square using the
codegen authoring path: the math lives in `pdemodel.txt` (a
small SymEngine DSL), and `text2code` produces `my_model.hpp`
plus the runtime kernel libraries.

## Files

- `pdemodel.txt` describes the PDE in the SymEngine DSL.
  `text2code` reads it and emits `my_model.hpp` with the
  `GeneratedModel` struct, plus `libpdemodelserial.{so,dylib}`,
  the run-time kernel library that the legacy ABI links against.
- `main.cpp` is the same three-line driver as path 01, with the
  templated model swapped from `Poisson2D` to `GeneratedModel`.
- `pdeapp.txt` configures the runtime (identical structure to
  path 01).
- `grid.bin` is the binary mesh consumed by the runtime.
- `CMakeLists.txt` builds the binary as `tutorial_03_codegen_cli`.
  The build requires `my_model.hpp` and `libpdemodelserial.*` to
  exist in this directory; the run script generates them with
  `text2code` before invoking `cmake --build`.

## Build and run

```bash
cd $EXASIM
$EXASIM/build/text2code --out-dir $EXASIM/tutorial/03-codegen-cli \
    $EXASIM/tutorial/03-codegen-cli/pdeapp.txt
cmake --build build_cpu --target tutorial_03_codegen_cli
cd $EXASIM/tutorial/03-codegen-cli
mkdir -p datain dataout
$EXASIM/build_cpu/tutorial_03_codegen_cli ./pdeapp.txt
```
