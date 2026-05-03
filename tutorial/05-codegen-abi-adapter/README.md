# 05 — Codegen Model with the AbiAdapter driver

This section solves Poisson 2D using the codegen authoring path
together with the legacy `cput2cEXASIM` family of binaries. The
runtime is `backend/Main/main.cpp`, which calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`. The
`AbiAdapter` marker type instructs the FEM templates to dispatch
every kernel through the `libpdemodel.hpp` ABI: linked symbols
inside `libpdemodel{serial,cuda,hip}.{so,dylib}` are resolved at
load time.

The AbiAdapter binaries are built by Exasim itself when the CMake
option `WITH_TEXT2CODE=ON` is set; this section reuses those
binaries rather than building one of its own.

## Files

- `pdemodel.txt` describes the PDE in the SymEngine DSL.
- `pdeapp.txt` is the runtime configuration consumed by the
  AbiAdapter binary.
- `grid.bin` is the binary mesh.

## Build and run

Configure Exasim with `-DWITH_TEXT2CODE=ON` so that
`cput2cEXASIM` is built. The AbiAdapter binary's `RPATH` is
hard-coded to `backend/Model/`, so this path runs `text2code`
without `--out-dir` to populate that location.

```bash
cd $EXASIM
cmake -S install -B build_cpu -DWITH_TEXT2CODE=ON ...
cmake --build build_cpu --target cput2cEXASIM
cd $EXASIM/tutorial/05-codegen-abi-adapter
$EXASIM/build/text2code ./pdeapp.txt          # writes backend/Model/libpdemodel*
mkdir -p datain dataout
$EXASIM/build_cpu/cput2cEXASIM ./pdeapp.txt
```

The AbiAdapter MPI and GPU variants are named `cpumpit2cEXASIM`,
`gput2cEXASIM`, and `gpumpit2cEXASIM`. Each is built when
`WITH_TEXT2CODE=ON` and the matching backend is enabled.
