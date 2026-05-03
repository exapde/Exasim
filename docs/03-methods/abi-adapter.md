# AbiAdapter — `cput2cEXASIM` and family

The original Exasim driving path. `backend/Main/main.cpp` calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`. The
`AbiAdapter` marker type instructs the FEM templates to dispatch
every kernel call through the legacy `libpdemodel.hpp` ABI: the
binary RPATH-links `libpdemodel{serial,cuda,hip}.{so,dylib}` from
`backend/Model/`.

text2code generates the library from `pdemodel.txt`; the binary
loads it at startup.

## Targets

The legacy targets are defined in `install/CMakeLists.txt` and gated
on `WITH_TEXT2CODE=ON`:

| Target | Backend | MPI |
|---|---|---|
| `cpuEXASIM` | CPU built-in models | no |
| `cput2cEXASIM` | CPU + text2code lib | no |
| `cpumpit2cEXASIM` | CPU + text2code lib | yes |
| `gpuEXASIM` | GPU built-in models | no |
| `gput2cEXASIM` | GPU + text2code lib | no |
| `gpumpit2cEXASIM` | GPU + text2code lib | yes |

Build:

```bash
cmake -S install -B build -DWITH_TEXT2CODE=ON ...
cmake --build build --target cput2cEXASIM
```

## Invocation

```bash
cd apps/library_example/poisson2d_codegen
$EXASIM/build/text2code ./pdeapp.txt       # writes backend/Model/libpdemodel*
$EXASIM/build/cput2cEXASIM ./pdeapp.txt
```

The binary's RPATH is `backend/Model/`. text2code without
`--out-dir` writes there. The codegen + facade tutorials use
`regenerate.sh` which passes `--out-dir <example_dir>` for parallel-
safe per-example libs; cput2cEXASIM does not honor that, so use the
plain text2code invocation above.

## Test coverage

`tutorial:run-all:cpu` runs `cput2cEXASIM` against `poisson2d`
([`../../tutorial/05-codegen-abi-adapter/`](../../tutorial/05-codegen-abi-adapter/README.md)).
The AbiAdapter MPI / GPU variants
(`cpumpit2cEXASIM`, `gput2cEXASIM`, `gpumpit2cEXASIM`) compile when
`WITH_TEXT2CODE=ON` but are not currently part of the test matrix.
