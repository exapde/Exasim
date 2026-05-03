# 05 — Codegen + AbiAdapter (`cput2cEXASIM` and family)

The original Exasim driving path. `backend/Main/main.cpp` calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`. The
`AbiAdapter` marker tells the FEM templates to dispatch every
kernel through the legacy `libpdemodel.hpp` ABI — i.e. linked
symbols in `libpdemodel{serial,cuda,hip}.{so,dylib}` resolved at
load time.

The binary's RPATH points at `backend/Model/`, so text2code must
write its output there (not the per-target dir).

## Generate + run

```bash
cd $EXASIM/apps/library_example/poisson2d_codegen
$EXASIM/build/text2code ./pdeapp.txt   # writes backend/Model/libpdemodel*.dylib
$EXASIM/build/cput2cEXASIM ./pdeapp.txt
```

For MPI / GPU / MPI+GPU: `cpumpit2cEXASIM`, `gput2cEXASIM`,
`gpumpit2cEXASIM`. These are gated on `WITH_TEXT2CODE=ON` at
configure time.

## Note on text2code output

The other tutorial paths use `regenerate.sh`, which calls
`text2code --out-dir <example_dir>`. The output goes to the
example's own directory and the example's binary RPATH-loads
from there. The AbiAdapter binaries
(`cput2cEXASIM`, etc.) do NOT honor that — their RPATH is
hardcoded to `backend/Model/`. Always run `text2code` without
`--out-dir` before running the AbiAdapter binaries.

Reference: [`../../docs/03-methods/abi-adapter.md`](../../docs/03-methods/abi-adapter.md).
