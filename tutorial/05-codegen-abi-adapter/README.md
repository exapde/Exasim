# 05 — Codegen + AbiAdapter (legacy `cput2cEXASIM`)

Authoring: SymEngine DSL in `pdemodel.txt`.
Driving: `cput2cEXASIM` (or `cpumpit2cEXASIM`, `gput2cEXASIM`,
`gpumpit2cEXASIM`). These binaries are built from
`backend/Main/main.cpp`, which calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`. The
`AbiAdapter` marker tells the FEM templates to dispatch every kernel
through `dlsym` against `libpdemodel{serial,cuda,hip}.{so,dylib}`.

## Status

This path is not currently exercised by the test matrix. A spot-check
with `cput2cEXASIM` on `poisson2d` produces NaN residual at Newton
iteration 0 — the dispatch through libpdemodel disagrees with the
templated runtime in some way. Treat as deprecated until a
regression test brings it back to green.

For new code use path 03 ([`03-codegen-cli/`](../03-codegen-cli/README.md))
or path 04 ([`04-codegen-facade/`](../04-codegen-facade/README.md)).

## Build + run (for reference)

```bash
cmake -S install -B build_cpu -DWITH_TEXT2CODE=ON ...
cmake --build build_cpu --target cput2cEXASIM
cd apps/library_example/poisson2d_codegen
bash ../regenerate.sh poisson2d
$EXASIM/build_cpu/cput2cEXASIM ./pdeapp.txt
```

See [`../../docs/03-methods/abi-adapter.md`](../../docs/03-methods/abi-adapter.md)
for the full target list.
