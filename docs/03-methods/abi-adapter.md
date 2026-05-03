# AbiAdapter — `cput2cEXASIM` and family

The original Exasim driving path. `backend/Main/main.cpp` calls
`exasim::run<exasim::detail::AbiAdapter>(argc, argv)`. The
`AbiAdapter` marker type instructs the FEM templates to dispatch
every kernel call through the legacy `libpdemodel.hpp` ABI — i.e.
`dlsym` lookup against `libpdemodel{serial,cuda,hip}.{so,dylib}`.

text2code generates the library from `pdemodel.txt`; the binary
loads it at runtime.

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
bash ../regenerate.sh poisson2d            # generate libpdemodel*
$EXASIM/build/cput2cEXASIM ./pdeapp.txt
```

The binary's `RPATH` points at `apps/library_example/<name>_codegen/`
(per-example output dir), so it picks up the freshly-generated
library.

## Status

The AbiAdapter path is **not exercised by the test matrix**. The 12
codegen examples in `apps/library_example/<name>_codegen/` use the
templated `GeneratedModel` path (their `main.cpp` calls
`exasim::run<GeneratedModel>(argc, argv)`), not the AbiAdapter.

A spot-check on `cput2cEXASIM` against `poisson2d` produces NaN
residual at Newton iteration 0 — the dispatch through libpdemodel
disagrees with the runtime in some way. Diagnosis pending; treat
this path as deprecated until a regression test brings it back to
green. The templated path
([`hand-written-model.md`](hand-written-model.md),
[`embedded-facade.md`](embedded-facade.md)) is the supported entry
point for new code.

## When you'd want this path

- Reproducing pre-templated-path Exasim runs.
- Migrating an old project that depended on `cput2cEXASIM`'s exact
  CLI surface.
- Bisecting whether a regression came from the runtime kernels or
  the templated dispatch — the AbiAdapter path bypasses the
  templated kernel<M> wrappers entirely.
