# Library examples

Each subdirectory here is a small consumer of Exasim as a header-only
template library. The runtime entry point is
`exasim::CSolution<MyModel>` (wrapped by `exasim::run<MyModel>(argc, argv)`
in `<exasim/run.hpp>`); `MyModel` is either hand-written or emitted
by `text2code` from a `pdemodel.txt`.

For a guided tutorial on writing a new PDE from scratch (the
hand-written path), see `doc/getting_started.md`. For the API
reference (every method, every signature), see `doc/header_only_api.md`.

## Two flavors

**Hand-written**: `poisson2d/` — the user wrote `my_model.hpp` directly.
This is the pattern you'd follow for a new PDE you're inventing or
porting.

**Auto-generated**: `<name>_codegen/` — the user wrote a `pdemodel.txt`
in the SymEngine DSL, ran the retargeted `text2code`, and got
`my_model.hpp` for free. This is the path most existing apps in
`apps/poisson/` and `apps/navierstokes/` take. As of this writing, 12
of 14 source apps work this way (see `install/CMakeLists.txt` for the
two unsupported cases and the reasons).

## Generating a codegen example from a source app

```bash
bash apps/library_example/port_codegen.sh <name> --build --validate
```

This:
1. Locates the source app in `apps/poisson/<name>/` or `apps/navierstokes/<name>/`.
2. Copies `main.cpp` (generic — same for every codegen example) and
   `pdeapp.txt` (with mpiprocs=1, gendatain=1, and any `repeat(<dt>, n)`
   in the dt schedule collapsed to a single `<dt>` so the validation
   gate stays under a few seconds per app).
3. Calls `regenerate.sh <name>` to copy `grid.bin` + `pdemodel.txt` from
   the source app, run `text2code`, and drop the resulting
   `my_model.hpp` into the codegen dir.
4. Builds the `<name>_codegen` cmake target.
5. Records a baseline by running legacy `cput2cEXASIM` on the same
   inputs, saving `outqoi.txt` and `md5 + bins` to
   `baseline/<name>_serial/`.
6. Validates the codegen output against the recorded baseline:
   - QoI byte-diff if the model defines QoIvolume; else
   - Bin md5 match; else
   - Numerical-close fallback (rel diff < 1e-6, the typical NewtonTol).

Once a port is committed (`main.cpp`, `pdeapp.txt`, `my_model.hpp`,
plus `add_codegen_example(<name>)` in `install/CMakeLists.txt`), the
gate runs as `bash apps/library_example/validate_codegen.sh <name>`.

`grid.bin`, `pdemodel.txt`, `datain/`, `dataout/`, and the runtime-only
companion bins (`xdg.bin`, `udg.bin`, `vdg.bin`, `wdg.bin`,
`partition.bin`) are gitignored under `*_codegen/` — they're
regenerable.

## What the generic main.cpp does

It's the same MPI/Kokkos init + multi-domain coupling boilerplate as
`backend/Main/main.cpp`, with the legacy `CSolution<>` (defaulting to
`AbiAdapter`) replaced by `CSolution<GeneratedModel>` (or
`CSolution<Poisson2D>` in the hand-written example). The user's `main`
references the templated FEM internals through `<exasim/solution.hpp>`
and friends; the templates inline the user's pointwise math from
`my_model.hpp` at compile time.

For external (out-of-tree) consumers, `find_package(Exasim REQUIRED)`
exposes the `Exasim::headers` INTERFACE target — see `doc/header_only_api.md`.
