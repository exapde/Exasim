# 2. Methods overview

There are three independent choices when running Exasim:

| Dimension | Options |
|---|---|
| **How you author the math** | Hand-written `Model` C++ struct, or symbolic `pdemodel.txt` + text2code |
| **How you drive the solver** | Legacy CLI (`exasim::run<M>(argc, argv)`) or embedded library (`ExasimSolver<M>`) |
| **What backend it runs on** | CPU, GPU (CUDA / HIP), MPI, MPI+GPU |

These dimensions are independent — every authoring × driving × backend
combination works.

## Decision tree

```
                      ┌─ Standalone binary, runtime-configured?
                      │     → Legacy CLI driving
                      │       (read pdeapp.txt, mesh from disk, write dataout/)
Need a CLI tool? ─────┤
                      │
                      └─ Embed in a containing C++ program?
                            → ExasimSolver<M> driving
                              (mesh as flat C arrays, params as method calls,
                               solution from accessors, no disk round-trip)


                              ┌─ Hand-derive the math?
                              │     → Hand-written Model struct
                              │       (full control, learn the contract)
Where does the math come from?─┤
                              │
                              └─ Symbolic / quick-iteration?
                                    → Codegen path
                                      (write pdemodel.txt in a small DSL,
                                       text2code emits the Jacobians,
                                       my_model.hpp is generated)
```

## The four binaries each combination produces

For any concrete `<name>` and authoring choice, the build produces
four binaries (one per backend):

```
build_cpu/<name>          # CPU, single rank
build_gpu/<name>_gpu      # GPU (CUDA or HIP), single rank
build_mpi/<name>_mpi      # CPU, multi-rank MPI
build_mpi_gpu/<name>_mpi_gpu  # GPU + MPI, multi-rank
```

These come from the same source files; the difference is which
Kokkos build prefix is linked and which `EXASIM_*` defines are set
at compile time.

## What's the same across all combinations

- The model contract (`Flux`, `Source`, `Ubou`, `Initu`, `Tau`, …).
  The signatures are identical whether you wrote them by hand or
  text2code emitted them.
- The runtime kernels (`<exasim/solution.hpp>` instantiates the
  templated FEM internals on whatever `Model` you give it).
- The mesh, boundary, and parameter contract (the `Model` struct
  declares only the math; the surrounding `pdeapp.txt` or
  `ExasimSolver<M>` calls supply the rest).

## Concrete walkthroughs

For each combination there's a focused walkthrough:

| Combination | Doc |
|---|---|
| Hand-written + legacy | [`03-methods/hand-written-model.md`](03-methods/hand-written-model.md) |
| Hand-written + embedded | [`03-methods/embedded-facade.md`](03-methods/embedded-facade.md) (also covers codegen + embedded) |
| Codegen + legacy | [`03-methods/codegen-text2code.md`](03-methods/codegen-text2code.md) |
| Codegen + embedded | [`03-methods/embedded-facade.md`](03-methods/embedded-facade.md) §`load_pdeapp` |
| Backend-specific build & run | [`03-methods/cpu-gpu-mpi-mpigpu.md`](03-methods/cpu-gpu-mpi-mpigpu.md) |

## Embedding pros and cons

The embedded library API (`ExasimSolver<M>`) is the recommended
default for new code. It is:

- **In-memory throughout** — no `pdeapp.txt`, `datain/*.bin`, or
  `dataout/*.bin` round-trips. Your program builds the mesh as
  C arrays, passes them to the solver, and reads `udg()` /
  `uhat()` / `wdg()` accessors when the solve completes.
- **Type-safe** — boundary classifiers are C++ lambdas, not
  tinyexpr strings parsed at runtime. Misuse fails at compile time.
- **MPI-aware** — `set_mesh_distributed` lets each rank give its
  own slice of the global mesh. ParMETIS partitions for load
  balance internally.

The legacy CLI driving (`exasim::run<M>(argc, argv)`) is preferred
when:

- You want a runnable binary configured by a text file (HPC
  job scripts, parameter studies).
- You're calling Exasim from Julia / Python / MATLAB via subprocess.
- You're working with an existing `pdeapp.txt` workflow you don't
  want to port.

Both produce identical numerics. Same `Model` struct, same FEM
kernels — just different ways to feed them inputs and read outputs.

## Hand-written vs codegen pros and cons

Hand-written:
- Full control over the math, especially Jacobian derivations.
- No SymEngine learning curve.
- Custom physics that don't fit the DSL (e.g. table lookups,
  external solvers in your `Source` term).
- Easier to debug when residuals misbehave — you know what every
  line is doing.

Codegen:
- Concise — write the math once, derivatives are auto-generated.
- Less boilerplate than hand-writing dozens of pointwise
  functions for a complex PDE.
- Integrated with `text2code`'s SymEngine CSE so output is fast.

Both produce the same `my_model.hpp` shape; `<exasim/run.hpp>`
and `ExasimSolver<M>` don't know which path created it.
