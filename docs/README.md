# Exasim documentation

Exasim is a C++ library for high-order discontinuous Galerkin (DG)
PDE solvers, with first-class support for hand-written model code,
text2code symbolic-DSL code generation, and four backends
(CPU / GPU / MPI / MPI+GPU).

This documentation is organized into five sections:

| Section | Audience | Contents |
|---|---|---|
| [Introduction](00-introduction.md) | Everyone | What Exasim is, capabilities, the authoring × driving × backend matrix |
| [Installation](01-installation.md) | First-time builders | Vendored + external deps, building text2code, building Exasim |
| [Methods overview](02-methods-overview.md) | New users | Decision tree for which authoring + driving combination to pick |
| [03. Methods](03-methods/) | Active users | One walkthrough per concrete combination |
| [04. Internals](04-internals/) | Contributors | Test harness, baseline format, preprocessing internals, known issues |

## The methods directory

| File | Combination |
|---|---|
| [`hand-written-model.md`](03-methods/hand-written-model.md) | Write a `Model` C++ struct directly + `<exasim/run.hpp>` CLI |
| [`embedded-facade.md`](03-methods/embedded-facade.md) | `ExasimSolver<M>` library API — embed Exasim in a containing program, no `pdeapp.txt` round-trip |
| [`codegen-text2code.md`](03-methods/codegen-text2code.md) | Symbolic `pdemodel.txt` → text2code → `my_model.hpp` |
| [`cpu-gpu-mpi-mpigpu.md`](03-methods/cpu-gpu-mpi-mpigpu.md) | The four backend variants and how to build / run each |

## The internals directory

| File | Topic |
|---|---|
| [`testing.md`](04-internals/testing.md) | ctest matrix, validate_codegen.sh harness, baseline gating |
| [`baselines.md`](04-internals/baselines.md) | Element-L2 baseline format, sidecar files, recording new baselines |
| [`architecture.md`](04-internals/architecture.md) | Preprocessing pipeline, DMD, ParMETIS, periodic node merging |
| [`known-divergences.md`](04-internals/known-divergences.md) | Open known-failures and their classification |

## Pre-existing references

The following older documents are still authoritative for legacy
content and have not been merged into the new structure:

- [`apps/library_example/README.md`](../apps/library_example/README.md) — the in-tree example apps and their port/regenerate harness.
- [`include/exasim/README.md`](../include/exasim/README.md) — public-header layout and contract reference.
- [`doc/Exasim.pdf`](../doc/Exasim.pdf) — original architecture paper.

The migration is incremental — older `doc/*.md` files redirect to
`docs/` equivalents.
