# Exasim documentation

Exasim is a C++ library for high-order discontinuous Galerkin (DG)
PDE solvers, with first-class support for handwritten model code,
text2code symbolic-DSL code generation, and four backends
(CPU / GPU / MPI / MPI+GPU).

| Section | Audience | Contents |
|---|---|---|
| [Introduction](00-introduction.md) | Everyone | What Exasim is, capabilities, the three-axis decision (PDE × Solver × Mesh) |
| [Installation](01-installation.md) | First-time builders | Vendored + external deps, building text2code, building Exasim, per-platform install (macOS, Linux CPU, Linux + NVIDIA, Linux + AMD) |
| [Model contract](02-model-contract.md) | Active users | Reference for the `Model` C++ struct: every method, every signature, defaults, indexing layouts |
| [Internals](03-internals/) | Contributors | Test harness, baseline format, preprocessing internals, known-divergences |

For walkthroughs that exercise each combination of (PDE × Solver
× Mesh) end-to-end on the same Poisson problem, see
[`tutorial/`](../tutorial/README.md). The tutorial is the
recommended entry point for new users; this doc tree is
reference-shaped.

## The internals directory

| File | Topic |
|---|---|
| [`testing.md`](03-internals/testing.md) | ctest matrix, validate_codegen.sh harness, baseline gating |
| [`baselines.md`](03-internals/baselines.md) | Element-L2 baseline format, sidecar files, recording new baselines |
| [`architecture.md`](03-internals/architecture.md) | Preprocessing pipeline, DMD, ParMETIS, periodic node merging |
| [`known-divergences.md`](03-internals/known-divergences.md) | Open known-failures and their classification |

## Pre-existing references

The following older documents are still authoritative for legacy
content and have not been merged into the new structure:

- [`apps/library_example/README.md`](../apps/library_example/README.md) — the in-tree example apps and their port/regenerate harness.
- [`include/exasim/README.md`](../include/exasim/README.md) — public-header layout and contract reference.
- [`doc/Exasim.pdf`](../doc/Exasim.pdf) — original architecture paper.
