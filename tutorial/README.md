# Tutorial

Solve the same problem (Poisson 2D on the unit square, manufactured
solution `u = sin(πx)·sin(πy)`) through every supported authoring
× driving combination, ordered from **most automated** (drop two
text files at a prebuilt binary) to **least automated** (manage
mesh, boundaries, and partitioning by hand from C++ across MPI
ranks).

Start with [`00-overview/README.md`](00-overview/README.md) for
the naming convention and a guide to picking a section.

| Section | Title | Walkthrough |
|---|---|---|
| 00 | Overview and naming convention | [`00-overview/`](00-overview/README.md) |
| 01 | Generated model with the prebuilt solver | [`01-generated-prebuilt/`](01-generated-prebuilt/README.md) |
| 02 | Generated model with the CLI wrapper | [`02-generated-cli/`](02-generated-cli/README.md) |
| 03 | Generated model with the embedded solver | [`03-generated-embedded/`](03-generated-embedded/README.md) |
| 04 | Handwritten model with the CLI wrapper | [`04-handwritten-cli/`](04-handwritten-cli/README.md) |
| 05 | Handwritten model with the embedded solver (in-memory mesh) | [`05-handwritten-embedded/`](05-handwritten-embedded/README.md) |
| 06 | Handwritten model with the embedded solver (distributed mesh) | [`06-handwritten-distributed/`](06-handwritten-distributed/README.md) |

Every section produces the same converged solution `udg`. The
`tutorial/run-all.sh` script runs each section against a chosen
backend variant (CPU, GPU, MPI, MPI+GPU) and verifies the result.

## Prerequisites

Build Exasim per [`../docs/01-installation.md`](../docs/01-installation.md).
Run the tutorial test with:

```bash
ctest --test-dir build -R tutorial
```

Each section's instructions show the build and run commands for
all four backend variants. See the
[overview](00-overview/README.md) for the variant naming
convention.
