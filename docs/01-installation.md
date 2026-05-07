# Installation

Pick the page that matches your platform.

| Platform | Page |
|---|---|
| macOS (Apple Silicon or Intel) | [`01-installation/macos.md`](01-installation/macos.md) |
| Linux x86_64, CPU only | [`01-installation/linux-cpu.md`](01-installation/linux-cpu.md) |
| Linux x86_64, NVIDIA GPU | [`01-installation/linux-nvidia.md`](01-installation/linux-nvidia.md) |
| Linux x86_64, AMD GPU | [`01-installation/linux-amd.md`](01-installation/linux-amd.md) |

All platforms share the dependency build steps in
[`01-installation/common.md`](01-installation/common.md). The per-platform
pages reference it.

## Repository folder name

The cloned folder must be named `Exasim`. Renaming it breaks the
Kokkos build paths.
