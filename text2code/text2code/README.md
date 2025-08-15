# text2code

`text2code` is a C++17 utility that **parses a PDE application description** (from `pdeapp.txt` and companion model files like `pdemodel.txt`), **builds mesh/master/connectivity data**, optionally **generates C++ code for model kernels**, and **writes binary inputs** that can be consumed by downstream solvers (e.g., Exasim backends).

The project is self-contained (no external parser framework) and uses a lightweight math expression engine (`tinyexpr`) to evaluate formulas found in the text configuration.

---

## Highlights

- **Plain-text → binary**: read `pdeapp.txt` / `pdemodel.txt` and emit compact binaries (e.g., `mesh.bin`, etc.) in `datain path`.
- **Mesh/master utilities**: routines for reading/building meshes (`readmesh.cpp`, `makemesh.cpp`), creating master data (`makemaster.cpp`), connectivity (`connectivity.cpp`), and domain decomposition (`domaindecomposition.cpp`).
- **Code generation hooks**: when `pde.gencode` is enabled, the generator (`CodeGenerator.*`) prepares model code; a companion compiler path (`CodeCompiler.cpp`) can build helper objects or dynamic libraries (depending on your CMake option `USE_CMAKE`).
- **Portable, minimal deps**: standard C++17 and CMake. Optional partitioning hooks can use METIS if available (guarded in sources).

---

## Project structure

```
text2code/
|-- CMakeLists.txt
|-- text2code.cpp            # main entry: parses pdeapp, orchestrates pipeline
|-- CodeGenerator.cpp/.hpp   # emit C++ for model kernels (enabled by pde.gencode)
|-- CodeCompiler.cpp         # (optional) compile generated sources (USE_CMAKE path)
|-- TextParser.hpp           # tiny text/field parser for key=val style inputs
|-- tinyexpr.cpp/.h          # expression evaluator used during parsing
|-- readpdeapp.cpp           # read pdeapp.txt and friends; populate pde struct
|-- pdemodel.txt             # example model description (symbols, sources, etc.)
|-- pdeapp.txt               # example application configuration
|-- readmesh.cpp             # read mesh files (or generate) into memory
|-- makemesh.cpp             # create/modify mesh (uniform grids, bounds, etc.)
|-- makemaster.cpp           # create master element/quad rules, shapes
|-- connectivity.cpp         # element/face connectivity, maps
|-- domaindecomposition.cpp  # split mesh into parts (METIS hooks)
|-- helpers.cpp              # common utilities (I/O, filesystem helpers, timers)
|-- writebinaryfiles.cpp     # write mesh/master/solution/etc. to *.bin
```

---

## Build

Requirements:

- **CMake 3.16 minimum**
- **A C++17 compiler** (Clang, GCC, or MSVC). On macOS, Xcode command line tools are sufficient.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces an executable (typically named **`text2code`**) under `build/`, which is moved to .

---

## Usage

At minimum, pass your **application file** (`pdeapp.txt`). The program will discover related files (e.g., `pdemodel.txt`) based on paths in `pdeapp.txt`.

```bash
./build/text2code pdeapp.txt
```

### Common knobs (driven by `pdeapp.txt`)

- **`pde.gencode`**: when true, `CodeGenerator` emits model C++ (e.g., flux/source kernels).
- **`USE_CMAKE`** (CMake option or env): select the path used by `CodeCompiler` to build generated code.
- **`datainpath`**: directory where binary outputs are written.

---

## Typical pipeline

1. **Parse app/model** — `readpdeapp.cpp` and `TextParser.hpp` read key-value sections; formulas are evaluated by **tinyexpr**.  
2. **Build geometry/master/connectivity** — `readmesh.cpp`, `makemesh.cpp`, `makemaster.cpp`, `connectivity.cpp` build mesh/master structures.  
3. **(Optional) Codegen** — if `pde.gencode=1`, `CodeGenerator` produces model code.  
4. **(Optional) Compile generated code** — via `CodeCompiler.cpp` if enabled.  
5. **Write binaries** — `writebinaryfiles.cpp` writes `.bin` files into `datainpath`.  

---

## Examples

Minimal example:

```bash
# Build:
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run:
./build/text2code pdeapp.txt
```

---

## Partitioning with METIS

`domaindecomposition.cpp` includes METIS hooks for partitioning. Disable partitioning for single-part runs or install METIS for multi-partition use.

---

## Troubleshooting

- **No outputs**: check `datainpath` in `pdeapp.txt`.  
- **Compile errors in generated code**: disable `pde.gencode` to isolate.  
- **Linker errors for METIS**: disable partitioning or link METIS.  
