# StaticLibraryABI

This is a minimal C++ example of:

- one precompiled core static library
- one uniform runtime provider ABI
- three providers for the same operation

The operation is:

```cpp
double Add(double a, double b)
{
    return a + b;
}
```

Provider variants:

1. `providers/dynamic/dynamic_provider.cpp`
   - implements `Add` in a shared library
2. `providers/source/source_provider.cpp`
   - implements `Add` in a plain C++ source file
3. `providers/template/template_provider.cpp`
   - implements:
     - `template <class Model> double AddTemplate(double a, double b)`
     - `double Add(double a, double b) { return AddTemplate<MyModel>(a, b); }`

Uniform ABI:

- `include/driver_abi.hpp`
- every provider exports:

```cpp
extern "C" const staticabi::DriverABI* GetDriverABI();
```

The core static library:

- lives in `src/core.cpp`
- knows only `DriverABI`
- never knows how the provider was implemented

Dynamic-library path handling:

- `run_dynamic` links directly against `dynamic_provider`
- CMake embeds the provider directory into the executable build rpath
- no explicit `dlopen` loader class is used

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run

```bash
./build/run_source
./build/run_template
./build/run_dynamic
```

On Linux, the shared library will be `.so`, but `run_dynamic` is still launched
the same way because the path is handled by the executable link settings.
