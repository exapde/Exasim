# multi_tu_consumer — multi-TU ODR regression test

Three translation units (`unit_a.cpp`, `unit_b.cpp`, `main.cpp`) each
include `<exasim/run.hpp>`. The link step is the test: any non-inline
free function defined inside the public-header transitive closure
will be present in two or more `.o` files, and the linker will refuse
to merge them.

This is the regression check for PR #73 review blocker B2 (header-only
ODR safety). Until B2 is addressed, the test is **expected to fail**;
the failure surfaces the audit list of symbols that need `inline`.

## Running

### In-tree (the usual case during development)

```
cmake -S install -B build \
      -D EXASIM_NOMPI=ON -D EXASIM_BUILD_TESTS=ON
cmake --build build --target multi_tu_consumer
```

The test passes when the `multi_tu_consumer` binary links cleanly. The
binary itself is not run — link success is the whole signal.

### Out-of-tree (post-install, after B1)

```
cmake --install build --prefix /tmp/exasim_install
cmake -S tests/multi_tu_consumer -B /tmp/multi_tu_build \
      -DCMAKE_PREFIX_PATH=/tmp/exasim_install
cmake --build /tmp/multi_tu_build
```

Out-of-tree mode also exercises B1 (install header layout); both
must be passing for this command to succeed.

## What "broken" looks like today

`ld` reports lines of the form

```
duplicate symbol 'parseExpression(...)' in:
    .../unit_a.cpp.o
    .../unit_b.cpp.o
    .../main.cpp.o
```

Each duplicate is a non-inline free function defined in a header
under `include/exasim/`'s transitive include set. Adding `inline`
(or moving the body to a `.cpp` carried by the static lib, or moving
to an anonymous namespace if internal linkage is acceptable) fixes
the violation.
