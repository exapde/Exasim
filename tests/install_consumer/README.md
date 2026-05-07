# install_consumer — out-of-tree install layout regression test

A minimal `find_package(Exasim)` consumer built **outside** the
Exasim source tree, against an installed copy. The build itself is
the test: it validates that

1. `find_package(Exasim REQUIRED)` finds the installed package config
   (`<prefix>/lib/cmake/Exasim/`).
2. `Exasim::headers` exposes both `<exasim/...>` and `<backend/...>`
   include resolutions correctly out of `<prefix>/include/`.
3. The transitive include chain from a public header back into
   `backend/...` resolves entirely against `<prefix>/include/`, not
   against the source tree. (PR #73 review B1.)
4. The public-header ODR set links cleanly. (PR #73 review B2.)

Test scope is intentionally narrow: it pulls `<exasim/model.hpp>`
and runs a `static_assert` on a downstream-defined `Poisson2D`. It
does **not** include `<exasim/run.hpp>` — that drags in
`backend/Model/KokkosDrivers.cpp`, which transitively `#include`s
per-model codegen kernels (`KokkosFlux.cpp`, …) that aren't part of
the install (they're text2code output, generated per model by the
consumer at runtime). A separate end-to-end consumer test that
exercises `exasim::run` is left for after the codegen-vs-install
contract is firmed up.

## Running

```
# 1. Install Exasim (uses the in-tree build dir)
cmake --install build --prefix /tmp/exasim_install

# 2. Configure the consumer against the installed package
cmake -S tests/install_consumer -B /tmp/install_consumer_build \
      -D "CMAKE_PREFIX_PATH=/tmp/exasim_install;/path/to/kokkos/buildserial"

# 3. Build (this is the test)
cmake --build /tmp/install_consumer_build

# 4. Run for sanity (returns 0)
/tmp/install_consumer_build/install_consumer
```

`CMAKE_PREFIX_PATH` needs Kokkos's build/install location alongside
Exasim's prefix (Exasim's installed package config does not vendor
Kokkos — consumers point at the same Kokkos build Exasim was
configured against).

## What "broken" looks like

If B1 regresses, configure or compile fails with errors like:

```
fatal error: 'backend/Common/common.h' file not found
```

(an installed wrapper resolved a `<backend/...>` include against the
source tree's relative layout instead of `<prefix>/include/backend/`),
or:

```
In file included from /Users/.../tests/install_consumer/../../include/exasim/...
error: redefinition of 'AbiAdapter'
```

(a backend header used `#include "../../include/exasim/..."` which
resolved through the consumer's source dir into the original source
tree, double-including a header that was already pulled in via
`<exasim/...>`).
