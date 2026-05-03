# Common build steps

Steps that are the same on every platform. The per-platform pages
list system packages, link to this page for the dependency builds,
then list the Exasim configure line.

## Sources

| Library | Vendored at | Upstream |
|---|---|---|
| Kokkos | `kokkos/` | <https://github.com/kokkos/kokkos> |
| GKlib | `metis/GKlib/` | <https://github.com/KarypisLab/GKlib> |
| METIS | `metis/METIS/` | <https://github.com/KarypisLab/METIS> |
| ParMETIS | `metis/ParMETIS/` | <https://github.com/KarypisLab/ParMETIS> |
| SymEngine | `text2code/symengine/` | <https://github.com/symengine/symengine> |

Vendored copies are already in the repo. The upstream URLs are for
reference only.

## Build GKlib + METIS + ParMETIS

```bash
EXASIM=$PWD                    # repo root

cd $EXASIM/metis/GKlib
make config prefix=$EXASIM/metis/GKlib
make -j8 install

cd $EXASIM/metis/METIS
make config gklib_path=$EXASIM/metis/GKlib \
            prefix=$EXASIM/metis/METIS \
            shared=1 r64=1 i64=1
make -j8 install
```

ParMETIS is required only for MPI builds:

```bash
cd $EXASIM/metis/ParMETIS
PATH=$MPI_PREFIX/bin:$PATH make config \
    gklib_path=$EXASIM/metis/GKlib \
    metis_path=$EXASIM/metis/METIS \
    prefix=$EXASIM/metis/ParMETIS \
    cc=mpicc cxx=mpicxx
PATH=$MPI_PREFIX/bin:$PATH make -j8 install
```

`$MPI_PREFIX` is the MPI install prefix from the per-platform page.

## Build Kokkos serial

Required on every platform. The CPU and MPI Exasim builds link against
this prefix.

```bash
cd $EXASIM/kokkos
mkdir buildserial && cd buildserial
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/kokkos/buildserial \
         -DKokkos_ENABLE_SERIAL=ON
make -j16 install
```

## Build SymEngine + text2code

Required for the codegen authoring path
([`docs/03-methods/codegen-text2code.md`](../03-methods/codegen-text2code.md)).
Skip if you only use the hand-written Model path.

```bash
cd $EXASIM/text2code/symengine
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$EXASIM/text2code/symengine/build \
         -DBUILD_SHARED_LIBS=OFF -DWITH_GMP=OFF -DINTEGER_CLASS=boostmp
make -j16 install

cd $EXASIM/text2code/text2code
cmake -S . -B build -DEXASIM_ROOT=$EXASIM
cmake --build build -j8
```

The text2code binary lands at `$EXASIM/build/text2code`.

## System-wide install

After building any backend variant:

```bash
cmake --install build_cpu --prefix /opt/exasim
```

Layout:

```
/opt/exasim/
├── include/exasim/          # public headers
├── include/backend/         # template-impl headers transitively included
├── lib/cmake/Exasim/        # find_package(Exasim) machinery
└── bin/text2code            # codegen tool (if text2code was built)
```

Consumer projects do `find_package(Exasim REQUIRED)` and link
`Exasim::headers Kokkos::kokkos MPI::MPI_CXX`. See
[`03-methods/hand-written-model.md`](../03-methods/hand-written-model.md)
for the consumer CMake template.
