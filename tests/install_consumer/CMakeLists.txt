# Out-of-tree install consumer — what a real downstream user writes
# against an installed Exasim. Validates PR #73 review blocker B1
# (find_package(Exasim) consumers can resolve the public headers and
# their transitive backend includes).
#
# Reviewer-runnable script:
#
#   cmake --install $EXASIM/build --prefix /tmp/exasim_install
#   cmake -S tests/install_consumer -B /tmp/install_consumer_build \
#         -DCMAKE_PREFIX_PATH=/tmp/exasim_install
#   cmake --build /tmp/install_consumer_build
#
# Build success means find_package(Exasim) found the installed
# package, the public header transitive include chain (B1) resolved
# correctly out-of-tree, and the public-header ODR set (B2) linked
# cleanly.
#
# The consumer uses a hand-written `Poisson2D` model and the
# templated `exasim::run<Poisson2D>(...)` path. We deliberately do
# NOT define `_TEXT2CODE`: that macro routes through the AbiAdapter
# else-branch which expects per-model codegen kernels (KokkosFlux,
# HdgFlux, etc.) provided by libpdemodel.{so,dylib}. Out-of-tree
# without text2code'd kernels available, the templated branch is the
# only buildable path — and exercising it is what the test is for.

cmake_minimum_required(VERSION 3.16)
project(install_consumer CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Exasim REQUIRED)
find_package(Kokkos REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(install_consumer main.cpp)
target_include_directories(install_consumer PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}")           # for #include "my_model.hpp"
target_link_libraries(install_consumer PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
