# Standalone consumer-style CMakeLists.txt — what an external,
# multi-file downstream user would write. Used by both modes:
#
#   In-tree (driven from ../CMakeLists.txt under
#   EXASIM_BUILD_TESTS=ON): not consumed; the parent file registers
#   the target directly.
#
#   Out-of-tree (post-install):
#     cmake --install $EXASIM/build --prefix /tmp/exasim_install
#     cmake -S . -B /tmp/multi_tu_build \
#           -DCMAKE_PREFIX_PATH=/tmp/exasim_install
#     cmake --build /tmp/multi_tu_build
#
# The link step IS the test. Any non-inline free function defined in
# a public Exasim header pulls into both unit_a.o and unit_b.o; the
# linker reports duplicate symbols. PR #73 review B2 expects this
# test to pass after the ODR sweep.

cmake_minimum_required(VERSION 3.16)
project(multi_tu_consumer CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Exasim REQUIRED)
find_package(Kokkos REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(multi_tu_consumer main.cpp unit_a.cpp unit_b.cpp)
target_compile_definitions(multi_tu_consumer PRIVATE _TEXT2CODE)
target_link_libraries(multi_tu_consumer PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
