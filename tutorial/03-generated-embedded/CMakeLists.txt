# Consumer-style CMakeLists for the "generated model + embedded
# solver" tutorial. This is what your CMakeLists.txt would look
# like if you wrote a real out-of-tree consumer of Exasim.
#
# This file is NOT used by the in-tree tutorial build (which
# registers targets from `tutorial/CMakeLists.txt` directly). Use
# it as a standalone build entry point:
#
#   cmake --install $EXASIM/build --prefix /opt/exasim
#   $EXASIM/build/text2code --out-dir . ./pdeapp.txt
#   cmake -S . -B build \
#         -DCMAKE_PREFIX_PATH=/opt/exasim \
#         -DKokkos_DIR=/opt/exasim/external/kokkos
#   cmake --build build
#   ./build/tutorial_03_generated_embedded ./pdeapp.txt

cmake_minimum_required(VERSION 3.16)
project(tutorial_03_generated_embedded CXX)
set(CMAKE_CXX_STANDARD 17)

find_package(Exasim REQUIRED)
find_package(Kokkos REQUIRED)
find_package(BLAS   REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)
target_compile_definitions(${PROJECT_NAME} PRIVATE _TEXT2CODE)
target_link_libraries(${PROJECT_NAME} PRIVATE
    Exasim::headers Kokkos::kokkos
    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_link_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${PROJECT_NAME} PRIVATE pdemodelserial)
set_target_properties(${PROJECT_NAME} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
