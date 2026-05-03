# Consumer-style CMakeLists for the "generated model + CLI wrapper"
# tutorial. This is what your CMakeLists.txt would look like if you
# wrote a real out-of-tree consumer of Exasim.
#
# This file is NOT used by the in-tree tutorial build (which
# registers targets from `tutorial/CMakeLists.txt` directly). Use
# it as a standalone build entry point:
#
#   # 1. Install Exasim once.
#   cmake --install $EXASIM/build --prefix /opt/exasim
#
#   # 2. Run text2code in this directory to produce
#   #    my_model.hpp + libpdemodelserial.{so,dylib}.
#   $EXASIM/build/text2code --out-dir . ./pdeapp.txt
#
#   # 3. Build this app.
#   cmake -S . -B build \
#         -DCMAKE_PREFIX_PATH=/opt/exasim \
#         -DKokkos_DIR=/opt/exasim/external/kokkos
#   cmake --build build
#
#   # 4. Run.
#   cd /path/to/this/dir
#   mkdir -p datain dataout
#   ./build/tutorial_02_generated_cli ./pdeapp.txt

cmake_minimum_required(VERSION 3.16)
project(tutorial_02_generated_cli CXX)
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

# The codegen library `libpdemodelserial.{so,dylib}` is emitted by
# text2code into this directory. Link against it.
target_link_directories(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${PROJECT_NAME} PRIVATE pdemodelserial)
set_target_properties(${PROJECT_NAME} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
