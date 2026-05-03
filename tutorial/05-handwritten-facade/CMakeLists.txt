# Tutorial path 02 — hand-written Model + embedded `ExasimSolver<M>`.
#
# Builds `tutorial_05_handwritten_facade`. The binary builds its
# own mesh as flat C arrays, drives the solver through the
# `ExasimSolver<Poisson2D>` API, and prints the maximum absolute
# value of the converged `udg`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_05_handwritten_facade)

add_executable(${_target} main.cpp)
target_compile_features(${_target} PRIVATE cxx_std_17)
target_compile_definitions(${_target} PRIVATE _TEXT2CODE)
target_include_directories(${_target} PRIVATE
    ${EXASIM_DIR}/include
    "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
target_link_libraries(${_target} PRIVATE
    Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
link_metis(${_target})
set_target_properties(${_target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${Model_LIB_DIR}")
