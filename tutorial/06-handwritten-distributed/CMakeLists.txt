# Tutorial section 06 — handwritten model with the distributed-mesh
# embedded solver.
#
# Builds `tutorial_06_handwritten_distributed`, which is meaningful
# only on MPI-enabled build dirs (`build_mpi`, `build_mpi_gpu`).
# On serial build dirs the binary still compiles and runs as a
# single rank with mpiprocs=1, but the path is exercised primarily
# under `mpirun -np N`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_06_handwritten_distributed)

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
