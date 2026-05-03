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
tutorial_configure_target(${_target})

target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")
