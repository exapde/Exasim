# Tutorial section 05 — handwritten model with the embedded solver.
#
# Builds `tutorial_05_handwritten_embedded`. The binary builds its
# own mesh as flat C arrays, drives the solver through the
# `ExasimSolver<Poisson2D>` API, and prints the maximum absolute
# value of the converged `udg`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_05_handwritten_embedded)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")
