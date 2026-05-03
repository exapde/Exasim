# Tutorial section 03 — generated model with the embedded solver.
#
# Builds `tutorial_03_generated_embedded` from `main.cpp`. The build
# requires `my_model.hpp` and `libpdemodelserial.{so,dylib}` to
# exist in this directory. The script `tutorial/run-all.sh`
# generates them with `text2code --out-dir <here>` before invoking
# `cmake --build`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_03_generated_embedded)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

target_link_directories(${_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
