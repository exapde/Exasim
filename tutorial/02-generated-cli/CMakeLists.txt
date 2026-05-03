# Tutorial section 02 — generated model with the CLI wrapper.
#
# Builds `tutorial_02_generated_cli` from `main.cpp`. The build
# requires `my_model.hpp` and `libpdemodelserial.{so,dylib}` to
# already exist in this directory; both are produced by running
# `text2code --out-dir <here> ./pdeapp.txt`. The script
# `tutorial/run-all.sh` runs `text2code` before invoking
# `cmake --build`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_02_generated_cli)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

# This section uses the codegen library that text2code wrote into
# its own source directory (not the global Model_LIB_DIR), so
# override the link search path and rpath accordingly.
target_link_directories(${_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
