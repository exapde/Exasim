# Tutorial section 04 — handwritten model with the CLI wrapper.
#
# Builds the binary `tutorial_04_handwritten_cli` against the
# Exasim header-only library. The binary takes `pdeapp.txt` as its
# single argument and writes its output to `dataout/`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_04_handwritten_cli)

add_executable(${_target} main.cpp)
tutorial_configure_target(${_target})

# The handwritten path links the global Model_LIB_DIR placeholder
# library so the legacy ABI plumbing has its symbols available.
target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
set_target_properties(${_target} PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")
