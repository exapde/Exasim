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
target_compile_features(${_target} PRIVATE cxx_std_17)
target_compile_definitions(${_target} PRIVATE _TEXT2CODE)
target_include_directories(${_target} PRIVATE
    ${EXASIM_DIR}/include
    "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_directories(${_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${_target} PRIVATE
    Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
link_metis(${_target})
set_target_properties(${_target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}")
