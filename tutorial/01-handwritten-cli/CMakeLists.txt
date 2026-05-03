# Tutorial path 01 — hand-written Model + legacy CLI driver.
#
# Builds the binary `tutorial_01_handwritten_cli` against the
# Exasim header-only library. The binary takes `pdeapp.txt` as its
# single argument and writes its output to `dataout/`.

cmake_minimum_required(VERSION 3.16)

set(_target tutorial_01_handwritten_cli)

add_executable(${_target} main.cpp)
target_compile_features(${_target} PRIVATE cxx_std_17)
target_compile_definitions(${_target} PRIVATE _TEXT2CODE)
target_include_directories(${_target} PRIVATE
    ${EXASIM_DIR}/include
    "${CMAKE_CURRENT_SOURCE_DIR}")

# The target links Kokkos for the device-portable kernels and the
# text2code-emitted model library that the legacy ABI expects to be
# present at link time. Both also feed the templated kernels.
target_link_directories(${_target} PRIVATE ${Model_LIB_DIR})
target_link_libraries(${_target} PRIVATE
    Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
link_metis(${_target})
set_target_properties(${_target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${Model_LIB_DIR}")
