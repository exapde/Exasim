# In-tree wiring for tests/ subprojects.
#
# Each subdirectory has its own consumer-style CMakeLists.txt for
# documentation and out-of-tree build. This file registers the same
# tests directly against the in-tree backend so a developer can run
# them without a separate `cmake --install` step.
#
# Gated by EXASIM_BUILD_TESTS (default OFF). Some tests are expected
# to FAIL until the corresponding fix lands — see comments below.

cmake_minimum_required(VERSION 3.16)

# multi_tu_consumer — surfaces multi-TU ODR violations in the public
# header set. Three .cpp files include <exasim/run.hpp>; if any
# transitive header defines a non-inline free function, the link
# fails. PR #73 review B2 will inline those; before that lands this
# target will not link, which is the point.
#
# Setup mirrors tutorial 04 (hand-written Poisson2D + libpdemodel.cpp
# compiled directly into the binary) so the only difference from a
# known-good single-TU build is the unit_a/unit_b/main split. That
# isolates ODR breakage from any other linkage concern.
set(_mtu_pdemodel_src "${EXASIM_DIR}/backend/Model/libpdemodel.cpp")
add_executable(multi_tu_consumer
    multi_tu_consumer/main.cpp
    multi_tu_consumer/unit_a.cpp
    multi_tu_consumer/unit_b.cpp
    ${_mtu_pdemodel_src})
target_compile_features(multi_tu_consumer PRIVATE cxx_std_17)
target_compile_definitions(multi_tu_consumer PRIVATE _TEXT2CODE)
target_include_directories(multi_tu_consumer PRIVATE
    ${EXASIM_DIR}/include
    ${EXASIM_DIR}/backend/Model
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_tu_consumer")
target_link_libraries(multi_tu_consumer PRIVATE
    Kokkos::kokkos ${lapackblas_libraries})

if(EXASIM_MPI AND MPI_FOUND)
    target_compile_definitions(multi_tu_consumer PRIVATE _MPI)
    target_link_libraries(multi_tu_consumer PRIVATE ${MPI_LIBRARIES})
endif()
if(EXASIM_CUDA)
    target_compile_definitions(multi_tu_consumer PRIVATE _CUDA)
    set_source_files_properties(${_mtu_pdemodel_src} PROPERTIES
        LANGUAGE CUDA
        COMPILE_OPTIONS "--expt-extended-lambda;--expt-relaxed-constexpr")
elseif(EXASIM_HIP)
    target_compile_definitions(multi_tu_consumer PRIVATE _HIP)
    set_source_files_properties(${_mtu_pdemodel_src} PROPERTIES
        LANGUAGE HIP)
endif()

link_metis(multi_tu_consumer)
set_target_properties(multi_tu_consumer PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
