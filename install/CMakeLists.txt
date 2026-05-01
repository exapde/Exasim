################################################################################
#                                     EXASIM
#
# Contributing authors: Ngoc Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
################################################################################

cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 70)
endif()

project(exasim LANGUAGES CXX)

get_filename_component(EXASIM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# Exasim's main program is located here
set(EXASIM_MAIN_DIR ${EXASIM_DIR}/backend/Main)
file(GLOB MAIN_SOURCES ${EXASIM_MAIN_DIR}/main.cpp)

# Phase HOT.1: backend/Common/common.cpp was deleted — its contents
# (BLAS sentinels, EXASIM_COMM_* MPI globals) are now C++17 `inline`
# variables in backend/Common/common.h, ODR-safe across multiple TUs
# without a definition TU. The MAIN_SOURCES list is back to just main.cpp.

# Text2Code shared libraries live here
set(Model_LIB_DIR   ${EXASIM_DIR}/backend/Model)

# BuiltinModel shared libraries live here
set(BuiltinModel_LIB_DIR   ${EXASIM_DIR}/backend/Model/BuiltIn)

# ---- Option: build EXASIM variants that link text2code-generated libs ----
option(WITH_TEXT2CODE "Build extra EXASIM executables that link text2code-generated dynamic libraries" OFF)

# logical names for text2code-produced shared libs (lib<name>.so/.dylib/.dll)
set(T2C_CPU_LIB  pdemodelserial)
set(T2C_CUDA_LIB pdemodelcuda)
set(T2C_HIP_LIB  pdemodelhip)

# ---- Option: build EXASIM variants that link builtinmodel libs ----
option(WITH_BUILTINMODEL "Build extra EXASIM executables that link builtin PDE model dynamic libraries" OFF)

# logical names for builtin shared libs (lib<name>.so/.dylib/.dll)
set(BM_CPU_LIB  builtinmodelserial)
set(BM_CUDA_LIB builtinmodelcuda)
set(BM_HIP_LIB  builtinmodelhip)

# ---- Option: link METIS / ParMETIS / GKlib -------------------------------
option(WITH_PARMETIS "Build extra EXASIM executables that link metis libraries" OFF)

# METIS: ParMETIS / METIS / GKlib roots (override if your layout differs)
set(PARMETIS_ROOT "${EXASIM_DIR}/metis/ParMETIS"
    CACHE PATH "Path to ParMETIS root")
set(METIS_ROOT    "${EXASIM_DIR}/metis/METIS"
    CACHE PATH "Path to METIS root")
set(GKLIB_ROOT    "${EXASIM_DIR}/metis/GKlib"
    CACHE PATH "Path to GKlib root")

# METIS: derived include/lib dirs (override if needed)
set(PARMETIS_INCLUDE_DIR "${PARMETIS_ROOT}/include"
    CACHE PATH "ParMETIS include dir")
set(PARMETIS_LIBRARY_DIR "${PARMETIS_ROOT}/lib"
    CACHE PATH "ParMETIS library dir")

# METIS layout: public headers under build/xinclude, library under build/libmetis
set(METIS_INCLUDE_DIR "${METIS_ROOT}/build/xinclude"
    CACHE PATH "METIS public include dir")
set(METIS_LIBRARY_DIR "${METIS_ROOT}/build/libmetis"
    CACHE PATH "METIS library dir")

set(GKLIB_INCLUDE_DIR "${GKLIB_ROOT}/include"
    CACHE PATH "GKlib include dir")
set(GKLIB_LIBRARY_DIR "${GKLIB_ROOT}/lib"
    CACHE PATH "GKlib library dir")

# METIS: helper to attach METIS / ParMETIS / GKlib to a target when WITH_PARMETIS=ON
function(link_metis target)
  if (NOT WITH_PARMETIS)
    return()
  endif()

  # Preprocessor macro so source can guard METIS-specific code
  target_compile_definitions(${target} PRIVATE HAVE_PARMETIS)

  # Includes
  target_include_directories(${target} PRIVATE
    "${PARMETIS_INCLUDE_DIR}"
    "${METIS_INCLUDE_DIR}"
    "${GKLIB_INCLUDE_DIR}"
  )

  # Library search paths (only if they exist, to avoid noisy warnings)
  if (EXISTS "${PARMETIS_LIBRARY_DIR}")
    target_link_directories(${target} PRIVATE "${PARMETIS_LIBRARY_DIR}")
  endif()
  if (EXISTS "${METIS_LIBRARY_DIR}")
    target_link_directories(${target} PRIVATE "${METIS_LIBRARY_DIR}")
  endif()
  if (EXISTS "${GKLIB_LIBRARY_DIR}")
    target_link_directories(${target} PRIVATE "${GKLIB_LIBRARY_DIR}")
  endif()

  # Link order: ParMETIS, METIS, GKlib
  target_link_libraries(${target} PRIVATE
    parmetis
    metis
    GKlib
  )
endfunction()

# BLAS / LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  message(FATAL_ERROR "Lapack and Blas libraries are not found")
endif()

# ----------------------------- CUDA branch ------------------------------
if(EXASIM_CUDA)
  set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildcuda)
  find_package(Kokkos REQUIRED)

  enable_language(CUDA)
  find_package(CUDA REQUIRED)

  include_directories(${CUDA_INCLUDE_DIRS})

  if(EXASIM_NOMPI)
    add_executable(gpuEXASIM ${MAIN_SOURCES})
    target_link_directories(gpuEXASIM PRIVATE ${Model_LIB_DIR})
    target_compile_definitions(gpuEXASIM PRIVATE _CUDA)
    target_compile_options(gpuEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_link_libraries(gpuEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    link_metis(gpuEXASIM)     # METIS

    if(WITH_TEXT2CODE)
      add_executable(gput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _CUDA _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
        ${T2C_CUDA_LIB})
      link_metis(gput2cEXASIM)   # METIS
      set_target_properties(gput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    if(WITH_BUILTINMODEL)
      add_executable(gpubuiltinEXASIM ${MAIN_SOURCES})
      target_link_directories(gpubuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
      target_compile_definitions(gpubuiltinEXASIM PRIVATE _CUDA HAVE_BUILTINMODEL _BUILTINMODEL)
      target_compile_options(gpubuiltinEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpubuiltinEXASIM PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
        ${BM_CUDA_LIB})
      link_metis(gpubuiltinEXASIM)
      set_target_properties(gpubuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_executable(gpumpiEXASIM ${MAIN_SOURCES})
      target_link_directories(gpumpiEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gpumpiEXASIM PRIVATE _MPI _CUDA)
      target_compile_options(gpumpiEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpumpiEXASIM PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES})
      link_metis(gpumpiEXASIM)   # METIS

      if(WITH_TEXT2CODE)
        add_executable(gpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _CUDA _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES}
          ${T2C_CUDA_LIB})
        link_metis(gpumpit2cEXASIM) # METIS
        set_target_properties(gpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()

      if(WITH_BUILTINMODEL)
        add_executable(gpumpibuiltinEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpibuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
        target_compile_definitions(gpumpibuiltinEXASIM PRIVATE _MPI _CUDA HAVE_BUILTINMODEL _BUILTINMODEL)
        target_compile_options(gpumpibuiltinEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpibuiltinEXASIM PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES}
          ${BM_CUDA_LIB})
        link_metis(gpumpibuiltinEXASIM)
        set_target_properties(gpumpibuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()

# ------------------------------ HIP branch ------------------------------
elseif(EXASIM_HIP)
  if(EXASIM_TUOLUMNE)
    set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildtuolumne)
  else()
    set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildhip)
  endif()

  message(${Kokkos_DIR})
  find_package(Kokkos REQUIRED)

  enable_language(HIP)
  find_package(HIP REQUIRED)
  find_package(rocBLAS REQUIRED)
  find_package(hipBLAS REQUIRED)

  include_directories(${HIP_INCLUDE_DIRS})

  if(EXASIM_NOMPI)
    add_executable(gpuEXASIM ${MAIN_SOURCES})
    target_link_directories(gpuEXASIM PRIVATE ${Model_LIB_DIR})
    target_compile_definitions(gpuEXASIM PRIVATE _HIP)
    target_compile_options(gpuEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_link_libraries(gpuEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries})
    link_metis(gpuEXASIM)     # METIS

    if(WITH_TEXT2CODE)
      add_executable(gput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _HIP _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${T2C_HIP_LIB})
      link_metis(gput2cEXASIM)   # METIS
      set_target_properties(gput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    if(WITH_BUILTINMODEL)
      add_executable(gpubuiltinEXASIM ${MAIN_SOURCES})
      target_link_directories(gpubuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
      target_compile_definitions(gpubuiltinEXASIM PRIVATE _HIP HAVE_BUILTINMODEL _BUILTINMODEL)
      target_compile_options(gpubuiltinEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpubuiltinEXASIM PRIVATE
        Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries}
        ${BM_HIP_LIB})
      link_metis(gpubuiltinEXASIM)
      set_target_properties(gpubuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_executable(gpumpiEXASIM ${MAIN_SOURCES})
      target_link_directories(gpumpiEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gpumpiEXASIM PRIVATE _MPI _HIP)
      target_compile_options(gpumpiEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpumpiEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES})
      link_metis(gpumpiEXASIM)   # METIS

      if(WITH_TEXT2CODE)
        add_executable(gpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _HIP _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${T2C_HIP_LIB})
        link_metis(gpumpit2cEXASIM) # METIS
        set_target_properties(gpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()

      if(WITH_BUILTINMODEL)
        add_executable(gpumpibuiltinEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpibuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
        target_compile_definitions(gpumpibuiltinEXASIM PRIVATE _MPI _HIP HAVE_BUILTINMODEL _BUILTINMODEL)
        target_compile_options(gpumpibuiltinEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpibuiltinEXASIM PRIVATE
          Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${BM_HIP_LIB})
        link_metis(gpumpibuiltinEXASIM)
        set_target_properties(gpumpibuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()

# ------------------------------ CPU branch ------------------------------
else()
  set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildserial)
  find_package(Kokkos REQUIRED)

  if(EXASIM_NOMPI)
    add_executable(cpuEXASIM ${MAIN_SOURCES})
    target_link_directories(cpuEXASIM PRIVATE ${Model_LIB_DIR})
    target_link_libraries(cpuEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries})
    target_compile_options(cpuEXASIM PRIVATE -std=c++17 -w -O2)
    link_metis(cpuEXASIM)    # METIS
    #target_compile_options(cpuEXASIM PRIVATE -std=c++17 -fsanitize=address,undefined -fno-omit-frame-pointer -g)
    #target_link_options(cpuEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)

    if(WITH_TEXT2CODE)
      add_executable(cput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(cput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(cput2cEXASIM PRIVATE _TEXT2CODE)
      target_link_libraries(cput2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
      target_compile_options(cput2cEXASIM PRIVATE -std=c++17 -w -O2)   
      link_metis(cput2cEXASIM)   # METIS
      set_target_properties(cput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    if(WITH_BUILTINMODEL)
      add_executable(cpubuiltinEXASIM ${MAIN_SOURCES})
      target_link_directories(cpubuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
      target_compile_definitions(cpubuiltinEXASIM PRIVATE HAVE_BUILTINMODEL _BUILTINMODEL)
      target_link_libraries(cpubuiltinEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${BM_CPU_LIB})
      target_compile_options(cpubuiltinEXASIM PRIVATE -std=c++17 -w -O2)
      link_metis(cpubuiltinEXASIM)
      set_target_properties(cpubuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(${MPI_INCLUDE_PATH})

      add_executable(cpumpiEXASIM ${MAIN_SOURCES})
      target_link_directories(cpumpiEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(cpumpiEXASIM PRIVATE _MPI)
      target_compile_options(cpumpiEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
      #target_compile_options(cpumpiEXASIM PRIVATE -std=c++17 -w -ffast-math -fsanitize=address,undefined -fno-omit-frame-pointer -g)
      #target_link_options(cpumpiEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)  
      target_link_libraries(cpumpiEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
      link_metis(cpumpiEXASIM)   # METIS

      if(WITH_TEXT2CODE)
        add_executable(cpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(cpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(cpumpit2cEXASIM PRIVATE _MPI _TEXT2CODE)
        target_compile_options(cpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
        target_link_libraries(cpumpit2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_CPU_LIB})
        link_metis(cpumpit2cEXASIM) # METIS
        set_target_properties(cpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()

      if(WITH_BUILTINMODEL)
        add_executable(cpumpibuiltinEXASIM ${MAIN_SOURCES})
        target_link_directories(cpumpibuiltinEXASIM PRIVATE ${BuiltinModel_LIB_DIR})
        target_compile_definitions(cpumpibuiltinEXASIM PRIVATE _MPI HAVE_BUILTINMODEL _BUILTINMODEL)
        target_compile_options(cpumpibuiltinEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
        target_link_libraries(cpumpibuiltinEXASIM PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${BM_CPU_LIB})
        link_metis(cpumpibuiltinEXASIM)
        set_target_properties(cpumpibuiltinEXASIM PROPERTIES BUILD_RPATH "${BuiltinModel_LIB_DIR}")
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()
endif()

# ======================================================================
#  HOT.5 — Header-only library install + find_package support
# ======================================================================
#
# Downstream consumers do `find_package(Exasim REQUIRED)` and link
# `Exasim::headers`. The target is INTERFACE-only — no compiled
# library is installed; only the include tree under `include/exasim/`
# (the templated kernel and driver headers) plus the backend headers
# the kernels include.
#
# Kokkos/MPI/BLAS are NOT propagated as `find_dependency` here on
# purpose. Downstream projects choose their own Kokkos build dir,
# their own MPI, and their own BLAS — Exasim only requires they're
# findable, and when they aren't, the user's link step says which
# symbol is missing. Adding `find_dependency(Kokkos)` here would lock
# consumers into the same Kokkos build path Exasim was built against.

add_library(exasim_headers INTERFACE)
add_library(Exasim::headers ALIAS exasim_headers)

# EXPORT_NAME picks the post-namespace target name in the installed
# config (Exasim::headers, not Exasim::exasim_headers).
set_target_properties(exasim_headers PROPERTIES EXPORT_NAME headers)

target_include_directories(exasim_headers INTERFACE
    $<BUILD_INTERFACE:${EXASIM_DIR}/include>
    $<BUILD_INTERFACE:${EXASIM_DIR}/backend>
    $<INSTALL_INTERFACE:include>)

target_compile_features(exasim_headers INTERFACE cxx_std_17)

# Install the public template headers + backend headers the kernels
# transitively include. Filtered to .h/.hpp/.cpp (the template-impl
# files included from .h headers, not standalone TUs).
install(DIRECTORY ${EXASIM_DIR}/include/exasim/
    DESTINATION include/exasim
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h")

# backend/ headers are #include'd by include/exasim/*.hpp via relative
# paths (e.g., `#include "../common.h"`). Install the .h / .hpp / .cpp
# files (the .cpp files here are template-impl bodies #include'd by .h
# headers, not standalone TUs). Skip backend/Model/* — that's the
# legacy ABI which downstream consumers don't need on the templated path.
install(DIRECTORY ${EXASIM_DIR}/backend/
    DESTINATION include/backend
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "*.cpp"
        REGEX "/Model/" EXCLUDE)

install(TARGETS exasim_headers
    EXPORT  ExasimTargets)

install(EXPORT ExasimTargets
    FILE        ExasimTargets.cmake
    NAMESPACE   Exasim::
    DESTINATION lib/cmake/Exasim)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${EXASIM_DIR}/cmake/ExasimConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfig.cmake
    INSTALL_DESTINATION lib/cmake/Exasim)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfigVersion.cmake
    VERSION       0.1.0
    COMPATIBILITY SameMajorVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfigVersion.cmake
    DESTINATION lib/cmake/Exasim)

# If text2code has already been built (text2code is configured as a
# separate CMake project — see text2code/text2code/CMakeLists.txt),
# pick up the binary and ship it. Users can invoke it as
# `${CMAKE_INSTALL_PREFIX}/bin/text2code <pdeapp.txt>`.
set(_text2code_bin "${EXASIM_DIR}/build/text2code")
if(EXISTS "${_text2code_bin}")
    install(PROGRAMS "${_text2code_bin}" DESTINATION bin)
endif()

# ======================================================================
#  HOT.3 — user-facing examples driven by CSolution<MyModel>
# ======================================================================
#
# Each example is a small directory under apps/library_example/<name>/
# containing my_model.hpp + main.cpp. Built only on the CPU/MPI legs
# (the GPU build matrix is left for a follow-up).

option(EXASIM_BUILD_LIBRARY_EXAMPLES
       "Build user-facing examples that consume Exasim as a header-only template library"
       ON)

if(EXASIM_BUILD_LIBRARY_EXAMPLES)
  # Even though CSolution<Poisson2D> never takes the AbiAdapter path
  # (and therefore never calls ::FluxDriver / KokkosFlux / etc.), the
  # `if constexpr` else-branch in EXASIM_DRIVER_CALL still references
  # those symbols. They have to be findable at compile time (declared
  # via Model/ModelDrivers.cpp) and resolvable at link time (provided
  # by libpdemodelserial / libpdemodelcuda / libpdemodelhip). The
  # `_TEXT2CODE` define + linking the matching `pdemodel*` lib handles
  # both. A future cleanup can elide the AbiAdapter branch entirely
  # when EXASIM_BUILD_LIBRARY_EXAMPLES is the only configuration.
  #
  # Per-config target naming convention:
  #   <name>_codegen      — EXASIM_NOMPI && CPU only
  #   <name>_codegen_gpu  — EXASIM_NOMPI && (EXASIM_CUDA || EXASIM_HIP).
  #                          Single name regardless of which Kokkos
  #                          backend; the user's templated kernels
  #                          are GPU-portable (Kokkos::sin etc.) so
  #                          the same source compiles for both.
  # MPI variants are not built here yet — when added they'll follow
  # the same `_codegen_mpi` / `_codegen_mpi_gpu` shape.

  if(EXASIM_NOMPI AND NOT (EXASIM_CUDA OR EXASIM_HIP))
    add_executable(poisson2d_template
        ${EXASIM_DIR}/apps/library_example/poisson2d/main.cpp)
    target_compile_definitions(poisson2d_template PRIVATE _TEXT2CODE)
    target_compile_options(poisson2d_template PRIVATE -std=c++17 -w -O2)
    target_include_directories(poisson2d_template PRIVATE
        ${EXASIM_DIR}/include)
    target_link_directories(poisson2d_template PRIVATE ${Model_LIB_DIR})
    target_link_libraries(poisson2d_template PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
    link_metis(poisson2d_template)
    set_target_properties(poisson2d_template PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")

    # solve_square — fully embedded example: programmatic mesh +
    # programmatic config (no pdeapp.txt, no pdemodel.txt). Tutorial:
    # doc/tutorial_cpu_gpu_mpi.md §5.
    add_executable(solve_square
        ${EXASIM_DIR}/apps/library_example/solve_square/solve_square.cpp)
    target_compile_definitions(solve_square PRIVATE _TEXT2CODE)
    target_compile_options(solve_square PRIVATE -std=c++17 -w -O2)
    target_include_directories(solve_square PRIVATE
        ${EXASIM_DIR}/include)
    target_link_directories(solve_square PRIVATE ${Model_LIB_DIR})
    target_link_libraries(solve_square PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
    link_metis(solve_square)
    set_target_properties(solve_square PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")

    # solve_square_facade — same problem driven through the
    # high-level ExasimSolver<M> wrapper (HOT.7.5).
    add_executable(solve_square_facade
        ${EXASIM_DIR}/apps/library_example/solve_square/solve_square_facade.cpp)
    target_compile_definitions(solve_square_facade PRIVATE _TEXT2CODE)
    target_compile_options(solve_square_facade PRIVATE -std=c++17 -w -O2)
    target_include_directories(solve_square_facade PRIVATE
        ${EXASIM_DIR}/include)
    target_link_directories(solve_square_facade PRIVATE ${Model_LIB_DIR})
    target_link_libraries(solve_square_facade PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
    link_metis(solve_square_facade)
    set_target_properties(solve_square_facade PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  # solve_square_facade_mpi — façade + multi-rank CPU.
  if(EXASIM_MPI AND MPI_FOUND AND NOT (EXASIM_CUDA OR EXASIM_HIP))
    add_executable(solve_square_facade_mpi
        ${EXASIM_DIR}/apps/library_example/solve_square/solve_square_facade.cpp)
    target_compile_definitions(solve_square_facade_mpi PRIVATE _TEXT2CODE _MPI)
    target_compile_options(solve_square_facade_mpi PRIVATE -std=c++17 -w -O2)
    target_include_directories(solve_square_facade_mpi PRIVATE ${EXASIM_DIR}/include)
    target_link_directories(solve_square_facade_mpi PRIVATE ${Model_LIB_DIR})
    target_link_libraries(solve_square_facade_mpi PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_CPU_LIB})
    link_metis(solve_square_facade_mpi)
    set_target_properties(solve_square_facade_mpi PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  # solve_square_facade_gpu — façade + single-rank CUDA/HIP.
  if(EXASIM_NOMPI AND (EXASIM_CUDA OR EXASIM_HIP))
    add_executable(solve_square_facade_gpu
        ${EXASIM_DIR}/apps/library_example/solve_square/solve_square_facade.cpp)
    target_compile_definitions(solve_square_facade_gpu PRIVATE _TEXT2CODE)
    target_compile_options(solve_square_facade_gpu PRIVATE
        -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_include_directories(solve_square_facade_gpu PRIVATE ${EXASIM_DIR}/include)
    target_link_directories(solve_square_facade_gpu PRIVATE ${Model_LIB_DIR})
    if(EXASIM_CUDA)
      target_compile_definitions(solve_square_facade_gpu PRIVATE _CUDA)
      target_link_libraries(solve_square_facade_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries}
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
    else()
      target_compile_definitions(solve_square_facade_gpu PRIVATE _HIP)
      target_link_libraries(solve_square_facade_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${T2C_HIP_LIB})
    endif()
    link_metis(solve_square_facade_gpu)
    set_target_properties(solve_square_facade_gpu PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  # solve_square_facade_mpi_gpu — façade + multi-rank CUDA/HIP.
  if(EXASIM_MPI AND MPI_FOUND AND (EXASIM_CUDA OR EXASIM_HIP))
    add_executable(solve_square_facade_mpi_gpu
        ${EXASIM_DIR}/apps/library_example/solve_square/solve_square_facade.cpp)
    target_compile_definitions(solve_square_facade_mpi_gpu PRIVATE _TEXT2CODE _MPI)
    target_compile_options(solve_square_facade_mpi_gpu PRIVATE
        -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_include_directories(solve_square_facade_mpi_gpu PRIVATE ${EXASIM_DIR}/include)
    target_link_directories(solve_square_facade_mpi_gpu PRIVATE ${Model_LIB_DIR})
    if(EXASIM_CUDA)
      target_compile_definitions(solve_square_facade_mpi_gpu PRIVATE _CUDA)
      target_link_libraries(solve_square_facade_mpi_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
    else()
      target_compile_definitions(solve_square_facade_mpi_gpu PRIVATE _HIP)
      target_link_libraries(solve_square_facade_mpi_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_HIP_LIB})
    endif()
    link_metis(solve_square_facade_mpi_gpu)
    set_target_properties(solve_square_facade_mpi_gpu PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  # Hand-written Poisson 2D, GPU variant. One target name regardless
  # of CUDA vs HIP — Kokkos picks the backend, the templated kernels
  # are portable.
  if(EXASIM_NOMPI AND (EXASIM_CUDA OR EXASIM_HIP))
    add_executable(poisson2d_template_gpu
        ${EXASIM_DIR}/apps/library_example/poisson2d/main.cpp)
    target_compile_definitions(poisson2d_template_gpu PRIVATE _TEXT2CODE)
    target_compile_options(poisson2d_template_gpu PRIVATE
        -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_include_directories(poisson2d_template_gpu PRIVATE
        ${EXASIM_DIR}/include)
    target_link_directories(poisson2d_template_gpu PRIVATE ${Model_LIB_DIR})
    if(EXASIM_CUDA)
      target_compile_definitions(poisson2d_template_gpu PRIVATE _CUDA)
      target_link_libraries(poisson2d_template_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries}
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
    else()
      target_compile_definitions(poisson2d_template_gpu PRIVATE _HIP)
      target_link_libraries(poisson2d_template_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${T2C_HIP_LIB})
    endif()
    link_metis(poisson2d_template_gpu)
    set_target_properties(poisson2d_template_gpu PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  if(EXASIM_MPI AND MPI_FOUND)
    add_executable(poisson2d_template_mpi
        ${EXASIM_DIR}/apps/library_example/poisson2d/main.cpp)
    target_compile_definitions(poisson2d_template_mpi PRIVATE _MPI _TEXT2CODE)
    target_compile_options(poisson2d_template_mpi PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
    target_include_directories(poisson2d_template_mpi PRIVATE
        ${EXASIM_DIR}/include)
    target_link_directories(poisson2d_template_mpi PRIVATE ${Model_LIB_DIR})
    target_link_libraries(poisson2d_template_mpi PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_CPU_LIB})
    link_metis(poisson2d_template_mpi)
    set_target_properties(poisson2d_template_mpi PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
  endif()

  # HOT.4 / HOT.4b: codegen examples — `my_model.hpp` produced by
  # retargeted text2code (not hand-written). Each gates against
  # the corresponding legacy baseline by `outqoi.txt` match.
  function(add_codegen_example name)
    set(_main "${EXASIM_DIR}/apps/library_example/${name}_codegen/main.cpp")
    # Skip apps that haven't been ported yet — `port_codegen.sh <name>`
    # creates the directory + main.cpp on demand.
    if(NOT EXISTS "${_main}")
      message(STATUS "Skipping ${name}_codegen — main.cpp not present (run apps/library_example/port_codegen.sh ${name})")
      return()
    endif()

    # CPU variant (single-rank, no GPU).
    if(EXASIM_NOMPI AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen "${_main}")
      target_compile_definitions(${name}_codegen PRIVATE _TEXT2CODE)
      target_compile_options(${name}_codegen PRIVATE -std=c++17 -w -O2)
      target_include_directories(${name}_codegen PRIVATE
          ${EXASIM_DIR}/include)
      target_link_directories(${name}_codegen PRIVATE ${Model_LIB_DIR})
      target_link_libraries(${name}_codegen PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
      link_metis(${name}_codegen)
      set_target_properties(${name}_codegen PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    # MPI variant (multi-rank CPU).
    if(EXASIM_MPI AND MPI_FOUND AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen_mpi "${_main}")
      target_compile_definitions(${name}_codegen_mpi PRIVATE _MPI _TEXT2CODE)
      target_compile_options(${name}_codegen_mpi PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
      target_include_directories(${name}_codegen_mpi PRIVATE
          ${EXASIM_DIR}/include)
      target_link_directories(${name}_codegen_mpi PRIVATE ${Model_LIB_DIR})
      target_link_libraries(${name}_codegen_mpi PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_CPU_LIB})
      link_metis(${name}_codegen_mpi)
      set_target_properties(${name}_codegen_mpi PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    # GPU variant — one target name covers both CUDA and HIP. Which
    # backend Kokkos uses is decided by how kokkos/buildcuda or
    # kokkos/buildhip was configured; the user's pointwise math is
    # GPU-portable (Kokkos:: qualified by HOT.6.1) so the same TU
    # compiles either way. The `_CUDA` / `_HIP` define routes the
    # device-selection branches in <exasim/run.hpp>.
    if(EXASIM_NOMPI AND (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen_gpu "${_main}")
      target_compile_definitions(${name}_codegen_gpu PRIVATE _TEXT2CODE)
      target_compile_options(${name}_codegen_gpu PRIVATE
          -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_include_directories(${name}_codegen_gpu PRIVATE
          ${EXASIM_DIR}/include)
      target_link_directories(${name}_codegen_gpu PRIVATE ${Model_LIB_DIR})
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_codegen_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_codegen_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
      else()
        target_compile_definitions(${name}_codegen_gpu PRIVATE _HIP)
        target_link_libraries(${name}_codegen_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${T2C_HIP_LIB})
      endif()
      link_metis(${name}_codegen_gpu)
      set_target_properties(${name}_codegen_gpu PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    # MPI + GPU variant — multi-rank GPU. Each MPI rank binds one GPU
    # (the device-selection branch in <exasim/run.hpp> picks
    # `shmrank % deviceCount`).
    if(EXASIM_MPI AND MPI_FOUND AND (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen_mpi_gpu "${_main}")
      target_compile_definitions(${name}_codegen_mpi_gpu PRIVATE _MPI _TEXT2CODE)
      target_compile_options(${name}_codegen_mpi_gpu PRIVATE
          -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_include_directories(${name}_codegen_mpi_gpu PRIVATE
          ${EXASIM_DIR}/include)
      target_link_directories(${name}_codegen_mpi_gpu PRIVATE ${Model_LIB_DIR})
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_codegen_mpi_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_codegen_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
      else()
        target_compile_definitions(${name}_codegen_mpi_gpu PRIVATE _HIP)
        target_link_libraries(${name}_codegen_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_HIP_LIB})
      endif()
      link_metis(${name}_codegen_mpi_gpu)
      set_target_properties(${name}_codegen_mpi_gpu PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()
  endfunction()

  add_codegen_example(poisson2d)
  add_codegen_example(poisson3d)
  add_codegen_example(periodic)
  add_codegen_example(naca0012steady)
  add_codegen_example(lshape)
  add_codegen_example(isoq3d)
  add_codegen_example(cone)
  add_codegen_example(orion)
  add_codegen_example(nsmach8)
  add_codegen_example(isoq)
  add_codegen_example(sharpb2)
  add_codegen_example(naca0012unsteady)
  # Apps not supported by the current codegen contract:
  #   - poisson1d: pdemodel.txt declares `vectors x(2)` but runtime has
  #     nd=1 — codegen bakes nd=2 from pdemodel, runtime asserts ncu/nd
  #     mismatch and aborts. (Latent bug in the source pdemodel — fix
  #     would be to read nd from grid.bin, not vector x size.)
  #   - reactingsharpb2: builtinmodelID=4 routes the runtime to a
  #     chemistry path that needs Mutation++ (HAVE_MPP) compiled in.
  #     Even legacy cput2cEXASIM doesn't converge on this app without
  #     Mutation++, so the codegen path can't either. Out of scope.
endif()

# ======================================================================
#  Summary of EXASIM configuration
# ======================================================================
message(STATUS "==================== EXASIM CONFIGURATION ====================")

message(STATUS "EXASIM_DIR                = ${EXASIM_DIR}")
message(STATUS "EXASIM_CUDA               = ${EXASIM_CUDA}")
message(STATUS "EXASIM_HIP                = ${EXASIM_HIP}")
message(STATUS "EXASIM_MPI                = ${EXASIM_MPI}")
message(STATUS "EXASIM_NOMPI              = ${EXASIM_NOMPI}")
message(STATUS "EXASIM_TUOLUMNE           = ${EXASIM_TUOLUMNE}")

message(STATUS "WITH_TEXT2CODE            = ${WITH_TEXT2CODE}")
message(STATUS "WITH_BUILTINMODEL         = ${WITH_BUILTINMODEL}")
message(STATUS "WITH_PARMETIS             = ${WITH_PARMETIS}")

message(STATUS "Model_LIB_DIR             = ${Model_LIB_DIR}")
message(STATUS "BuiltinModel_LIB_DIR      = ${BuiltinModel_LIB_DIR}")

message(STATUS "PARMETIS_ROOT             = ${PARMETIS_ROOT}")
message(STATUS "METIS_ROOT                = ${METIS_ROOT}")
message(STATUS "GKLIB_ROOT                = ${GKLIB_ROOT}")

message(STATUS "PARMETIS_INCLUDE_DIR      = ${PARMETIS_INCLUDE_DIR}")
message(STATUS "PARMETIS_LIBRARY_DIR      = ${PARMETIS_LIBRARY_DIR}")

message(STATUS "METIS_INCLUDE_DIR         = ${METIS_INCLUDE_DIR}")
message(STATUS "METIS_LIBRARY_DIR         = ${METIS_LIBRARY_DIR}")

message(STATUS "GKLIB_INCLUDE_DIR         = ${GKLIB_INCLUDE_DIR}")
message(STATUS "GKLIB_LIBRARY_DIR         = ${GKLIB_LIBRARY_DIR}")

message(STATUS "BLAS_LIBRARIES            = ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LIBRARIES          = ${LAPACK_LIBRARIES}")

message(STATUS "T2C_CPU_LIB               = ${T2C_CPU_LIB}")
message(STATUS "T2C_CUDA_LIB              = ${T2C_CUDA_LIB}")
message(STATUS "T2C_HIP_LIB               = ${T2C_HIP_LIB}")

message(STATUS "BM_CPU_LIB                = ${BM_CPU_LIB}")
message(STATUS "BM_CUDA_LIB               = ${BM_CUDA_LIB}")
message(STATUS "BM_HIP_LIB                = ${BM_HIP_LIB}")

message(STATUS "CMAKE_CUDA_ARCHITECTURES  = ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CMAKE_BUILD_TYPE          = ${CMAKE_BUILD_TYPE}")

# ---------------- C++ Compiler Info --------------------
message(STATUS "C++ COMPILER              = ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ COMPILER ID           = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ COMPILER VERSION      = ${CMAKE_CXX_COMPILER_VERSION}")

# ---------------- MPI Info -----------------------------
if(MPI_FOUND)
  message(STATUS "MPI FOUND                 = YES")
  message(STATUS "MPI CXX COMPILER          = ${MPI_CXX_COMPILER}")
  message(STATUS "MPI INCLUDE PATH          = ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "MPI LIBRARIES             = ${MPI_CXX_LIBRARIES}")
  message(STATUS "MPIEXEC_EXECUTABLE        = ${MPIEXEC_EXECUTABLE}")
else()
  message(STATUS "MPI FOUND                 = NO")
endif()

# ---------------- CUDA Info ----------------------------
if(EXASIM_CUDA)
  message(STATUS "CUDA FOUND                = YES")
  message(STATUS "CUDA COMPILER             = ${CMAKE_CUDA_COMPILER}")
  message(STATUS "CUDA TOOLKIT ROOT         = ${CUDAToolkit_ROOT}")
  message(STATUS "CUDA LIBRARIES            = ${CUDA_LIBRARIES}")
  message(STATUS "CUDA CUBLAS LIBRARIES     = ${CUDA_CUBLAS_LIBRARIES}")
else()
  message(STATUS "CUDA FOUND                = NO")
endif()

# ---------------- HIP Info -----------------------------
if(EXASIM_HIP)
  message(STATUS "HIP FOUND                 = YES")
  message(STATUS "HIP COMPILER              = ${HIP_HIPCC_EXECUTABLE}")
  message(STATUS "HIP ROOT DIRECTORY        = ${HIP_ROOT_DIR}")
  message(STATUS "HIP INCLUDE DIRS          = ${HIP_INCLUDE_DIRS}")
  message(STATUS "HIP LIBRARIES             = ${HIP_LIBRARIES}")

  # rocBLAS + hipBLAS info if found
  if(rocBLAS_FOUND)
    message(STATUS "rocBLAS FOUND             = YES")
    message(STATUS "rocBLAS LIBRARIES         = ${rocBLAS_LIBRARIES}")
  else()
    message(STATUS "rocBLAS FOUND             = NO")
  endif()

  if(hipBLAS_FOUND)
    message(STATUS "hipBLAS FOUND             = YES")
    message(STATUS "hipBLAS LIBRARIES         = ${hipBLAS_LIBRARIES}")
  else()
    message(STATUS "hipBLAS FOUND             = NO")
  endif()

else()
  message(STATUS "HIP FOUND                 = NO")
endif()

# ---------------- Kokkos Info --------------------------
if(TARGET Kokkos::kokkos)
  message(STATUS "Kokkos TARGET             = Kokkos::kokkos (FOUND)")
else()
  message(STATUS "Kokkos TARGET             = NOT FOUND")
endif()

message(STATUS "Kokkos_DIR                = ${Kokkos_DIR}")

# Some Kokkos packages export a version variable
if(DEFINED Kokkos_VERSION)
  message(STATUS "Kokkos_VERSION            = ${Kokkos_VERSION}")
elseif(DEFINED Kokkos_VERSION_STRING)
  message(STATUS "Kokkos_VERSION            = ${Kokkos_VERSION_STRING}")
endif()

# If Kokkos exported backend info from its config, show it
if(DEFINED Kokkos_ENABLE_SERIAL)
  message(STATUS "Kokkos_ENABLE_SERIAL      = ${Kokkos_ENABLE_SERIAL}")
endif()
if(DEFINED Kokkos_ENABLE_OPENMP)
  message(STATUS "Kokkos_ENABLE_OPENMP      = ${Kokkos_ENABLE_OPENMP}")
endif()
if(DEFINED Kokkos_ENABLE_CUDA)
  message(STATUS "Kokkos_ENABLE_CUDA        = ${Kokkos_ENABLE_CUDA}")
endif()
if(DEFINED Kokkos_ENABLE_HIP)
  message(STATUS "Kokkos_ENABLE_HIP         = ${Kokkos_ENABLE_HIP}")
endif()

message(STATUS "===============================================================")