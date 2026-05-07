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

# Global include path so the public wrappers' angle-form
# `<backend/...>` includes resolve to ${EXASIM_DIR}/backend/. Out-of-
# tree consumers go through Exasim::headers (which sets the same
# resolution against <prefix>/include); this directory-scope rule
# covers in-tree targets that don't link Exasim::headers explicitly
# (cpuEXASIM, the library examples, the tutorial subprojects, etc.).
include_directories(${EXASIM_DIR} ${EXASIM_DIR}/include)

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
# — historically each binary linked one of these as a separately
# loaded .so. We now compile `backend/Model/libpdemodel.cpp` inline
# instead (see `exasim_link_pdemodel()` below). text2code still
# emits the .so files for legacy callers, but in-tree binaries no
# longer load them.
set(T2C_CPU_LIB  pdemodelserial)
set(T2C_CUDA_LIB pdemodelcuda)
set(T2C_HIP_LIB  pdemodelhip)

# Compile `backend/Model/libpdemodel.cpp` directly into `target`,
# replacing the legacy pattern `target_link_libraries(... ${T2C_*_LIB})`
# which used to load text2code's
# `libpdemodel{serial,cuda,hip}.{so,dylib}` as a separately loaded
# shared library.
#
# Why: with the `.so` pattern, Kokkos ends up statically linked into
# both the main exe and the loaded `libpdemodel*.so`, and Kokkos's
# static class members (`SerialInternal::all_instances`,
# `CudaInternal::constantMemHostStagingPerDevice`, ...) get duplicate
# `__cxa_atexit` destructor registrations — one per DSO. At program
# exit both fire and double-free the same heap storage, aborting
# with `free(): double free detected in tcache 2` /
# `Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()`.
# Math is correct; only shutdown is broken.
#
# Compiling `libpdemodel.cpp` as a TU of the binary keeps the
# process single-DSO with respect to Kokkos, so the destructor is
# registered exactly once.
# `libpdemodel.cpp` does `#include "cpuInitu.cpp"` etc — quote-include
# search rules look in the directory of the *including file* first,
# before the -I list. backend/Model/ ships zero-fill stubs (cpuInitu
# returns 0 for all components), so compiling `backend/Model/libpdemodel.cpp`
# always picks up the stubs even when the per-example -I has the real
# codegen output. NS apps with non-zero free-stream initial conditions
# (e.g. naca0012*) then crash with `Residual norm is nan` at iter 0
# because density=0 → 1/r=∞ in Flux. Per-example codegen targets must
# instead compile the per-example libpdemodel.cpp (regenerate.sh writes
# one into each <name>_codegen/ dir) so the relative #includes resolve
# to the real codegen Symbolic*/Kokkos*/cpuInit* files in that dir.
function(exasim_link_pdemodel target)
  if(ARGC GREATER 1)
    set(_pdemodel_src "${ARGV1}")
  else()
    set(_pdemodel_src "${EXASIM_DIR}/backend/Model/libpdemodel.cpp")
  endif()
  target_sources(${target} PRIVATE ${_pdemodel_src})
  target_include_directories(${target} PRIVATE
    "${EXASIM_DIR}/backend/Model")
  if(EXASIM_CUDA)
    set_source_files_properties(${_pdemodel_src} PROPERTIES
      LANGUAGE CUDA
      COMPILE_OPTIONS "--expt-extended-lambda;--expt-relaxed-constexpr")
  elseif(EXASIM_HIP)
    set_source_files_properties(${_pdemodel_src} PROPERTIES
      LANGUAGE HIP)
  endif()
endfunction()

# ---- Option: build EXASIM variants that link builtinmodel libs ----
option(WITH_BUILTINMODEL "Build extra EXASIM executables that link builtin PDE model dynamic libraries" OFF)

# logical names for builtin shared libs (lib<name>.so/.dylib/.dll)
set(BM_CPU_LIB  builtinmodelserial)
set(BM_CUDA_LIB builtinmodelcuda)
set(BM_HIP_LIB  builtinmodelhip)

# ---- Option: AddressSanitizer (ASAN) -------------------------------------
#
# Use to chase host-side memory bugs (use-after-free, double-free, OOB).
# When ON, every Exasim target gets `-fsanitize=address
# -fno-omit-frame-pointer` on host compiles and `-fsanitize=address` at
# link. CUDA TUs also get the flags via `-Xcompiler`.
#
# Recommended runtime envs to pair with this:
#   ASAN_OPTIONS=protect_shadow_gap=0:detect_leaks=0:abort_on_error=1
#   LSAN_OPTIONS=suppressions=/dev/null
#
# `protect_shadow_gap=0` is required for CUDA (cuBLAS/cuRAND etc. mmap
# into the shadow gap region by default and ASAN flags those as bugs).
#
# Caveat: ASAN does not catch errors *inside* CUDA kernels — use
# compute-sanitizer for those.
option(EXASIM_ASAN "Build with AddressSanitizer (-fsanitize=address)" OFF)

if(EXASIM_ASAN)
  message(STATUS "EXASIM_ASAN              = ON  (host-side memory error detection)")
  add_compile_options(
      $<$<COMPILE_LANGUAGE:CXX,C>:-fsanitize=address>
      $<$<COMPILE_LANGUAGE:CXX,C>:-fno-omit-frame-pointer>
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fsanitize=address>
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fno-omit-frame-pointer>)
  add_link_options(-fsanitize=address)
endif()

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
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    exasim_link_pdemodel(gput2cEXASIM)
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
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES})
    exasim_link_pdemodel(gpumpit2cEXASIM)
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
      target_link_libraries(gput2cEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries})
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
          Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES})
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
      target_link_libraries(cput2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries})
    exasim_link_pdemodel(cput2cEXASIM)
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
        target_link_libraries(cpumpit2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
    exasim_link_pdemodel(cpumpit2cEXASIM)
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
    # Public include forms (post-B1):
    #   <exasim/...>  resolves under <prefix>/include/exasim/
    #                 (in-tree: ${EXASIM_DIR}/include/exasim/)
    #   <backend/...> resolves under <prefix>/include/backend/
    #                 (in-tree: ${EXASIM_DIR}/backend/)
    # In-tree we expose ${EXASIM_DIR} itself so <backend/...> hits the
    # source tree directly; installed we expose <prefix>/include so it
    # hits the installed copy under include/backend/.
    $<BUILD_INTERFACE:${EXASIM_DIR}/include>
    $<BUILD_INTERFACE:${EXASIM_DIR}>
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

# backend/ headers are #include'd by include/exasim/*.hpp via the
# angle-form `<backend/...>` (B1). Install the .h / .hpp / .cpp files
# (.cpp files here are template-impl bodies #include'd by .h headers,
# not standalone TUs). Skip backend/Model/* — that's the legacy ABI
# which downstream consumers don't need on the templated path. The
# only exception is libpdemodel.{hpp,cpp}, which downstream consumers
# DO need: the AbiAdapter else-branch in EXASIM_DRIVER_CALL pulls
# *Driver shims that call the per-model kernels (KokkosFlux, etc.),
# and libpdemodel.cpp provides those kernel TUs.
install(DIRECTORY ${EXASIM_DIR}/backend/
    DESTINATION include/backend
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "*.cpp"
        REGEX "/Model/" EXCLUDE)
install(FILES
        ${EXASIM_DIR}/backend/Model/libpdemodel.hpp
        ${EXASIM_DIR}/backend/Model/libpdemodel.cpp
        ${EXASIM_DIR}/backend/Model/KokkosDrivers.cpp
        ${EXASIM_DIR}/backend/Model/ModelDrivers.cpp
    DESTINATION include/backend/Model)

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

# Build + install text2code as part of the main install target.
# text2code is its own CMake sub-project (text2code/text2code/) that
# emits its binary into ${EXASIM_DIR}/build/text2code. We use
# ExternalProject_Add so the main build/install drive it transparently.
include(ExternalProject)
set(_t2c_src     "${EXASIM_DIR}/text2code/text2code")
set(_t2c_build   "${CMAKE_CURRENT_BINARY_DIR}/text2code-build")
set(_t2c_binary  "${EXASIM_DIR}/build/text2code")
ExternalProject_Add(text2code_ext
    SOURCE_DIR        ${_t2c_src}
    BINARY_DIR        ${_t2c_build}
    CMAKE_ARGS        -DEXASIM_ROOT=${EXASIM_DIR}
                      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    BUILD_BYPRODUCTS  ${_t2c_binary}
    INSTALL_COMMAND   ""
    BUILD_ALWAYS      OFF)
install(PROGRAMS ${_t2c_binary} DESTINATION bin
        OPTIONAL)  # OPTIONAL so the install doesn't fail if text2code
                   # was disabled (the rest of Exasim still works).

# ======================================================================
#  HOT.3 — user-facing examples driven by CSolution<MyModel>
# ======================================================================
#
# Each example is a small directory under apps/library_example/<name>/
# containing my_model.hpp + main.cpp. Built only on the CPU/MPI legs
# (the GPU build matrix is left for a follow-up).

option(EXASIM_BUILD_LIBRARY_EXAMPLES
       "Build user-facing examples that consume Exasim as a header-only template library"
       OFF)

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
        Kokkos::kokkos ${lapackblas_libraries})
    exasim_link_pdemodel(poisson2d_template)
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
        Kokkos::kokkos ${lapackblas_libraries})
    exasim_link_pdemodel(solve_square)
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
        Kokkos::kokkos ${lapackblas_libraries})
    exasim_link_pdemodel(solve_square_facade)
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
        Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
    exasim_link_pdemodel(solve_square_facade_mpi)
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
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    exasim_link_pdemodel(solve_square_facade_gpu)
    else()
      target_compile_definitions(solve_square_facade_gpu PRIVATE _HIP)
      target_link_libraries(solve_square_facade_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries})
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
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    exasim_link_pdemodel(solve_square_facade_mpi_gpu)
    else()
      target_compile_definitions(solve_square_facade_mpi_gpu PRIVATE _HIP)
      target_link_libraries(solve_square_facade_mpi_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
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
          ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
    exasim_link_pdemodel(poisson2d_template_gpu)
    else()
      target_compile_definitions(poisson2d_template_gpu PRIVATE _HIP)
      target_link_libraries(poisson2d_template_gpu PRIVATE
          Kokkos::kokkos ${lapackblas_libraries})
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
        Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
    exasim_link_pdemodel(poisson2d_template_mpi)
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

    # HOT.7.14 — every codegen + facade target for `name` now picks up
    # its `my_model.hpp` and dlopens its `libpdemodel*.{so,dylib}`
    # from this per-example dir, not the shared backend/Model/. Drops
    # the test-matrix RESOURCE_LOCK contention so gates can run -j8
    # within a build dir AND in parallel across build dirs.
    set(_codegen_dir "${EXASIM_DIR}/apps/library_example/${name}_codegen")

    # Tell make that main.cpp.o (and main_facade.cpp.o for this
    # example) depend on the per-example `my_model.hpp` produced by
    # `regenerate.sh` during ctest. The compiler-emitted `.d` only
    # records `backend/Model/my_model.hpp` (the include-path
    # fallback that exists at initial cmake-build time, before any
    # regen has run), so without an explicit OBJECT_DEPENDS the
    # subsequent regen flips a file that make doesn't watch and
    # `cmake --build --target` is a no-op — the binary keeps
    # text2code's default ModelD physics regardless of which
    # example it was supposed to be. RESOURCE_LOCK above
    # serializes regen vs build; this dep makes the build actually
    # rebuild.
    #
    # `_facade_main` is shared across all examples, so each call
    # appends its own per-example dep; the resulting list triggers
    # over-recompilation of facade `.o`s when any example's
    # `my_model.hpp` changes, but each recompile picks up its own
    # per-example file via the per-target -I order, so the result
    # is still correct.
    set(_per_example_my_model "${_codegen_dir}/my_model.hpp")
    set_property(SOURCE "${EXASIM_DIR}/apps/library_example/${name}_codegen/main.cpp"
                 APPEND PROPERTY OBJECT_DEPENDS "${_per_example_my_model}")
    set_property(SOURCE "${EXASIM_DIR}/apps/library_example/main_facade.cpp"
                 APPEND PROPERTY OBJECT_DEPENDS "${_per_example_my_model}")

    # Per-example libpdemodel.cpp (written by regenerate.sh). Compiling
    # *from this path* makes `#include "cpuInitu.cpp"` etc resolve to
    # the per-example codegen output rather than backend/Model/'s
    # zero-fill stubs (see exasim_link_pdemodel comment above).
    set(_codegen_pdemodel "${_codegen_dir}/libpdemodel.cpp")

    # CPU variant (single-rank, no GPU).
    if(EXASIM_NOMPI AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen "${_main}")
      target_compile_definitions(${name}_codegen PRIVATE _TEXT2CODE)
      target_compile_options(${name}_codegen PRIVATE -std=c++17 -w -O2)
      target_include_directories(${name}_codegen PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_codegen PRIVATE "${_codegen_dir}")
      target_link_libraries(${name}_codegen PRIVATE
          Kokkos::kokkos ${lapackblas_libraries})
      exasim_link_pdemodel(${name}_codegen "${_codegen_pdemodel}")
      link_metis(${name}_codegen)
      set_target_properties(${name}_codegen PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()

    # MPI variant (multi-rank CPU).
    if(EXASIM_MPI AND MPI_FOUND AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_codegen_mpi "${_main}")
      target_compile_definitions(${name}_codegen_mpi PRIVATE _MPI _TEXT2CODE)
      target_compile_options(${name}_codegen_mpi PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
      target_include_directories(${name}_codegen_mpi PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_codegen_mpi PRIVATE "${_codegen_dir}")
      target_link_libraries(${name}_codegen_mpi PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
      exasim_link_pdemodel(${name}_codegen_mpi "${_codegen_pdemodel}")
      link_metis(${name}_codegen_mpi)
      set_target_properties(${name}_codegen_mpi PROPERTIES BUILD_RPATH "${_codegen_dir}")
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
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_codegen_gpu PRIVATE "${_codegen_dir}")
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_codegen_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_codegen_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
      exasim_link_pdemodel(${name}_codegen_gpu "${_codegen_pdemodel}")
      else()
        target_compile_definitions(${name}_codegen_gpu PRIVATE _HIP)
        target_link_libraries(${name}_codegen_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries})
      exasim_link_pdemodel(${name}_codegen_gpu "${_codegen_pdemodel}")
      endif()
      link_metis(${name}_codegen_gpu)
      set_target_properties(${name}_codegen_gpu PROPERTIES BUILD_RPATH "${_codegen_dir}")
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
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_codegen_mpi_gpu PRIVATE "${_codegen_dir}")
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_codegen_mpi_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_codegen_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
      exasim_link_pdemodel(${name}_codegen_mpi_gpu "${_codegen_pdemodel}")
      else()
        target_compile_definitions(${name}_codegen_mpi_gpu PRIVATE _HIP)
        target_link_libraries(${name}_codegen_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
      exasim_link_pdemodel(${name}_codegen_mpi_gpu "${_codegen_pdemodel}")
      endif()
      link_metis(${name}_codegen_mpi_gpu)
      set_target_properties(${name}_codegen_mpi_gpu PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()

    # ---------------- HOT.7.9: facade-driven variants -------------------
    # Same model + pdeapp.txt + grid.bin as the codegen targets above,
    # but the `main` is `apps/library_example/main_facade.cpp` which
    # routes through `ExasimSolver<GeneratedModel>`. Each example's
    # `my_model.hpp` is included via `target_include_directories` of
    # its codegen dir (so the same file is consumed by both mains).
    set(_facade_main "${EXASIM_DIR}/apps/library_example/main_facade.cpp")

    if(EXASIM_NOMPI AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_facade "${_facade_main}")
      target_compile_definitions(${name}_facade PRIVATE _TEXT2CODE)
      target_compile_options(${name}_facade PRIVATE -std=c++17 -w -O2)
      target_include_directories(${name}_facade PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_facade PRIVATE "${_codegen_dir}")
      target_link_libraries(${name}_facade PRIVATE
          Kokkos::kokkos ${lapackblas_libraries})
      exasim_link_pdemodel(${name}_facade "${_codegen_pdemodel}")
      link_metis(${name}_facade)
      set_target_properties(${name}_facade PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()

    if(EXASIM_MPI AND MPI_FOUND AND NOT (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_facade_mpi "${_facade_main}")
      target_compile_definitions(${name}_facade_mpi PRIVATE _MPI _TEXT2CODE)
      target_compile_options(${name}_facade_mpi PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
      target_include_directories(${name}_facade_mpi PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_facade_mpi PRIVATE "${_codegen_dir}")
      target_link_libraries(${name}_facade_mpi PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
      exasim_link_pdemodel(${name}_facade_mpi "${_codegen_pdemodel}")
      link_metis(${name}_facade_mpi)
      set_target_properties(${name}_facade_mpi PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()

    if(EXASIM_NOMPI AND (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_facade_gpu "${_facade_main}")
      target_compile_definitions(${name}_facade_gpu PRIVATE _TEXT2CODE)
      target_compile_options(${name}_facade_gpu PRIVATE
          -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_include_directories(${name}_facade_gpu PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_facade_gpu PRIVATE "${_codegen_dir}")
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_facade_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_facade_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
      exasim_link_pdemodel(${name}_facade_gpu "${_codegen_pdemodel}")
      else()
        target_compile_definitions(${name}_facade_gpu PRIVATE _HIP)
        target_link_libraries(${name}_facade_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries})
      exasim_link_pdemodel(${name}_facade_gpu "${_codegen_pdemodel}")
      endif()
      link_metis(${name}_facade_gpu)
      set_target_properties(${name}_facade_gpu PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()

    if(EXASIM_MPI AND MPI_FOUND AND (EXASIM_CUDA OR EXASIM_HIP))
      add_executable(${name}_facade_mpi_gpu "${_facade_main}")
      target_compile_definitions(${name}_facade_mpi_gpu PRIVATE _MPI _TEXT2CODE)
      target_compile_options(${name}_facade_mpi_gpu PRIVATE
          -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_include_directories(${name}_facade_mpi_gpu PRIVATE
          ${EXASIM_DIR}/include "${_codegen_dir}")
      target_link_directories(${name}_facade_mpi_gpu PRIVATE "${_codegen_dir}")
      if(EXASIM_CUDA)
        target_compile_definitions(${name}_facade_mpi_gpu PRIVATE _CUDA)
        target_link_libraries(${name}_facade_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES}
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
      exasim_link_pdemodel(${name}_facade_mpi_gpu "${_codegen_pdemodel}")
      else()
        target_compile_definitions(${name}_facade_mpi_gpu PRIVATE _HIP)
        target_link_libraries(${name}_facade_mpi_gpu PRIVATE
            Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
      exasim_link_pdemodel(${name}_facade_mpi_gpu "${_codegen_pdemodel}")
      endif()
      link_metis(${name}_facade_mpi_gpu)
      set_target_properties(${name}_facade_mpi_gpu PROPERTIES BUILD_RPATH "${_codegen_dir}")
    endif()
  endfunction()

  # Tutorial subprojects. Building these always (when the
  # library examples are on) gives a self-contained walkthrough.
  add_subdirectory(${EXASIM_DIR}/tutorial ${CMAKE_CURRENT_BINARY_DIR}/tutorial)

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

  # ======================================================================
  # CTest integration (HOT.7.9)
  #
  # `ctest` invokes apps/library_example/validate_codegen.sh per
  # backend × example × family (codegen vs facade). Each test:
  #   1. Runs `regenerate.sh <example>` to populate datain / my_model.hpp.
  #   2. Runs the binary with the app's pdeapp.txt.
  #   3. Diffs dataout/outqoi.txt against baseline/<example>_serial/.
  #
  # CTest fixtures wire the dependencies so a fresh `ctest` invocation
  # builds text2code, builds the binary for the test, and runs the
  # gate — all in the right order:
  #
  #   _build_text2code              (FIXTURES_SETUP    text2code)
  #     ↓
  #   _build_<example>_<binary>     (FIXTURES_SETUP    build_<example>_<binary>)
  #     ↓
  #   <family>:<example>:<variant>  (FIXTURES_REQUIRED both above)
  #     └── runs validate_codegen.sh which itself does
  #         regenerate.sh + binary + diff.
  #
  # HOT.7.14: gates can run -j8 in parallel — text2code writes per-
  # example libpdemodel*.{so,dylib} into apps/library_example/
  # <name>_codegen/, and each binary's BUILD_RPATH points there
  # (set by add_codegen_example()). No more shared backend/Model/
  # libpdemodel.so → no RESOURCE_LOCK contention.
  #
  # `ctest` resolves the DAG; users just say `ctest -j` and the
  # binaries get built automatically.
  # ======================================================================
  enable_testing()

  set(_codegen_examples poisson2d poisson3d periodic naca0012steady
                        lshape isoq3d cone orion nsmach8 isoq sharpb2
                        naca0012unsteady)

  # Variant + binary suffix derived from the build's configure flags.
  set(_validate_extra_args "")
  if (EXASIM_NOMPI AND NOT (EXASIM_CUDA OR EXASIM_HIP))
    set(_variant "cpu")
    set(_target_suffix "")
  elseif (EXASIM_NOMPI AND (EXASIM_CUDA OR EXASIM_HIP))
    set(_variant "gpu")
    set(_target_suffix "_gpu")
    set(_validate_extra_args "--variant" "gpu")
  elseif (EXASIM_MPI AND MPI_FOUND AND NOT (EXASIM_CUDA OR EXASIM_HIP))
    set(_variant "mpi")
    set(_target_suffix "_mpi")
    set(_validate_extra_args "--variant" "mpi" "--np" "2")
  elseif (EXASIM_MPI AND MPI_FOUND AND (EXASIM_CUDA OR EXASIM_HIP))
    set(_variant "mpi_gpu")
    set(_target_suffix "_mpi_gpu")
    set(_validate_extra_args "--variant" "mpi_gpu" "--np" "2")
  endif()

  # ----- RESOURCE_LOCK: backend_model_files ----------------------------
  #
  # text2code's POST_BUILD and regenerate.sh both rewrite shared files
  # under backend/Model/ — chiefly Symbolic{Functions,ScalarsVectors}
  # .{cpp,hpp}, libpdemodel.cpp, and my_model.hpp. backend/Model/ is on
  # the include path of every codegen/facade target's libpdemodel.cpp
  # compile (via the per-target -I lookup that falls back to
  # backend/Model when the per-example dir is missing a header).
  #
  # Without serialization, parallel ctest gates produce stale .o files:
  # fixture _regen_X writes Symbolic*.cpp (X physics) while another
  # target's libpdemodel.cpp.o is mid-compile, baking torn content into
  # the binary. The HOT.7.14 per-target dataout dirs fixed dataout
  # contention but didn't cover the source-file race; this lock does.
  #
  # Tests that hold the lock: _build_text2code, every _regen_*, every
  # _build_*_codegen / _build_*_facade, and tutorial:run-all:* (it
  # invokes text2code + cmake --build itself). The validation gates
  # (codegen:X:V / facade:X:V) only execute pre-built binaries and
  # write per-example dataout — they don't need the lock.
  set(_backend_lock "backend_model_files")

  # Setup fixture: build text2code. text2code lives in a separate
  # cmake project (text2code/text2code/) but its binary lands at
  # ${EXASIM_DIR}/build/text2code regardless of where its own build
  # dir is. The fixture runs cmake --build on its build dir if
  # present; otherwise the test errors with a helpful message at
  # its first invocation rather than silently using a stale binary.
  add_test(NAME _build_text2code
           COMMAND ${CMAKE_COMMAND} --build "${EXASIM_DIR}/text2code/text2code/build"
                                    --parallel)
  set_tests_properties(_build_text2code PROPERTIES
      FIXTURES_SETUP text2code
      RESOURCE_LOCK  ${_backend_lock}
      TIMEOUT        300)

  foreach (_ex IN LISTS _codegen_examples)
    # Per-example regenerate fixture — runs ONCE per example to produce
    # libpdemodel*.{so,dylib} + my_model.hpp in the example dir, which
    # the codegen + facade binaries link against. The per-target dirs
    # cover dataout, but text2code still writes shared backend/Model/
    # source files, so the lock above serializes regens against builds.
    add_test(NAME _regen_${_ex}
             COMMAND bash "${EXASIM_DIR}/apps/library_example/regenerate.sh" "${_ex}")
    set_tests_properties(_regen_${_ex} PROPERTIES
        FIXTURES_SETUP    regen_${_ex}
        FIXTURES_REQUIRED text2code
        RESOURCE_LOCK     ${_backend_lock}
        TIMEOUT           300)

    foreach (_stem codegen facade)
      set(_target  "${_ex}_${_stem}${_target_suffix}")
      set(_fixture "build_${_target}")

      # Per-target build. The link step needs the example's
      # libpdemodelserial.dylib (or .so / libpdemodelcuda.so); ensure
      # the regen fixture has populated it before cmake --build.
      add_test(NAME _build_${_target}
               COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}"
                                        --target ${_target} --parallel)
      set_tests_properties(_build_${_target} PROPERTIES
          FIXTURES_SETUP    ${_fixture}
          FIXTURES_REQUIRED regen_${_ex}
          RESOURCE_LOCK     ${_backend_lock}
          TIMEOUT           600)

      # Gate: regenerate already ran, binary already linked, just
      # run the binary and diff. validate_codegen.sh skips its
      # internal regen via --no-regen.
      set(_test_args "${_validate_extra_args}" "--no-regen" "--build" "${CMAKE_BINARY_DIR}")
      if (_stem STREQUAL "facade")
        list(PREPEND _test_args "--facade")
      endif()

      add_test(NAME "${_stem}:${_ex}:${_variant}"
               COMMAND bash "${EXASIM_DIR}/apps/library_example/validate_codegen.sh"
                       ${_test_args} "${_ex}"
               WORKING_DIRECTORY "${EXASIM_DIR}")
      # 1800s gate timeout — Mac runs cone (104k elems) and
      # sharpb2 (15k) at 134/144s, but lanka CPU under -j8 contention
      # is ~10x slower per-thread on these heavy NS meshes. 1200s
      # was still tripping. GPU / MPI variants are much faster and
      # finish well under this.
      set_tests_properties("${_stem}:${_ex}:${_variant}" PROPERTIES
          FIXTURES_REQUIRED       "regen_${_ex};${_fixture}"
          PASS_REGULAR_EXPRESSION "\\[ OK \\]"
          TIMEOUT                 1800)
    endforeach()
  endforeach()

  # End-to-end tutorial smoke for the current backend variant. The
  # six sections (codegen+AbiAdapter / codegen+CLI / codegen+facade /
  # hand-written+CLI / hand-written+facade / hand-written+facade-MPI)
  # are exercised against poisson2d.
  set(_tut_args "--variant" "${_variant}" "--build" "${CMAKE_BINARY_DIR}")
  if (EXASIM_MPI AND MPI_FOUND)
    list(APPEND _tut_args "--np" "2")
  endif()
  add_test(NAME "tutorial:run-all:${_variant}"
           COMMAND bash "${EXASIM_DIR}/tutorial/run-all.sh" ${_tut_args}
           WORKING_DIRECTORY "${EXASIM_DIR}")
  # tutorial/run-all.sh invokes text2code (writing shared backend/Model
  # files) and `cmake --build` on per-section targets, so it needs the
  # same lock as _regen_* / _build_*_*.
  set_tests_properties("tutorial:run-all:${_variant}" PROPERTIES
      FIXTURES_REQUIRED "regen_poisson2d;build_poisson2d_codegen${_target_suffix};build_poisson2d_facade${_target_suffix}"
      RESOURCE_LOCK     ${_backend_lock}
      TIMEOUT           600)
endif()

# ======================================================================
#  Reusable-library tests (PR #73 review B-Test1, B-Test2)
# ======================================================================
#
# Off by default so the standard build is unaffected. Some of these
# targets are expected to FAIL until the matching review blocker is
# fixed — the failure is the regression check.
option(EXASIM_BUILD_TESTS
       "Build reusable-library smoke tests under tests/ (multi-TU consumer, install consumer)"
       OFF)
if(EXASIM_BUILD_TESTS)
  add_subdirectory(${EXASIM_DIR}/tests ${CMAKE_CURRENT_BINARY_DIR}/tests)
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
message(STATUS "EXASIM_ASAN               = ${EXASIM_ASAN}")

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