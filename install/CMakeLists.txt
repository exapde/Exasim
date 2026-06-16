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

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

get_filename_component(EXASIM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)
set(EXASIM_INCLUDE_DIR ${EXASIM_DIR}/include)
# All build outputs live in the binary dir; the source tree is never written.
set(EXASIM_LIB_DIR ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${EXASIM_LIB_DIR})
include_directories(${EXASIM_INCLUDE_DIR})

# ---- Runtime data location (replaces exasimpath for master/gauss node files) ----
# Install the precomputed node data and bake its install location into the build
# so the runtime finds it without a user-supplied exasimpath. $EXASIM_DATA_DIR
# overrides at runtime for a relocated install (see backend/Common/exasim_paths.h).
set(EXASIM_DATA_DEST "${CMAKE_INSTALL_DATAROOTDIR}/exasim")
set(EXASIM_INSTALL_DATADIR "${CMAKE_INSTALL_PREFIX}/${EXASIM_DATA_DEST}")
set(EXASIM_GENINC_DIR "${CMAKE_BINARY_DIR}/exasim_gen")
configure_file("${EXASIM_DIR}/backend/Common/exasim_paths.h.in"
               "${EXASIM_GENINC_DIR}/exasim_paths.h" @ONLY)
include_directories("${EXASIM_GENINC_DIR}")
install(FILES "${EXASIM_DIR}/backend/Preprocessing/masternodes.bin"
              "${EXASIM_DIR}/backend/Preprocessing/gaussnodes.bin"
        DESTINATION "${EXASIM_DATA_DEST}")
# Install the generated header so out-of-tree templated consumers that compile
# the FEM internals (which call exasim_data_dir()) resolve it from include/.
install(FILES "${EXASIM_GENINC_DIR}/exasim_paths.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
set(EXASIM_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/Exasim")
set(EXASIM_PACKAGE_VERSION "0.0.0")

add_library(ExasimHeaders INTERFACE)
add_library(Exasim::headers ALIAS ExasimHeaders)
set_target_properties(ExasimHeaders PROPERTIES EXPORT_NAME headers)
target_compile_features(ExasimHeaders INTERFACE cxx_std_20)
target_include_directories(ExasimHeaders INTERFACE
  $<BUILD_INTERFACE:${EXASIM_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Exasim's main program is located here  
set(EXASIM_MAIN_DIR ${EXASIM_DIR}/backend/Main)
file(GLOB MAIN_SOURCES ${EXASIM_MAIN_DIR}/main.cpp)
set(T2C_MAIN_SOURCES ${EXASIM_MAIN_DIR}/maint2c.cpp)
set(KOKKOSKERNEL_MAIN_SOURCES ${EXASIM_MAIN_DIR}/mainkernel.cpp)
set(EXASIM_SOLVER_SOURCES ${EXASIM_MAIN_DIR}/ExasimSolver.cpp)

function(add_exasim_static_library target output_name)
  add_library(${target} STATIC ${EXASIM_SOLVER_SOURCES})
  target_link_libraries(${target} PUBLIC ExasimHeaders)
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${EXASIM_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  set_target_properties(${target} PROPERTIES
    OUTPUT_NAME ${output_name}
    ARCHIVE_OUTPUT_DIRECTORY ${EXASIM_LIB_DIR})
endfunction()

function(add_exasim_executable_from_static target library_target)
  add_executable(${target} ${MAIN_SOURCES})
  target_include_directories(${target} PRIVATE ${EXASIM_INCLUDE_DIR})
  target_link_libraries(${target} PRIVATE ${library_target})
endfunction()

# Create a model-linked solver executable (text2code or built-in model variant),
# consolidating the per-backend t2c/builtin blocks. Backend-specific flags
# (OPTIONS, EXTRA_LIBS) are passed by the caller so behavior is preserved exactly.
#   add_exasim_model_executable(<name>
#     MAIN_SOURCES <...> LIB_DIR <dir> MODEL_LIB <lib>
#     DEFINES <...> OPTIONS <...> [EXTRA_LIBS <...>])
function(add_exasim_model_executable name)
  cmake_parse_arguments(E "" "LIB_DIR;MODEL_LIB" "MAIN_SOURCES;DEFINES;OPTIONS;EXTRA_LIBS" ${ARGN})
  add_executable(${name} ${E_MAIN_SOURCES} ${EXASIM_SOLVER_SOURCES})
  target_link_directories(${name} PRIVATE ${E_LIB_DIR})
  target_compile_definitions(${name} PRIVATE ${E_DEFINES})
  target_compile_options(${name} PRIVATE ${E_OPTIONS})
  target_link_libraries(${name} PRIVATE
    Kokkos::kokkos ${lapackblas_libraries} ${E_EXTRA_LIBS} ${E_MODEL_LIB})
  link_metis(${name})
  set_target_properties(${name} PROPERTIES BUILD_RPATH "${E_LIB_DIR}")
endfunction()

function(add_exasim_preprocessing_executable_from_static target library_target)
  add_executable(${target} ${KOKKOSKERNEL_MAIN_SOURCES})
  target_include_directories(${target} PRIVATE ${EXASIM_INCLUDE_DIR})
  target_link_directories(${target} PRIVATE ${Model_LIB_DIR})
  target_compile_definitions(${target} PRIVATE _KOKKOSKERNEL)
  if(KOKKOSKERNEL_HEADER)
    target_compile_definitions(${target} PRIVATE
      "EXASIM_KOKKOSKERNELS_HEADER=\"${KOKKOSKERNEL_HEADER}\"")
  endif()
  target_link_libraries(${target} PRIVATE ${library_target})
  link_metis(${target})
  set_target_properties(${target} PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
endfunction()

# Text2Code shared libraries live here
set(Model_LIB_DIR   ${EXASIM_LIB_DIR})

# BuiltinModel shared libraries live here

# ---- Option: build EXASIM variants that link text2code-generated libs ----
option(EXASIM_LIB "Build and install static libraries only; skip all EXASIM executables" OFF)

# ---- Option: build EXASIM variants that link text2code-generated libs ----
option(WITH_TEXT2CODE "Build extra EXASIM executables that link text2code-generated dynamic libraries" OFF)

# logical names for text2code-produced shared libs (lib<name>.so/.dylib/.dll)
set(T2C_CPU_LIB  t2cmodelserial)
set(T2C_CUDA_LIB t2cmodelcuda)
set(T2C_HIP_LIB  t2cmodelhip)

# Check which text2code shared libraries already exist so we can skip
# the *t2c* executables when the user has not yet run text2code on their
# PDE model to produce them.  Without this guard the linker would fail
# with "cannot find -lt2cmodelserial".
set(EXASIM_HAS_T2C_CPU_LIB  FALSE)
set(EXASIM_HAS_T2C_CUDA_LIB FALSE)
set(EXASIM_HAS_T2C_HIP_LIB  FALSE)
find_library(_T2C_CPU_LIB_PATH  ${T2C_CPU_LIB}  PATHS ${Model_LIB_DIR} NO_DEFAULT_PATH)
find_library(_T2C_CUDA_LIB_PATH ${T2C_CUDA_LIB} PATHS ${Model_LIB_DIR} NO_DEFAULT_PATH)
find_library(_T2C_HIP_LIB_PATH  ${T2C_HIP_LIB}  PATHS ${Model_LIB_DIR} NO_DEFAULT_PATH)
if(_T2C_CPU_LIB_PATH)
  set(EXASIM_HAS_T2C_CPU_LIB TRUE)
endif()
if(_T2C_CUDA_LIB_PATH)
  set(EXASIM_HAS_T2C_CUDA_LIB TRUE)
endif()
if(_T2C_HIP_LIB_PATH)
  set(EXASIM_HAS_T2C_HIP_LIB TRUE)
endif()

# logical names for builtin static libs (lib<name>.a)
set(BM_CPU_LIB  builtinmodelserial)
set(BM_CUDA_LIB builtinmodelcuda)
set(BM_HIP_LIB  builtinmodelhip)

function(exasim_set_default_kokkos_dir backend)
  if(Kokkos_DIR)
    return()
  endif()

  set(_kokkos_candidates
    "${EXASIM_DIR}/build_local/deps/kokkos/build${backend}"
    "${EXASIM_DIR}/build/deps/kokkos/build${backend}"
    "${EXASIM_DIR}/deps/kokkos/build${backend}"
    "${EXASIM_DIR}/kokkos/build${backend}")

  foreach(_cand IN LISTS _kokkos_candidates)
    if(EXISTS "${_cand}/KokkosConfig.cmake")
      set(Kokkos_DIR "${_cand}" CACHE PATH "Path to KokkosConfig.cmake" FORCE)
      return()
    endif()
  endforeach()
endfunction()

# ---- Kokkos: resolve before the BuiltIn subdirectory ----
# exasim_finalize_builtin_models() (called from the BuiltIn subdirectory) has an
# `if(NOT TARGET Kokkos::kokkos)` guard; it skips find_package when the target
# already exists. If the guard is NOT skipped, the function calls find_package itself
# — but for CUDA builds KokkosConfig.cmake calls FIND_DEPENDENCY(CUDAToolkit REQUIRED)
# which requires CUDAToolkit to be locatable, and that can fail before the parent's
# if(EXASIM_CUDA) block has had a chance to set things up. Detect Kokkos here first
# so the guard fires and the BuiltIn subdirectory never calls find_package on its own.
if(EXASIM_CUDA)
  if(NOT Kokkos_DIR)
    exasim_set_default_kokkos_dir(cuda)
  endif()
  find_package(Kokkos REQUIRED)
elseif(EXASIM_HIP)
  if(NOT Kokkos_DIR)
    if(EXASIM_TUOLUMNE)
      exasim_set_default_kokkos_dir(tuolumne)
    else()
      exasim_set_default_kokkos_dir(hip)
    endif()
  endif()
  find_package(Kokkos REQUIRED)
else()
  if(NOT Kokkos_DIR)
    find_package(Kokkos QUIET)
  endif()
  if(NOT Kokkos_FOUND)
    if(NOT Kokkos_DIR)
      exasim_set_default_kokkos_dir(serial)
    endif()
    find_package(Kokkos REQUIRED)
  endif()
endif()

# The built-in model library is always required: cpuEXASIM/gpuEXASIM link it
# directly. Defines the CMake target builtinmodel{serial,cuda,hip}.
add_subdirectory(${EXASIM_DIR}/backend/Model/BuiltIn
                 ${CMAKE_BINARY_DIR}/builtinmodel)

# ---- Option: build EXASIM variants that include Kokkos kernel providers ----
set(KOKKOSKERNEL_HEADER "" CACHE STRING "Header defining PdeModel for the Kokkos kernel provider")
option(WITH_KOKKOSKERNEL "Build extra EXASIM executables that use header-only Kokkos kernel providers" OFF)
if(KOKKOSKERNEL_HEADER)
  if(NOT EXISTS "${KOKKOSKERNEL_HEADER}")
    message(FATAL_ERROR "KOKKOSKERNEL_HEADER does not exist: ${KOKKOSKERNEL_HEADER}")
  endif()
  set(WITH_KOKKOSKERNEL ON)
endif()

set(EXASIM_BUILD_PREPROCESSING_LIBS ${WITH_KOKKOSKERNEL})
if(EXASIM_LIB)
  set(EXASIM_BUILD_PREPROCESSING_LIBS ON)
endif()

# ---- Option: link METIS / ParMETIS / GKlib -------------------------------
option(WITH_PARMETIS "Build extra EXASIM executables that link metis libraries" OFF)

# METIS: ParMETIS / METIS / GKlib roots (override if your layout differs)
set(PARMETIS_ROOT "${EXASIM_DIR}/deps/metis/ParMETIS"
    CACHE PATH "Path to ParMETIS root")
set(METIS_ROOT    "${EXASIM_DIR}/deps/metis/METIS"
    CACHE PATH "Path to METIS root")
set(GKLIB_ROOT    "${EXASIM_DIR}/deps/metis/GKlib"
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

# Resolve METIS / ParMETIS / GKlib: prefer an installed copy (e.g.
# `apt install libmetis-dev libparmetis-dev`, or any package-manager install on
# CMAKE_PREFIX_PATH), and fall back to the vendored build under Exasim/metis via
# HINTS. METIS ships no upstream CMake config, so we locate the library + header
# directly with find_library/find_path -- no third-party Find-module is adopted.
if(WITH_PARMETIS)
  find_path(EXASIM_METIS_INCLUDE       metis.h    HINTS "${METIS_INCLUDE_DIR}")
  find_library(EXASIM_METIS_LIBRARY    metis      HINTS "${METIS_LIBRARY_DIR}")
  find_path(EXASIM_PARMETIS_INCLUDE    parmetis.h HINTS "${PARMETIS_INCLUDE_DIR}")
  find_library(EXASIM_PARMETIS_LIBRARY parmetis   HINTS "${PARMETIS_LIBRARY_DIR}")
  # GKlib is a separate library in the vendored build but is bundled into libmetis
  # on some installs (e.g. apt); link/include it only when present.
  find_path(EXASIM_GKLIB_INCLUDE    GKlib.h HINTS "${GKLIB_INCLUDE_DIR}")
  find_library(EXASIM_GKLIB_LIBRARY GKlib   HINTS "${GKLIB_LIBRARY_DIR}")
  if(NOT EXASIM_METIS_LIBRARY OR NOT EXASIM_PARMETIS_LIBRARY)
    message(FATAL_ERROR
      "WITH_PARMETIS=ON but METIS/ParMETIS were not found. Install them (e.g. "
      "apt install libmetis-dev libparmetis-dev) or build the vendored copy "
      "(cd Exasim/metis && make metis).")
  endif()
  message(STATUS "METIS                     = ${EXASIM_METIS_LIBRARY}")
  message(STATUS "ParMETIS                  = ${EXASIM_PARMETIS_LIBRARY}")
  message(STATUS "GKlib                     = ${EXASIM_GKLIB_LIBRARY}")
endif()

# METIS: helper to attach METIS / ParMETIS / GKlib to a target when WITH_PARMETIS=ON
function(link_metis target)
  if (NOT WITH_PARMETIS)
    return()
  endif()

  # Preprocessor macro so source can guard METIS-specific code
  target_compile_definitions(${target} PUBLIC HAVE_PARMETIS)

  # Includes
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${EXASIM_PARMETIS_INCLUDE}>
    $<BUILD_INTERFACE:${EXASIM_METIS_INCLUDE}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  if (EXASIM_GKLIB_INCLUDE)
    target_include_directories(${target} PUBLIC
      $<BUILD_INTERFACE:${EXASIM_GKLIB_INCLUDE}>)
  endif()

  # Link order: ParMETIS, METIS, then GKlib (only if it is a separate library)
  target_link_libraries(${target} PUBLIC
      ${EXASIM_PARMETIS_LIBRARY}
      ${EXASIM_METIS_LIBRARY}
  )
  if (EXASIM_GKLIB_LIBRARY)
    target_link_libraries(${target} PUBLIC ${EXASIM_GKLIB_LIBRARY})
  endif()
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
  # Respect a user-provided installed Kokkos (-DKokkos_DIR); else vendored CUDA build.
  if(NOT Kokkos_DIR)
    exasim_set_default_kokkos_dir(cuda)
  endif()
  find_package(Kokkos REQUIRED)

  enable_language(CUDA)
  find_package(CUDA REQUIRED)

  include_directories(${CUDA_INCLUDE_DIRS})

  if(EXASIM_NOMPI)
    add_exasim_static_library(gpuEXASIM_static gpuEXASIM)
    target_link_directories(gpuEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
    target_compile_definitions(gpuEXASIM_static PUBLIC _CUDA)
    target_compile_options(gpuEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_link_libraries(gpuEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${BM_CUDA_LIB})
    link_metis(gpuEXASIM_static)     # METIS
    if(NOT EXASIM_LIB)
      add_exasim_executable_from_static(gpuEXASIM gpuEXASIM_static)
    endif()

    if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_CUDA_LIB)
      add_executable(gput2cEXASIM ${T2C_MAIN_SOURCES} ${EXASIM_SOLVER_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _CUDA _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
        ${T2C_CUDA_LIB})
      link_metis(gput2cEXASIM)   # METIS
      set_target_properties(gput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    if(EXASIM_BUILD_PREPROCESSING_LIBS)
      add_exasim_static_library(gpupreEXASIM_static gpupreEXASIM)
      target_link_directories(gpupreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(gpupreEXASIM_static PUBLIC _KOKKOSKERNEL _CUDA)
      target_compile_options(gpupreEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpupreEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
      link_metis(gpupreEXASIM_static)
      if(NOT EXASIM_LIB)
        add_exasim_preprocessing_executable_from_static(gpupreEXASIM gpupreEXASIM_static)
        target_compile_definitions(gpupreEXASIM PRIVATE _CUDA)
        target_compile_options(gpupreEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      endif()
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_exasim_static_library(gpumpiEXASIM_static gpumpiEXASIM)
      target_link_directories(gpumpiEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(gpumpiEXASIM_static PUBLIC _MPI _CUDA)
      target_compile_options(gpumpiEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpumpiEXASIM_static PUBLIC
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES} ${BM_CUDA_LIB})
      link_metis(gpumpiEXASIM_static)   # METIS
      if(NOT EXASIM_LIB)
        add_exasim_executable_from_static(gpumpiEXASIM gpumpiEXASIM_static)
      endif()

      if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_CUDA_LIB)
        add_executable(gpumpit2cEXASIM ${T2C_MAIN_SOURCES} ${EXASIM_SOLVER_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _CUDA _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES}
          ${T2C_CUDA_LIB})
        link_metis(gpumpit2cEXASIM) # METIS
        set_target_properties(gpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()

      if(EXASIM_BUILD_PREPROCESSING_LIBS)
        add_exasim_static_library(gpumpipreEXASIM_static gpumpipreEXASIM)
        target_link_directories(gpumpipreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
        target_compile_definitions(gpumpipreEXASIM_static PUBLIC _KOKKOSKERNEL _MPI _CUDA)
        target_compile_options(gpumpipreEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpipreEXASIM_static PUBLIC
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES})
        link_metis(gpumpipreEXASIM_static)
        if(NOT EXASIM_LIB)
          add_exasim_preprocessing_executable_from_static(gpumpipreEXASIM gpumpipreEXASIM_static)
          target_compile_definitions(gpumpipreEXASIM PRIVATE _MPI _CUDA)
          target_compile_options(gpumpipreEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        endif()
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()

# ------------------------------ HIP branch ------------------------------
elseif(EXASIM_HIP)
  # Respect a user-provided installed Kokkos (-DKokkos_DIR); else vendored HIP build.
  if(NOT Kokkos_DIR)
    if(EXASIM_TUOLUMNE)
      exasim_set_default_kokkos_dir(tuolumne)
    else()
      exasim_set_default_kokkos_dir(hip)
    endif()
  endif()
  find_package(Kokkos REQUIRED)

  enable_language(HIP)
  find_package(HIP REQUIRED)
  find_package(rocBLAS REQUIRED)
  find_package(hipBLAS REQUIRED)

  include_directories(${HIP_INCLUDE_DIRS})

  if(EXASIM_NOMPI)
    add_exasim_static_library(gpuEXASIM_static gpuEXASIM)
    target_link_directories(gpuEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
    target_compile_definitions(gpuEXASIM_static PUBLIC _HIP)
    target_compile_options(gpuEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
    target_link_libraries(gpuEXASIM_static PUBLIC Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${BM_HIP_LIB})
    link_metis(gpuEXASIM_static)     # METIS
    if(NOT EXASIM_LIB)
      add_exasim_executable_from_static(gpuEXASIM gpuEXASIM_static)
    endif()

    if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_HIP_LIB)
      add_executable(gput2cEXASIM ${T2C_MAIN_SOURCES} ${EXASIM_SOLVER_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _HIP _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${T2C_HIP_LIB})
      link_metis(gput2cEXASIM)   # METIS
      set_target_properties(gput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
    endif()

    if(EXASIM_BUILD_PREPROCESSING_LIBS)
      add_exasim_static_library(gpupreEXASIM_static gpupreEXASIM)
      target_link_directories(gpupreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(gpupreEXASIM_static PUBLIC _KOKKOSKERNEL _HIP)
      target_compile_options(gpupreEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpupreEXASIM_static PUBLIC Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries})
      link_metis(gpupreEXASIM_static)
      if(NOT EXASIM_LIB)
        add_exasim_preprocessing_executable_from_static(gpupreEXASIM gpupreEXASIM_static)
        target_compile_definitions(gpupreEXASIM PRIVATE _HIP)
        target_compile_options(gpupreEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      endif()
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_exasim_static_library(gpumpiEXASIM_static gpumpiEXASIM)
      target_link_directories(gpumpiEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(gpumpiEXASIM_static PUBLIC _MPI _HIP)
      target_compile_options(gpumpiEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gpumpiEXASIM_static PUBLIC Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES} ${BM_HIP_LIB})
      link_metis(gpumpiEXASIM_static)   # METIS
      if(NOT EXASIM_LIB)
        add_exasim_executable_from_static(gpumpiEXASIM gpumpiEXASIM_static)
      endif()

      if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_HIP_LIB)
        add_executable(gpumpit2cEXASIM ${T2C_MAIN_SOURCES} ${EXASIM_SOLVER_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _HIP _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${T2C_HIP_LIB})
        link_metis(gpumpit2cEXASIM) # METIS
        set_target_properties(gpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()

      if(EXASIM_BUILD_PREPROCESSING_LIBS)
        add_exasim_static_library(gpumpipreEXASIM_static gpumpipreEXASIM)
        target_link_directories(gpumpipreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
        target_compile_definitions(gpumpipreEXASIM_static PUBLIC _KOKKOSKERNEL _MPI _HIP)
        target_compile_options(gpumpipreEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpipreEXASIM_static PUBLIC Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES})
        link_metis(gpumpipreEXASIM_static)
        if(NOT EXASIM_LIB)
          add_exasim_preprocessing_executable_from_static(gpumpipreEXASIM gpumpipreEXASIM_static)
          target_compile_definitions(gpumpipreEXASIM PRIVATE _MPI _HIP)
          target_compile_options(gpumpipreEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        endif()
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()

# ------------------------------ CPU branch ------------------------------
else()
  # Prefer an installed Kokkos (apt/conda/spack/brew, -DKokkos_DIR, or
  # CMAKE_PREFIX_PATH); fall back to the vendored serial build.
  if(NOT Kokkos_DIR)
    find_package(Kokkos QUIET)
  endif()
  if(NOT Kokkos_FOUND)
    if(NOT Kokkos_DIR)
      exasim_set_default_kokkos_dir(serial)
    endif()
    find_package(Kokkos REQUIRED)
  endif()

  if(EXASIM_NOMPI)
    add_exasim_static_library(cpuEXASIM_static cpuEXASIM)
    target_link_directories(cpuEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
    target_link_libraries(cpuEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries} ${BM_CPU_LIB})
    target_compile_options(cpuEXASIM_static PRIVATE -std=c++20 -w -O2)
    link_metis(cpuEXASIM_static)    # METIS
    if(NOT EXASIM_LIB)
      add_exasim_executable_from_static(cpuEXASIM cpuEXASIM_static)
    endif()
    #target_compile_options(cpuEXASIM_static PRIVATE -std=c++20 -fsanitize=address,undefined -fno-omit-frame-pointer -g)
    #target_link_options(cpuEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)

    if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_CPU_LIB)
      add_exasim_model_executable(cput2cEXASIM
        MAIN_SOURCES ${T2C_MAIN_SOURCES} LIB_DIR ${Model_LIB_DIR} MODEL_LIB ${T2C_CPU_LIB}
        DEFINES _TEXT2CODE OPTIONS -std=c++20 -w -O2)
    endif()

    if(EXASIM_BUILD_PREPROCESSING_LIBS)
      add_exasim_static_library(cpupreEXASIM_static cpupreEXASIM)
      target_link_directories(cpupreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(cpupreEXASIM_static PUBLIC _KOKKOSKERNEL)
      target_link_libraries(cpupreEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries})
      target_compile_options(cpupreEXASIM_static PRIVATE -std=c++20 -w -O2)
      link_metis(cpupreEXASIM_static)
      if(NOT EXASIM_LIB)
        add_exasim_preprocessing_executable_from_static(cpupreEXASIM cpupreEXASIM_static)
        target_compile_options(cpupreEXASIM PRIVATE -std=c++20 -w -O2)
      endif()
    endif()
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if (MPI_FOUND)
      message("An MPI library is found")
      include_directories(${MPI_INCLUDE_PATH})

      add_exasim_static_library(cpumpiEXASIM_static cpumpiEXASIM)
      target_link_directories(cpumpiEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
      target_compile_definitions(cpumpiEXASIM_static PUBLIC _MPI)
      target_compile_options(cpumpiEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG)
      #target_compile_options(cpumpiEXASIM PRIVATE -std=c++20 -w -ffast-math -fsanitize=address,undefined -fno-omit-frame-pointer -g)
      #target_link_options(cpumpiEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)  
      target_link_libraries(cpumpiEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${BM_CPU_LIB})
      link_metis(cpumpiEXASIM_static)   # METIS
      if(NOT EXASIM_LIB)
        add_exasim_executable_from_static(cpumpiEXASIM cpumpiEXASIM_static)
      endif()

      if(WITH_TEXT2CODE AND EXASIM_HAS_T2C_CPU_LIB)
        add_exasim_model_executable(cpumpit2cEXASIM
          MAIN_SOURCES ${T2C_MAIN_SOURCES} LIB_DIR ${Model_LIB_DIR} MODEL_LIB ${T2C_CPU_LIB}
          DEFINES _MPI _TEXT2CODE OPTIONS -std=c++20 -w -ffast-math -O3 -DNDEBUG
          EXTRA_LIBS ${MPI_LIBRARIES})
      endif()

      if(EXASIM_BUILD_PREPROCESSING_LIBS)
        add_exasim_static_library(cpumpipreEXASIM_static cpumpipreEXASIM)
        target_link_directories(cpumpipreEXASIM_static PUBLIC $<BUILD_INTERFACE:${Model_LIB_DIR}> $<INSTALL_INTERFACE:${CMAKE_INSTALL_LIBDIR}>)
        target_compile_definitions(cpumpipreEXASIM_static PUBLIC _KOKKOSKERNEL _MPI)
        target_compile_options(cpumpipreEXASIM_static PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG)
        target_link_libraries(cpumpipreEXASIM_static PUBLIC Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES})
        link_metis(cpumpipreEXASIM_static)
        if(NOT EXASIM_LIB)
          add_exasim_preprocessing_executable_from_static(cpumpipreEXASIM cpumpipreEXASIM_static)
          target_compile_definitions(cpumpipreEXASIM PRIVATE _MPI)
          target_compile_options(cpumpipreEXASIM PRIVATE -std=c++20 -w -ffast-math -O3 -DNDEBUG)
        endif()
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()
endif()

set(EXASIM_EXPORT_TARGETS ExasimHeaders)
set(EXASIM_PACKAGE_HAS_MPI OFF)

macro(register_exasim_export_target target export_name requires_mpi)
  if(TARGET ${target})
    set_target_properties(${target} PROPERTIES EXPORT_NAME ${export_name})
    list(APPEND EXASIM_EXPORT_TARGETS ${target})
    if(${requires_mpi})
      target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
      set(EXASIM_PACKAGE_HAS_MPI ON)
    endif()
  endif()
endmacro()

register_exasim_export_target(cpuEXASIM_static cpulib FALSE)
register_exasim_export_target(cpumpiEXASIM_static cpumpilib TRUE)
register_exasim_export_target(cpupreEXASIM_static cpuprelib FALSE)
register_exasim_export_target(cpumpipreEXASIM_static cpumpiprelib TRUE)
register_exasim_export_target(gpuEXASIM_static gpulib FALSE)
register_exasim_export_target(gpumpiEXASIM_static gpumpilib TRUE)
register_exasim_export_target(gpupreEXASIM_static gpuprelib FALSE)
register_exasim_export_target(gpumpipreEXASIM_static gpumpiprelib TRUE)

# Built-in model shared libraries (built above via add_subdirectory) ship as
# part of the package, exported as Exasim::builtinmodel{serial,cuda,hip}.
register_exasim_export_target(builtinmodelserial builtinmodelserial FALSE)
register_exasim_export_target(builtinmodelcuda   builtinmodelcuda   FALSE)
register_exasim_export_target(builtinmodelhip    builtinmodelhip    FALSE)

list(LENGTH EXASIM_EXPORT_TARGETS EXASIM_NUM_EXPORT_TARGETS)
if(EXASIM_NUM_EXPORT_TARGETS EQUAL 1)
  message(FATAL_ERROR "No Exasim static solver library target was built for packaging.")
endif()

configure_package_config_file(
  ${EXASIM_DIR}/cmake/ExasimConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfig.cmake
  INSTALL_DESTINATION ${EXASIM_CMAKE_INSTALL_DIR})

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfigVersion.cmake
  VERSION ${EXASIM_PACKAGE_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(TARGETS ${EXASIM_EXPORT_TARGETS}
  EXPORT ExasimTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(DIRECTORY ${EXASIM_INCLUDE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Header-only / templated authoring path: the public <exasim/*.hpp> headers pull
# in the FEM internals under backend/ (e.g. backend/Common/common.h). Install the
# backend headers so out-of-tree consumers can compile against the templates.
# Installed under include/backend/, so the existing include/ search dir resolves
# both <exasim/...> and backend/... .
install(DIRECTORY ${EXASIM_DIR}/backend/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/backend
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

# The templated path's discretization.hpp includes one model-driver
# implementation (.cpp) selected by HAVE_TEXT2CODE / HAVE_BUILTINMODEL / else.
# Install those so out-of-tree templated consumers can compile.
install(FILES
  ${EXASIM_DIR}/backend/Model/KokkosDrivers.cpp
  ${EXASIM_DIR}/backend/Model/ModelDrivers.cpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/backend/Model)
install(FILES
  ${EXASIM_DIR}/backend/Model/BuiltIn/BuiltinModelDrivers.cpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/backend/Model/BuiltIn)

# ---- Façade headers under the exasim/ prefix -------------------------------
# The C++ façade headers (ExasimSolver.hpp, ExasimSolverSetup.hpp,
# builtinlibprovider.hpp, ...) live at include/ top-level, but consumers include
# them as <exasim/...> (matching the templated-path headers already under
# include/exasim/). Generate forwarding headers so both spellings resolve.
file(GLOB _exasim_top_headers RELATIVE ${EXASIM_INCLUDE_DIR}
  ${EXASIM_INCLUDE_DIR}/*.hpp)
set(_exasim_fwd_dir ${CMAKE_CURRENT_BINARY_DIR}/exasim_forward_headers)
file(MAKE_DIRECTORY ${_exasim_fwd_dir})
foreach(_h ${_exasim_top_headers})
  if(NOT EXISTS ${EXASIM_INCLUDE_DIR}/exasim/${_h})
    file(WRITE ${_exasim_fwd_dir}/${_h} "#pragma once\n#include \"../${_h}\"\n")
  endif()
endforeach()
install(DIRECTORY ${_exasim_fwd_dir}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/exasim)

# (Built-in model libraries are exported as real targets above, via the
#  add_subdirectory(backend/Model/BuiltIn) build and the export set.)

if(WITH_PARMETIS)
  if(EXISTS "${METIS_INCLUDE_DIR}/metis.h")
    install(FILES "${METIS_INCLUDE_DIR}/metis.h"
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  if(EXISTS "${PARMETIS_INCLUDE_DIR}/parmetis.h")
    install(FILES "${PARMETIS_INCLUDE_DIR}/parmetis.h"
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  foreach(_metis_archive
      "${PARMETIS_LIBRARY_DIR}/libparmetis.a"
      "${METIS_LIBRARY_DIR}/libmetis.a"
      "${GKLIB_LIBRARY_DIR}/libGKlib.a")
    if(EXISTS "${_metis_archive}")
      install(FILES "${_metis_archive}"
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
  endforeach()
endif()

install(EXPORT ExasimTargets
  NAMESPACE Exasim::
  DESTINATION ${EXASIM_CMAKE_INSTALL_DIR})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ExasimConfigVersion.cmake
  DESTINATION ${EXASIM_CMAKE_INSTALL_DIR})

# ---- External built-in model helpers ---------------------------------------
# ExasimExternalModel.cmake and ExternalModelProvider.cpp.in are installed
# alongside ExasimConfig.cmake so consumers can call
# exasim_add_external_builtin_model() after find_package(Exasim).
install(FILES
  ${EXASIM_DIR}/cmake/ExasimExternalModel.cmake
  ${EXASIM_DIR}/cmake/ExternalModelProvider.cpp.in
  DESTINATION ${EXASIM_CMAKE_INSTALL_DIR})

# model.hpp is installed by the backend DIRECTORY rule (*.hpp pattern), but
# model.cpp / libt2cmodel.cpp are not (the pattern excludes .cpp). Install them
# explicitly so exasim_add_external_builtin_model(PDEMODEL) can instantiate the
# builtin templates and the installed text2code can build libt2cmodel*.{so,dylib}.
install(FILES
  ${EXASIM_DIR}/backend/Model/BuiltIn/model.cpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/backend/Model/BuiltIn)
install(FILES
  ${EXASIM_DIR}/backend/Model/Text2codeGenerated/libt2cmodel.cpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/backend/Model/Text2codeGenerated)

# Install the text2code binary if the superbuild built it (it forwards the
# location via EXASIM_TEXT2CODE; out-of-tree, never under the source tree).
if(DEFINED EXASIM_TEXT2CODE AND EXISTS "${EXASIM_TEXT2CODE}")
  install(PROGRAMS "${EXASIM_TEXT2CODE}"
    DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# ---- Language frontends (Python / Julia / MATLAB) ---------------------------
# The frontends are installed as proper packages; their cmakecompile step
# renders the frontend-app templates below into a hidden per-app build dir and
# builds against this install via find_package(Exasim).
option(EXASIM_FRONTENDS "Install the Python/Julia/MATLAB frontends" ON)

if(EXASIM_FRONTENDS)
  # App project templates, rendered by each frontend's cmakecompile step.
  install(DIRECTORY ${EXASIM_DIR}/cmake/frontend-app
    DESTINATION ${EXASIM_CMAKE_INSTALL_DIR})

  # Record the absolute path of the cmake that built/installed Exasim, so the
  # frontends' cmakecompile step can invoke it directly. This matters when a
  # frontend runs from a GUI (e.g. the MATLAB app) whose PATH does not include
  # the cmake used to build Exasim (Homebrew/conda dirs are not on the launchd
  # PATH). Falls back to a bare "cmake" lookup when the file is absent or its
  # path no longer exists (relocated install).
  install(CODE "file(WRITE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${EXASIM_CMAKE_INSTALL_DIR}/cmake_command.txt\" \"${CMAKE_COMMAND}\n\")")

  # Python: install the exasim package into the prefix's site-packages. When
  # the prefix is on the interpreter's sys.path (e.g. a virtualenv prefix),
  # `import exasim` then works without any path setup; otherwise add the
  # directory to PYTHONPATH.
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_FOUND AND EXISTS ${EXASIM_DIR}/frontends/Python/exasim)
    set(EXASIM_PYTHON_SITEDIR
      lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages
      CACHE STRING "Install dir for the exasim Python package (relative to the prefix)")
    install(DIRECTORY ${EXASIM_DIR}/frontends/Python/exasim
      DESTINATION ${EXASIM_PYTHON_SITEDIR}
      # Exclude ONLY the top-level package __init__.py (regenerated from the
      # .in template below). A bare PATTERN "__init__.py" would match the
      # trailing component of every path and drop the subpackage __init__.py
      # files too, turning each subpackage into an empty namespace package.
      REGEX "/exasim/__init__\\.py$" EXCLUDE
      PATTERN "__pycache__" EXCLUDE)
    install(CODE "
      file(READ \"${EXASIM_DIR}/frontends/Python/exasim/__init__.py.in\" _exasim_python_init_template)
      string(CONFIGURE \"\${_exasim_python_init_template}\" _exasim_python_init_contents @ONLY)
      file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${EXASIM_PYTHON_SITEDIR}/exasim\")
      file(WRITE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${EXASIM_PYTHON_SITEDIR}/exasim/__init__.py\"
        \"\${_exasim_python_init_contents}\")
    ")
    install(CODE "message(STATUS \"exasim Python package: add \${CMAKE_INSTALL_PREFIX}/${EXASIM_PYTHON_SITEDIR} to PYTHONPATH unless that prefix is already on the interpreter's sys.path\")")
  endif()

  # Optionally pip-install the package into the configured interpreter at
  # install time, so `import exasim` works with no PYTHONPATH at all. Off by
  # default because it mutates the interpreter's environment; requires a
  # pip-writable interpreter (a virtualenv or similar).
  option(EXASIM_PIP_INSTALL
    "pip-install the exasim Python package into the configured interpreter at install time" OFF)
  if(EXASIM_PIP_INSTALL AND Python3_FOUND)
    install(CODE "
      message(STATUS \"pip install ${EXASIM_DIR}/frontends/Python (into ${Python3_EXECUTABLE})\")
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pip install --quiet ${EXASIM_DIR}/frontends/Python
        RESULT_VARIABLE _pipres)
      if(NOT _pipres EQUAL 0)
        message(FATAL_ERROR \"pip install failed (is the interpreter pip-writable?)\")
      endif()
      # Bake the install prefix into the pip-installed copy so it can locate
      # this C++ install without EXASIM_PREFIX.
      execute_process(
        COMMAND ${Python3_EXECUTABLE} -c \"import exasim, os; open(os.path.join(os.path.dirname(exasim.__file__), '_installed_prefix.py'), 'w').write('PREFIX = ' + repr('\${CMAKE_INSTALL_PREFIX}') + chr(10))\"
        RESULT_VARIABLE _bakeres)
      if(NOT _bakeres EQUAL 0)
        message(WARNING \"could not bake the install prefix into the pip-installed exasim package; set EXASIM_PREFIX when using it\")
      endif()
      # setuptools may drop build artifacts next to the package source; the
      # source tree stays pristine.
      file(REMOVE_RECURSE \"${EXASIM_DIR}/frontends/Python/build\"
                          \"${EXASIM_DIR}/frontends/Python/exasim.egg-info\")
    ")
  endif()

  # Julia: Exasim.jl package; use it via
  #   push!(LOAD_PATH, "<prefix>/share/exasim/julia") ; using Exasim
  # or Pkg.develop(path="<prefix>/share/exasim/julia/Exasim").
  if(EXISTS ${EXASIM_DIR}/frontends/Julia/Exasim)
    install(DIRECTORY ${EXASIM_DIR}/frontends/Julia/Exasim
      DESTINATION share/exasim/julia
      PATTERN "src/Exasim.jl" EXCLUDE)
    install(CODE "
      file(READ \"${EXASIM_DIR}/frontends/Julia/Exasim/src/Exasim.jl.in\" _exasim_julia_module_template)
      string(CONFIGURE \"\${_exasim_julia_module_template}\" _exasim_julia_module_contents @ONLY)
      file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/share/exasim/julia/Exasim/src\")
      file(WRITE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/share/exasim/julia/Exasim/src/Exasim.jl\"
        \"\${_exasim_julia_module_contents}\")
    ")

    # Optionally Pkg.develop the installed package into the default Julia
    # environment at install time, so `using Exasim` works with no LOAD_PATH
    # setup. Off by default because it mutates the user's Julia environment.
    option(EXASIM_JULIA_DEVELOP
      "Pkg.develop the installed Exasim.jl into the default Julia environment at install time" OFF)
    if(EXASIM_JULIA_DEVELOP)
      find_program(EXASIM_JULIA_EXECUTABLE julia)
      if(EXASIM_JULIA_EXECUTABLE)
        install(CODE "
          message(STATUS \"Pkg.develop \${CMAKE_INSTALL_PREFIX}/share/exasim/julia/Exasim\")
          execute_process(
            COMMAND ${EXASIM_JULIA_EXECUTABLE} -e \"using Pkg; Pkg.develop(path=raw\\\"\${CMAKE_INSTALL_PREFIX}/share/exasim/julia/Exasim\\\")\"
            RESULT_VARIABLE _jlres)
          if(NOT _jlres EQUAL 0)
            message(FATAL_ERROR \"julia Pkg.develop failed\")
          endif()
        ")
      else()
        message(WARNING "EXASIM_JULIA_DEVELOP=ON but no julia executable found")
      endif()
    endif()
  endif()

  # MATLAB: function tree plus the exasim_setup.m addpath entry point;
  #   run("<prefix>/share/exasim/matlab/exasim_setup.m") in MATLAB.
  if(EXISTS ${EXASIM_DIR}/frontends/Matlab)
    install(DIRECTORY ${EXASIM_DIR}/frontends/Matlab/
      DESTINATION share/exasim/matlab
      PATTERN "exasim_setup.m" EXCLUDE)
    install(CODE "
      file(READ \"${EXASIM_DIR}/frontends/Matlab/exasim_setup.m.in\" _exasim_setup_template)
      string(CONFIGURE \"\${_exasim_setup_template}\" _exasim_setup_contents @ONLY)
      file(MAKE_DIRECTORY \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/share/exasim/matlab\")
      file(WRITE \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/share/exasim/matlab/exasim_setup.m\"
        \"\${_exasim_setup_contents}\")
    ")
  endif()
endif()

# NOTE: the frontend entry files (exasim_setup.m / __init__.py / Exasim.jl) are
# rendered from their .in templates into the INSTALL tree above. They are
# deliberately NOT written back into the source tree: doing so bakes the
# machine-specific CMAKE_INSTALL_PREFIX into tracked files (git churn, and a
# leaked absolute path if committed). A source checkout resolves the prefix at
# runtime via config.install_prefix() / exasim_install_prefix(), which already
# probe <repo>/local and <repo>-build/install, or honor EXASIM_PREFIX.

export(EXPORT ExasimTargets
  NAMESPACE Exasim::
  FILE ${CMAKE_CURRENT_BINARY_DIR}/ExasimTargets.cmake)

# ======================================================================
#  Summary of EXASIM configuration
# ======================================================================
message(STATUS "==================== EXASIM CONFIGURATION ====================")

message(STATUS "EXASIM_DIR                = ${EXASIM_DIR}")
message(STATUS "EXASIM_INCLUDE_DIR        = ${EXASIM_INCLUDE_DIR}")
message(STATUS "EXASIM_LIB_DIR            = ${EXASIM_LIB_DIR}")
message(STATUS "EXASIM_LIB                = ${EXASIM_LIB}")
message(STATUS "EXASIM_EXPORT_TARGETS       = ${EXASIM_EXPORT_TARGETS}")
message(STATUS "EXASIM_CUDA               = ${EXASIM_CUDA}")
message(STATUS "EXASIM_HIP                = ${EXASIM_HIP}")
message(STATUS "EXASIM_MPI                = ${EXASIM_MPI}")
message(STATUS "EXASIM_NOMPI              = ${EXASIM_NOMPI}")
message(STATUS "EXASIM_TUOLUMNE           = ${EXASIM_TUOLUMNE}")

message(STATUS "WITH_TEXT2CODE            = ${WITH_TEXT2CODE}")
message(STATUS "T2C CPU library found     = ${EXASIM_HAS_T2C_CPU_LIB}")
message(STATUS "T2C CUDA library found    = ${EXASIM_HAS_T2C_CUDA_LIB}")
message(STATUS "T2C HIP library found     = ${EXASIM_HAS_T2C_HIP_LIB}")
message(STATUS "WITH_KOKKOSKERNEL         = ${WITH_KOKKOSKERNEL}")
message(STATUS "EXASIM_BUILD_PREPROCESSING_LIBS = ${EXASIM_BUILD_PREPROCESSING_LIBS}")
message(STATUS "KOKKOSKERNEL_HEADER       = ${KOKKOSKERNEL_HEADER}")
if(KOKKOSKERNEL_HEADER)
  message(STATUS "EXASIM_KOKKOSKERNELS_HEADER = ${KOKKOSKERNEL_HEADER}")
else()
  message(STATUS "EXASIM_KOKKOSKERNELS_HEADER = kokkoskernels.hpp")
endif()
message(STATUS "WITH_PARMETIS             = ${WITH_PARMETIS}")

message(STATUS "Model_LIB_DIR             = ${Model_LIB_DIR}")

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

# ---- Test suite (opt-in) ---------------------------------------------------
# Out-of-tree consumer tests (tests/consumers/) build each consumption path ONLY
# via find_package(Exasim) against an install prefix, catching export regressions.
option(EXASIM_BUILD_TESTS "Build/register the Exasim test suite (ctest)" OFF)
if(EXASIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(${EXASIM_DIR}/tests ${CMAKE_BINARY_DIR}/tests)
endif()
