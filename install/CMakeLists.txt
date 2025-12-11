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

set(EXASIM_MAIN_DIR ${EXASIM_DIR}/backend/Main)
set(Model_LIB_DIR   ${EXASIM_DIR}/backend/Model)

file(GLOB MAIN_SOURCES ${EXASIM_MAIN_DIR}/main.cpp)

# ---- Option: build EXASIM variants that link text2code-generated libs ----
option(WITH_TEXT2CODE "Build extra EXASIM executables that link text2code-generated dynamic libraries" OFF)

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

# logical names for text2code-produced shared libs (lib<name>.so/.dylib/.dll)
set(T2C_CPU_LIB  pdemodelserial)
set(T2C_CUDA_LIB pdemodelcuda)
set(T2C_HIP_LIB  pdemodelhip)

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
      #target_compile_options(cpumpiEXASIM PRIVATE -std=c++17 -ffast-math -fsanitize=address,undefined -fno-omit-frame-pointer -g)
      #target_link_options(cpumpiEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)         
      set_target_properties(cput2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
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
    else()
      message("MPI library is not found")
    endif()
  endif()
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
message(STATUS "WITH_PARMETIS             = ${WITH_PARMETIS}")

message(STATUS "PARMETIS_ROOT             = ${PARMETIS_ROOT}")
message(STATUS "METIS_ROOT                = ${METIS_ROOT}")
message(STATUS "GKLIB_ROOT                = ${GKLIB_ROOT}")

message(STATUS "PARMETIS_INCLUDE_DIR      = ${PARMETIS_INCLUDE_DIR}")
message(STATUS "PARMETIS_LIBRARY_DIR      = ${PARMETIS_LIBRARY_DIR}")

message(STATUS "METIS_INCLUDE_DIR         = ${METIS_INCLUDE_DIR}")
message(STATUS "METIS_LIBRARY_DIR         = ${METIS_LIBRARY_DIR}")

message(STATUS "GKLIB_INCLUDE_DIR         = ${GKLIB_INCLUDE_DIR}")
message(STATUS "GKLIB_LIBRARY_DIR         = ${GKLIB_LIBRARY_DIR}")

message(STATUS "T2C_CPU_LIB               = ${T2C_CPU_LIB}")
message(STATUS "T2C_CUDA_LIB              = ${T2C_CUDA_LIB}")
message(STATUS "T2C_HIP_LIB               = ${T2C_HIP_LIB}")

message(STATUS "BLAS_LIBRARIES            = ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LIBRARIES          = ${LAPACK_LIBRARIES}")

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