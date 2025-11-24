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
file(GLOB LIB_SOURCES ${EXASIM_MAIN_DIR}/libexasim.cpp)
file(GLOB MAIN_LIB_SOURCES ${EXASIM_MAIN_DIR}/exasimmain.cpp)
file(GLOB MAIN_LIB_HEADER ${EXASIM_MAIN_DIR}/libexasim.h)

# ---- Option: build EXASIM variants that link text2code-generated libs ----
option(WITH_TEXT2CODE "Build extra EXASIM executables that link text2code-generated dynamic libraries" OFF)
option(WITH_ASSUME_GENERATOR "Assume exasim does not need a generator" OFF)



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


set(gpu_t2c_args)
if(NOT WITH_ASSUME_GENERATOR)
  list(APPEND gpu_t2c_args
    WITH_TEXT2CODE
    TEXT2CODE_LIB ${T2C_CUDA_LIB}
  )
endif()


set(hip_t2c_args)
if(NOT WITH_ASSUME_GENERATOR)
  list(APPEND hip_t2c_args
    WITH_TEXT2CODE
    TEXT2CODE_LIB ${T2C_HIP_LIB}
  )
endif()


set(cup_t2c_args)
if(NOT WITH_ASSUME_GENERATOR)
  list(APPEND cpu_t2c_args
    WITH_TEXT2CODE
    TEXT2CODE_LIB ${T2C_CPU_LIB}
  )
endif()

# Create EXASIM executable + shared library + mainlib-wrapper executable
#
# Usage examples:
#   # MPI variant
#   add_exasim_variant(cpumpi
#       EXTRA_LIBS      ${MPI_LIBRARIES}
#       WITH_TEXT2CODE
#       TEXT2CODE_LIB   ${T2C_CPU_LIB}
#   )
#
#   # Pure CPU variant (no MPI, no text2code)
#   add_exasim_variant(cpu
#       EXTRA_LIBS
#   )
#
function(add_exasim_variant variant_name)
    set(options WITH_TEXT2CODE)
    set(oneValueArgs TEXT2CODE_LIB)
    set(multiValueArgs EXTRA_LIBS)
    cmake_parse_arguments(EXASIM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Target names
    set(exe_target        "${variant_name}EXASIM")
    set(lib_target        "${variant_name}EXASIMLIB")
    set(mainlib_exe       "${variant_name}EXASIMMAINLIB")
    set(text2code_target  "${variant_name}t2cEXASIM")

    # Common compile options
    set(common_compile_opts
        -std=c++17
        -w
        -ffast-math
        -O3
        -DNDEBUG
    )

    # Decide if we add _MPI based on the variant name
    set(mpi_def)
    if(variant_name MATCHES "mpi")
        set(mpi_def _MPI)
    endif()

    # --------------------
    # Main executable
    # --------------------
    add_executable(${exe_target} ${MAIN_SOURCES})
    target_link_directories(${exe_target} PRIVATE ${Model_LIB_DIR})
    if(mpi_def)
        target_compile_definitions(${exe_target} PRIVATE ${mpi_def})
    endif()
    target_compile_options(${exe_target} PRIVATE ${common_compile_opts})
    target_link_libraries(${exe_target}
        PRIVATE
            Kokkos::kokkos
            ${lapackblas_libraries}
            ${EXASIM_EXTRA_LIBS}
    )

    # --------------------
    # Shared library
    # --------------------
    add_library(${lib_target} SHARED ${LIB_SOURCES})
    target_link_directories(${lib_target} PRIVATE ${Model_LIB_DIR})
    if(mpi_def)
        target_compile_definitions(${lib_target} PRIVATE ${mpi_def})
    endif()
    target_compile_options(${lib_target} PRIVATE ${common_compile_opts})
    target_link_libraries(${lib_target}
        PRIVATE
            Kokkos::kokkos
            ${lapackblas_libraries}
            ${EXASIM_EXTRA_LIBS}
    )

    target_include_directories(${lib_target}
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/backend/Main/
    )

    # --------------------
    # Executable that links to the shared lib
    # --------------------
    add_executable(${mainlib_exe} ${MAIN_LIB_SOURCES})
    target_link_libraries(${mainlib_exe} PRIVATE ${lib_target})

    # --------------------
    # Optional text2code executable
    # --------------------
    if(EXASIM_WITH_TEXT2CODE AND EXASIM_TEXT2CODE_LIB)
        add_executable(${text2code_target} ${MAIN_SOURCES})
        target_link_directories(${text2code_target} PRIVATE ${Model_LIB_DIR})

        # _TEXT2CODE plus optional _MPI
        if(mpi_def)
            target_compile_definitions(${text2code_target} PRIVATE ${mpi_def} _TEXT2CODE)
        else()
            target_compile_definitions(${text2code_target} PRIVATE _TEXT2CODE)
        endif()

        target_compile_options(${text2code_target} PRIVATE ${common_compile_opts})

        target_link_libraries(${text2code_target}
            PRIVATE
                Kokkos::kokkos
                ${lapackblas_libraries}
                ${EXASIM_EXTRA_LIBS}
                ${EXASIM_TEXT2CODE_LIB}
        )

        set_target_properties(${text2code_target}
            PROPERTIES
                BUILD_RPATH "${Model_LIB_DIR}"
        )
    endif()
endfunction()


# Helper for GPU (CUDA/HIP) EXASIM executables (+ optional t2c variant)
#
# Usage:
#   add_exasim_gpu_variant(gpu
#       BACKEND       CUDA|HIP
#       EXTRA_LIBS    <backend libs, MPI libs, etc.>
#       [USE_MPI]
#       [WITH_TEXT2CODE]
#       [TEXT2CODE_LIB <t2c backend lib>]
#   )
#
# Uses:
#   - MAIN_SOURCES
#   - Model_LIB_DIR
#   - lapackblas_libraries
function(add_exasim_gpu_variant variant_name)
    set(options WITH_TEXT2CODE USE_MPI)
    set(oneValueArgs BACKEND TEXT2CODE_LIB)
    set(multiValueArgs EXTRA_LIBS)
    cmake_parse_arguments(GPU "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT GPU_BACKEND)
        message(FATAL_ERROR "add_exasim_gpu_variant: BACKEND must be set to CUDA or HIP")
    endif()

    if(GPU_BACKEND STREQUAL "CUDA")
        set(base_def _CUDA)
    elseif(GPU_BACKEND STREQUAL "HIP")
        set(base_def _HIP)
    else()
        message(FATAL_ERROR "add_exasim_gpu_variant: unknown BACKEND '${GPU_BACKEND}'")
    endif()

    set(exe_target  "${variant_name}EXASIM")
    set(t2c_target  "${variant_name}t2cEXASIM")

    set(common_compile_opts
        -std=c++17
        -w
        -ffast-math
        -O3
        -DNDEBUG
        -fno-unroll-loops
    )

    # ---- main EXASIM executable ----
    add_executable(${exe_target} ${MAIN_SOURCES})
    target_link_directories(${exe_target} PRIVATE ${Model_LIB_DIR})

    set(defs ${base_def})
    if(GPU_USE_MPI)
        list(APPEND defs _MPI)
    endif()
    if(defs)
        target_compile_definitions(${exe_target} PRIVATE ${defs})
    endif()

    target_compile_options(${exe_target} PRIVATE ${common_compile_opts})
    target_link_libraries(${exe_target} PRIVATE
        Kokkos::kokkos
        ${lapackblas_libraries}
        ${GPU_EXTRA_LIBS}
    )

    # ---- optional text2code variant ----
    if(GPU_WITH_TEXT2CODE)
        if(NOT GPU_TEXT2CODE_LIB)
            message(FATAL_ERROR "add_exasim_gpu_variant: WITH_TEXT2CODE set but no TEXT2CODE_LIB provided")
        endif()

        add_executable(${t2c_target} ${MAIN_SOURCES})
        target_link_directories(${t2c_target} PRIVATE ${Model_LIB_DIR})

        set(t2c_defs ${base_def} _TEXT2CODE)
        if(GPU_USE_MPI)
            list(APPEND t2c_defs _MPI)
        endif()
        target_compile_definitions(${t2c_target} PRIVATE ${t2c_defs})
        target_compile_options(${t2c_target} PRIVATE ${common_compile_opts})

        target_link_libraries(${t2c_target} PRIVATE
            Kokkos::kokkos
            ${lapackblas_libraries}
            ${GPU_EXTRA_LIBS}
            ${GPU_TEXT2CODE_LIB}
        )

        set_target_properties(${t2c_target}
            PROPERTIES
                BUILD_RPATH "${Model_LIB_DIR}"
        )
    endif()
endfunction()



# ----------------------------- CUDA branch ------------------------------
if(EXASIM_CUDA)
  set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildcuda)
  find_package(Kokkos REQUIRED)

  enable_language(CUDA)
  find_package(CUDA REQUIRED)

  include_directories(${CUDA_INCLUDE_DIRS})

  if(EXASIM_NOMPI)
    add_exasim_gpu_variant(gpu
      BACKEND      CUDA
      EXTRA_LIBS   ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
      ${gpu_t2c_args}
    )
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if(MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_exasim_gpu_variant(gpumpi
        BACKEND      CUDA
        USE_MPI
        EXTRA_LIBS   ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES}
        ${gpu_t2c_args}
      )
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
    add_exasim_gpu_variant(gpu
      BACKEND      HIP
      EXTRA_LIBS   hipblas rocblas ${HIP_LIBRARIES}
      ${hip_t2c_args}
    )
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if(MPI_FOUND)
      message("An MPI library is found")
      include_directories(SYSTEM ${MPI_INCLUDE_PATH})

      add_exasim_gpu_variant(gpumpi
        BACKEND      HIP
        USE_MPI
        EXTRA_LIBS   hipblas rocblas ${HIP_LIBRARIES} ${MPI_LIBRARIES}
        ${hip_t2c_args}
      )
    else()
      message("MPI library is not found")
    endif()
  endif()

# ------------------------------ CPU branch ------------------------------
else()
  set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildserial)
  find_package(Kokkos REQUIRED)

  if(EXASIM_NOMPI)
    add_exasim_variant(cpu
      EXTRA_LIBS          # no extra libs beyond Kokkos + lapackblas
      ${cpu_t2c_args}
    )
  endif()

  if(EXASIM_MPI)
    find_package(MPI)
    if(MPI_FOUND)
      message("An MPI library is found")
      include_directories(${MPI_INCLUDE_PATH})

      add_exasim_variant(cpumpi
        EXTRA_LIBS    ${MPI_LIBRARIES}
        ${cpu_t2c_args}
      )
    else()
      message("MPI library is not found")
    endif()
  endif()
endif()
