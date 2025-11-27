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

    if(WITH_TEXT2CODE)
      add_executable(gput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _CUDA _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE
        Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES}
        ${T2C_CUDA_LIB})
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

      if(WITH_TEXT2CODE)
        add_executable(gpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _CUDA _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES}
          ${T2C_CUDA_LIB})
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

    if(WITH_TEXT2CODE)
      add_executable(gput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(gput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gput2cEXASIM PRIVATE _HIP _TEXT2CODE)
      target_compile_options(gput2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
      target_link_libraries(gput2cEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${T2C_HIP_LIB})
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

      if(WITH_TEXT2CODE)
        add_executable(gpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpit2cEXASIM PRIVATE _MPI _HIP _TEXT2CODE)
        target_compile_options(gpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG -fno-unroll-loops)
        target_link_libraries(gpumpit2cEXASIM PRIVATE
          Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES}
          ${T2C_HIP_LIB})
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
    #target_compile_options(cpuEXASIM PRIVATE -std=c++17 -fsanitize=address,undefined -fno-omit-frame-pointer -g)
    #target_link_options(cpuEXASIM PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)

    if(WITH_TEXT2CODE)
      add_executable(cput2cEXASIM ${MAIN_SOURCES})
      target_link_directories(cput2cEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(cput2cEXASIM PRIVATE _TEXT2CODE)
      target_link_libraries(cput2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${T2C_CPU_LIB})
      target_compile_options(cput2cEXASIM PRIVATE -std=c++17 -w -O2)   
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

      if(WITH_TEXT2CODE)
        add_executable(cpumpit2cEXASIM ${MAIN_SOURCES})
        target_link_directories(cpumpit2cEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(cpumpit2cEXASIM PRIVATE _MPI _TEXT2CODE)
        target_compile_options(cpumpit2cEXASIM PRIVATE -std=c++17 -w -ffast-math -O3 -DNDEBUG)
        target_link_libraries(cpumpit2cEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} ${T2C_CPU_LIB})
        set_target_properties(cpumpit2cEXASIM PROPERTIES BUILD_RPATH "${Model_LIB_DIR}")
      endif()
    else()
      message("MPI library is not found")
    endif()
  endif()
endif()
