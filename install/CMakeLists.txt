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

set(EXASIM_MAIN_DIR       ${EXASIM_DIR}/backend/Main)
set(Model_LIB_DIR ${EXASIM_DIR}/backend/Model)

#message(${EXASIM_MAIN_DIR})
#message(${EXASIM_NOMPI})
#message(${EXASIM_MPI})

file(GLOB MAIN_SOURCES ${EXASIM_MAIN_DIR}/main.cpp)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
   message("Lapack and Blas libraries are not found")
endif()

if(EXASIM_CUDA)
    set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildcuda)
    find_package(Kokkos REQUIRED)

    enable_language(CUDA)
    find_package(CUDA REQUIRED)

    # Specify the include directories
    include_directories(${CUDA_INCLUDE_DIRS})

    if(EXASIM_NOMPI)
      add_executable(gpuEXASIM ${MAIN_SOURCES})
      target_link_directories(gpuEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gpuEXASIM PRIVATE -D_CUDA)    
      target_compile_options(gpuEXASIM PRIVATE -std=c++17 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)   
      target_link_libraries(gpuEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} pdemodelcuda)          
    endif()

    if(EXASIM_MPI)  
      find_package(MPI)
      if (MPI_FOUND)
        message("An MPI library is found")
      else()
        message("MPI library is not found")
      endif()

      if(MPI_FOUND)
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})
        add_executable(gpumpiEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpiEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpiEXASIM PRIVATE -D_MPI -D_CUDA)    
        target_compile_options(gpumpiEXASIM PRIVATE -std=c++17 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)             
        target_link_libraries(gpumpiEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES} pdemodelcuda)            
      endif()
    endif()
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

    # Specify the include directories
    include_directories(${HIP_INCLUDE_DIRS})

    if(EXASIM_NOMPI)
      add_executable(gpuEXASIM ${MAIN_SOURCES})
      target_link_directories(gpuEXASIM PRIVATE ${Model_LIB_DIR})
      target_compile_definitions(gpuEXASIM PRIVATE -D_HIP)    
      target_compile_options(gpuEXASIM PRIVATE -std=c++17 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)   
      target_link_libraries(gpuEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} pdemodelhip)          
    endif()

    if(EXASIM_MPI)  
      find_package(MPI)
      if (MPI_FOUND)
        message("An MPI library is found")
      else()
        message("MPI library is not found")
      endif()

      if(MPI_FOUND)
        include_directories(SYSTEM ${MPI_INCLUDE_PATH})
        add_executable(gpumpiEXASIM ${MAIN_SOURCES})
        target_link_directories(gpumpiEXASIM PRIVATE ${Model_LIB_DIR})
        target_compile_definitions(gpumpiEXASIM PRIVATE -D_MPI -D_HIP)    
        target_compile_options(gpumpiEXASIM PRIVATE -std=c++17 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)             
        target_link_libraries(gpumpiEXASIM PRIVATE Kokkos::kokkos hipblas rocblas ${HIP_LIBRARIES} ${lapackblas_libraries} ${MPI_LIBRARIES} pdemodelhip)            
      endif()
    endif()
else()
    set(Kokkos_DIR ${EXASIM_DIR}/kokkos/buildserial)

    #message(${EXASIM_DIR})
    #message(${Kokkos_DIR})
    #message(${Model_LIB_DIR})

    find_package(Kokkos REQUIRED)

    if(EXASIM_NOMPI)
      add_executable(cpuEXASIM ${MAIN_SOURCES})
      target_link_directories(cpuEXASIM PRIVATE ${Model_LIB_DIR})
      target_link_libraries(cpuEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} pdemodelserial)
      target_compile_options(cpuEXASIM PRIVATE -std=c++17 -Wall -Og -ggdb3 -g) 
    endif()

    if(EXASIM_MPI)  

      find_package(MPI)
      if (MPI_FOUND)
        message("An MPI library is found")
      else()
        message("MPI library is not found")
      endif()

      if(MPI_FOUND)
        include_directories(${MPI_INCLUDE_PATH})
        add_executable(cpumpiEXASIM ${MAIN_SOURCES})
        target_link_directories(cpumpiEXASIM PRIVATE ${Model_LIB_DIR}) 
        target_compile_definitions(cpumpiEXASIM PRIVATE -D_MPI)    
        target_compile_options(cpumpiEXASIM PRIVATE -std=c++17 -ffast-math -O3 -DNDEBUG)   
        target_link_libraries(cpumpiEXASIM PRIVATE Kokkos::kokkos ${lapackblas_libraries} ${MPI_LIBRARIES} pdemodelserial)                      
      endif()    
    endif()    
endif()


