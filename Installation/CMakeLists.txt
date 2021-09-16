################################################################################                   
#                    Molecular Dynamics Potentials (MDP)
#                           CESMIX-MIT Project  
# 
# Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
################################################################################
# mkdir exec
# cd exec
# cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_VERSION01=ON ../Installation 
# cmake -D EXASIM_APP=ON -D EXASIM_CORES=ON -D EXASIM_EXECUTABLES=ON -D EXASIM_CUDA=ON -D EXASIM_VERSION01=ON ../Installation  
# cmake --build .
################################################################################

cmake_minimum_required(VERSION 3.12)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 70)
endif()

project(exasim LANGUAGES CXX)

get_filename_component(EXASIM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

#set(EXASIM_SOURCE_DIR     ${EXASIM_DIR}/Version0.1/Kernel)
if(EXASIM_VERSION02)
    set(EXASIM_MAIN_DIR       ${EXASIM_DIR}/Version0.2/Kernel/Main)
elseif(EXASIM_VERSION03)
    set(EXASIM_MAIN_DIR       ${EXASIM_DIR}/Version0.3/Kernel/Main)
else()    
    set(EXASIM_MAIN_DIR       ${EXASIM_DIR}/Version0.1/Kernel/Main)
endif()

set(EXASIM_CORE_DIR       ${EXASIM_DIR}/Core)

if(EXASIM_APPDIR)
    set(EXASIM_APP_DIR        ${EXASIM_APPDIR})
else()
    set(EXASIM_APP_DIR        ${EXASIM_DIR}/Applications/App)
endif()
#message(${EXASIM_APP_DIR})

file(GLOB MAIN_SOURCES ${EXASIM_MAIN_DIR}/main.cpp)
file(GLOB COMMONCORE_SOURCES ${EXASIM_CORE_DIR}/commonCore.cpp)
file(GLOB CPUCORE_SOURCES ${EXASIM_CORE_DIR}/opuCore.cpp)
file(GLOB GPUCORE_SOURCES ${EXASIM_CORE_DIR}/gpuCore.cu)
file(GLOB CPUAPP_SOURCES ${EXASIM_APP_DIR}/opuApp.cpp)
file(GLOB GPUAPP_SOURCES ${EXASIM_APP_DIR}/gpuApp.cu)    

if(EXASIM_APP)
    add_library(cpuApp SHARED ${CPUAPP_SOURCES})
    target_compile_options(cpuApp PRIVATE -ffast-math -O3 -DNDEBUG) 

    if(EXASIM_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)        

        cuda_add_library(gpuApp SHARED ${GPUAPP_SOURCES} OPTIONS -arch sm_70 -std c++11 -O3 -w -D_FORCE_INLINES -DNDEBUG --use_fast_math)    
        target_compile_options(gpuApp PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++11 --use_fast_math -arch=sm_70 -D_FORCE_INLINES -O3 -w -DNDEBUG>) 
    endif()
endif()

if(EXASIM_CORES)
    add_library(cpuCore STATIC ${CPUCORE_SOURCES})
    target_compile_options(cpuCore PRIVATE -ffast-math -O3 -DNDEBUG) 

    add_library(commonCore STATIC ${COMMONCORE_SOURCES})
    target_compile_options(commonCore PRIVATE -ffast-math -O3 -DNDEBUG) 

    if(EXASIM_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)        

        cuda_add_library(gpuCore STATIC ${GPUCORE_SOURCES} OPTIONS -arch sm_70 -std c++11 -O3 -w -D_FORCE_INLINES -DNDEBUG --use_fast_math)
        target_compile_options(gpuCore PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++11 --use_fast_math -arch=sm_70 -D_FORCE_INLINES -O3 -w -DNDEBUG>) 
    endif()
endif()

if(EXASIM_EXECUTABLES)

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND AND BLAS_FOUND)
       set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()

    find_package(MPI)

    if(EXASIM_CUDA)
        enable_language(CUDA)
        find_package(CUDA REQUIRED)

        add_executable(gpuEXASIM ${MAIN_SOURCES})
        target_compile_definitions(gpuEXASIM PRIVATE -D_CUDA)    
        target_compile_options(gpuEXASIM PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)   
        target_link_libraries(gpuEXASIM PRIVATE commonCore cpuCore cpuApp gpuCore gpuApp ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})          

        if(MPI_FOUND)
            include_directories(SYSTEM ${MPI_INCLUDE_PATH})
            add_executable(gpumpiEXASIM ${MAIN_SOURCES})
            target_compile_definitions(gpumpiEXASIM PRIVATE -D_MPI -D_CUDA)    
            target_compile_options(gpumpiEXASIM PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG -fno-unroll-loops)             
            target_link_libraries(gpumpiEXASIM PRIVATE commonCore cpuCore cpuApp gpuCore gpuApp ${lapackblas_libraries} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${MPI_LIBRARIES})            
        endif()
    else()
        add_executable(cpuEXASIM ${MAIN_SOURCES})
        target_link_libraries(cpuEXASIM PRIVATE commonCore cpuCore cpuApp ${lapackblas_libraries})
        target_compile_options(cpuEXASIM PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG) 

        if(MPI_FOUND)
            include_directories(${MPI_INCLUDE_PATH})
            add_executable(cpumpiEXASIM ${MAIN_SOURCES})
            target_compile_definitions(cpumpiEXASIM PRIVATE -D_MPI)    
            target_compile_options(cpumpiEXASIM PRIVATE -std=c++11 -ffast-math -O3 -DNDEBUG)   
            target_link_libraries(cpumpiEXASIM PRIVATE commonCore cpuCore cpuApp ${lapackblas_libraries} ${MPI_LIBRARIES})                      
        endif()    
    endif()
endif()

