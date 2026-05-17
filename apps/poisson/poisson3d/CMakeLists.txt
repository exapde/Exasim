cmake_minimum_required(VERSION 3.16)
project(exasim_app CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(EXASIM_MPI "Link the MPI-enabled Exasim preprocessing library" ON)
option(EXASIM_GPU "Link the GPU-enabled Exasim preprocessing library" OFF)

find_package(Exasim REQUIRED)

set(EXASIM_APP_LIBRARY_TARGET Exasim::cpuprelib)
if(EXASIM_MPI)
  set(EXASIM_APP_LIBRARY_TARGET Exasim::cpumpiprelib)
endif()
if(EXASIM_GPU)
  if(EXASIM_MPI)
    set(EXASIM_APP_LIBRARY_TARGET Exasim::gpumpiprelib)
  else()
    set(EXASIM_APP_LIBRARY_TARGET Exasim::gpuprelib)
  endif()
endif()

if(NOT TARGET ${EXASIM_APP_LIBRARY_TARGET})
  message(FATAL_ERROR
    "Requested Exasim target ${EXASIM_APP_LIBRARY_TARGET} is not available in the "
    "current Exasim package. Reconfigure Exasim to build that library variant.")
endif()

add_executable(exasimapp main.cpp)
target_include_directories(exasimapp PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(exasimapp PRIVATE _KOKKOSKERNEL)
if(EXASIM_MPI)
  target_compile_definitions(exasimapp PRIVATE _MPI)
endif()
target_link_libraries(exasimapp PRIVATE Exasim::headers ${EXASIM_APP_LIBRARY_TARGET})
