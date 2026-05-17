cmake_minimum_required(VERSION 3.16)
project(poisson3d_app CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Exasim REQUIRED)

add_executable(poisson3d main.cpp)
target_include_directories(poisson3d PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(poisson3d PRIVATE _MPI _KOKKOSKERNEL)
target_link_libraries(poisson3d PRIVATE Exasim::headers Exasim::cpumpiprelib)
