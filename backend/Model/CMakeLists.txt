cmake_minimum_required(VERSION 3.16)
project(PdeModelLib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Determine shared library extension
if(WIN32)
    set(SHARED_EXT ".dll")
elseif(APPLE)
    set(SHARED_EXT ".dylib")
else()
    set(SHARED_EXT ".so")
endif()

# Path to Exasim root (update this or make it configurable)
set(EXASIM_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../..)  # Adjust as needed

# Kokkos backend build paths
set(KOKKOS_SERIAL_PATH ${EXASIM_PATH}/kokkos/buildserial)
set(KOKKOS_CUDA_PATH   ${EXASIM_PATH}/kokkos/buildcuda)
set(KOKKOS_HIP_PATH    ${EXASIM_PATH}/kokkos/buildhip)

# Shared source
set(PDEMODEL_SRC ${EXASIM_PATH}/backend/Model/libpdemodel.cpp)

# Helper macro to create shared library
macro(add_backend_library backend path)
    if(EXISTS ${path})
        add_library(libpdemodel${backend} SHARED ${PDEMODEL_SRC})
        target_include_directories(libpdemodel${backend} PRIVATE ${path}/include)
        target_link_libraries(libpdemodel${backend} 
            PRIVATE 
            ${path}/lib/libkokkoscore.a 
            ${path}/lib/libkokkoscontainers.a
        )
        if(EXISTS ${path}/lib/libkokkossimd.a)
            target_link_libraries(libpdemodel${backend} PRIVATE ${path}/lib/libkokkossimd.a)
        endif()
        set_target_properties(libpdemodel${backend} PROPERTIES
            OUTPUT_NAME pdemodel${backend}
            LIBRARY_OUTPUT_DIRECTORY ${EXASIM_PATH}/backend/Model
        )
    endif()
endmacro()

add_backend_library(serial ${KOKKOS_SERIAL_PATH})
add_backend_library(cuda   ${KOKKOS_CUDA_PATH})
add_backend_library(hip    ${KOKKOS_HIP_PATH})
