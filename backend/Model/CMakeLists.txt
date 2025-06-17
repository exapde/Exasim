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

# Shared source
set(PDEMODEL_SRC ${EXASIM_PATH}/backend/Model/libpdemodel.cpp)

if(EXASIM_CUDA)
    set(Kokkos_DIR ${EXASIM_PATH}/kokkos/buildcuda)
    find_package(Kokkos REQUIRED)

    enable_language(CUDA)
    find_package(CUDA REQUIRED)

    # Specify the include directories
    include_directories(${CUDA_INCLUDE_DIRS})
    
    add_library(libpdemodelcuda SHARED ${PDEMODEL_SRC})
    target_link_libraries(libpdemodelcuda PRIVATE Kokkos::kokkos)
    set_target_properties(libpdemodelcuda PROPERTIES
        OUTPUT_NAME pdemodelcuda
        LIBRARY_OUTPUT_DIRECTORY ${EXASIM_PATH}/backend/Model
    )
elseif(EXASIM_HIP)
    if(EXASIM_TUOLUMNE)
      set(Kokkos_DIR ${EXASIM_PATH}/kokkos/buildtuolumne)
    else()
      set(Kokkos_DIR ${EXASIM_PATH}/kokkos/buildhip)
    endif()

    find_package(Kokkos REQUIRED)

    enable_language(HIP)
    find_package(HIP REQUIRED)

    include_directories(${HIP_INCLUDE_DIRS})

    add_library(libpdemodelhip SHARED ${PDEMODEL_SRC})
    target_link_libraries(libpdemodelhip PRIVATE Kokkos::kokkos)
    set_target_properties(libpdemodelhip PROPERTIES
        OUTPUT_NAME pdemodelhip
        LIBRARY_OUTPUT_DIRECTORY ${EXASIM_PATH}/backend/Model
    )
else()
    set(Kokkos_DIR ${EXASIM_PATH}/kokkos/buildserial)
    find_package(Kokkos REQUIRED)
    
    add_library(libpdemodelserial SHARED ${PDEMODEL_SRC})
    target_link_libraries(libpdemodelserial PRIVATE Kokkos::kokkos)
    set_target_properties(libpdemodelserial PROPERTIES
        OUTPUT_NAME pdemodelserial
        LIBRARY_OUTPUT_DIRECTORY ${EXASIM_PATH}/backend/Model
    )
endif()
