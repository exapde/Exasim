cmake_minimum_required(VERSION 3.16)
project(StaticLibraryABI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(core STATIC
    src/core.cpp
)
target_include_directories(core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(dynamic_provider SHARED
    providers/dynamic/dynamic_provider.cpp
)
target_include_directories(dynamic_provider PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
set_target_properties(dynamic_provider PROPERTIES
    OUTPUT_NAME dynamic_provider
)

add_executable(run_dynamic
    apps/run_dynamic.cpp
)
target_link_libraries(run_dynamic PRIVATE core dynamic_provider)
set_target_properties(run_dynamic PROPERTIES
    BUILD_RPATH "$<TARGET_FILE_DIR:dynamic_provider>"
)

add_executable(run_source
    apps/run_source.cpp
    providers/source/source_provider.cpp
)
target_link_libraries(run_source PRIVATE core)

add_executable(run_template
    apps/run_template.cpp
    providers/template/template_provider.cpp
)
target_link_libraries(run_template PRIVATE core)
