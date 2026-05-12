cmake_minimum_required(VERSION 3.16)
project(text2code_builder LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Options & paths ----
option(WITH_METIS "Build text2code with METIS/GKlib" ON)

# Path to the Exasim root (so we install into <EXASIM_ROOT>/build)
# Default guesses two levels up from this CMakeLists.txt; override if needed.
set(EXASIM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../" CACHE PATH "Path to Exasim root")

# Put executables directly under Exasim/build at build time
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${EXASIM_ROOT}/build")

# Adjust these if your layout differs; they default to your examples.
#set(METIS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../METIS" CACHE PATH "Path to METIS root")
#set(GKLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../GKlib" CACHE PATH "Path to GKlib root")
set(METIS_ROOT    "${EXASIM_ROOT}/metis/METIS" CACHE PATH "Path to METIS root")
set(GKLIB_ROOT    "${EXASIM_ROOT}/metis/GKlib" CACHE PATH "Path to GKlib root")

# Derived include/lib dirs (override if needed)
set(METIS_INCLUDE_DIR "${METIS_ROOT}/build/xinclude" CACHE PATH "METIS public include dir")
set(METIS_LIBRARY_DIR "${METIS_ROOT}/build/libmetis"  CACHE PATH "METIS library dir")
set(GKLIB_INCLUDE_DIR "${GKLIB_ROOT}/include"         CACHE PATH "GKlib include dir")
set(GKLIB_LIBRARY_DIR "${GKLIB_ROOT}/lib"             CACHE PATH "GKlib library dir")

# ---- BLAS / LAPACK ----
# Try to find proper BLAS/LAPACK; fall back to -lblas -llapack if not found.
find_package(BLAS)
find_package(LAPACK)
if (NOT BLAS_FOUND)
  message(WARNING "CMake could not find BLAS; will link '-lblas' explicitly.")
endif()
if (NOT LAPACK_FOUND)
  message(WARNING "CMake could not find LAPACK; will link '-llapack' explicitly.")
endif()

# Helper to add -O2 (or /O2 on MSVC)
function(add_o2 target)
  if (MSVC)
    target_compile_options(${target} PRIVATE /O2)
  else()
    target_compile_options(${target} PRIVATE -O2)
  endif()
endfunction()

# Helper to link BLAS/LAPACK (found targets or raw names)
function(link_blas_lapack target)
  if (BLAS_FOUND)
    target_link_libraries(${target} PRIVATE BLAS::BLAS)
  else()
    target_link_libraries(${target} PRIVATE blas)
  endif()
  if (LAPACK_FOUND)
    target_link_libraries(${target} PRIVATE LAPACK::LAPACK)
  else()
    target_link_libraries(${target} PRIVATE lapack)
  endif()
endfunction()

set(MODEL_FILE "${CMAKE_CURRENT_SOURCE_DIR}/pdeapp.txt" CACHE FILEPATH "Input for text2code")

# ============================================================
# Target: text2code (with METIS)  (equivalent to your METIS line)
# Only if WITH_METIS=ON and the expected dirs exist
# ============================================================
if (WITH_METIS)
  if (EXISTS "${METIS_INCLUDE_DIR}" AND EXISTS "${METIS_LIBRARY_DIR}"
      AND EXISTS "${GKLIB_INCLUDE_DIR}" AND EXISTS "${GKLIB_LIBRARY_DIR}")

    add_executable(text2code text2code.cpp)
    add_o2(text2code)

    target_compile_definitions(text2code PRIVATE HAVE_METIS)

    target_include_directories(text2code PRIVATE
      "${METIS_INCLUDE_DIR}"
      "${GKLIB_INCLUDE_DIR}"
    )

    # Prefer modern per-target link directories
    target_link_directories(text2code PRIVATE
      "${METIS_LIBRARY_DIR}"
      "${GKLIB_LIBRARY_DIR}"
    )

    # Link order: your command links metis, GKlib, then BLAS/LAPACK
    target_link_libraries(text2code PRIVATE
      metis
      GKlib
    )
    link_blas_lapack(text2code)

    # Auto-run after the executable links
    if (EXISTS "${MODEL_FILE}")
      add_custom_command(
        TARGET text2code POST_BUILD
        COMMAND $<TARGET_FILE:text2code> "${MODEL_FILE}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"   # or wherever the program expects to run
        COMMENT "Running text2code on ${MODEL_FILE}"
      )
    endif()

    # -------- t2cusecmake (HAVE_METIS + USE_CMAKE) --------
    add_executable(t2cusecmake text2code.cpp)
    add_o2(t2cusecmake)
    target_compile_definitions(t2cusecmake PRIVATE HAVE_METIS USE_CMAKE)

    target_include_directories(t2cusecmake PRIVATE
      "${METIS_INCLUDE_DIR}"
      "${GKLIB_INCLUDE_DIR}"
    )
    target_link_directories(t2cusecmake PRIVATE
      "${METIS_LIBRARY_DIR}"
      "${GKLIB_LIBRARY_DIR}"
    )
    target_link_libraries(t2cusecmake PRIVATE
      metis
      GKlib
    )
    link_blas_lapack(t2cusecmake)

  else()
    message(WARNING
      "WITH_METIS=ON but one or more METIS/GKlib dirs do not exist:\n"
      "  METIS_INCLUDE_DIR='${METIS_INCLUDE_DIR}'\n"
      "  METIS_LIBRARY_DIR='${METIS_LIBRARY_DIR}'\n"
      "  GKLIB_INCLUDE_DIR='${GKLIB_INCLUDE_DIR}'\n"
      "  GKLIB_LIBRARY_DIR='${GKLIB_LIBRARY_DIR}'\n"
      "Skipping 'text2code' target.")
  endif()
else()
  # ============================================================
  # Target: t2cnometis  (equivalent to: g++ -O2 -std=c++17 text2code.cpp -o t2cnometis -lblas -llapack)
  # ============================================================
  add_executable(t2cnometis text2code.cpp)
  add_o2(t2cnometis)
  link_blas_lapack(t2cnometis)

  # Auto-run after the executable links
  if (EXISTS "${MODEL_FILE}")
    add_custom_command(
      TARGET t2cnometis POST_BUILD
      COMMAND $<TARGET_FILE:t2cnometis> "${MODEL_FILE}"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      COMMENT "Running t2cnometis on ${MODEL_FILE}"
    )
  endif()

  # -------- t2cnometisusecmake (USE_CMAKE) --------
  add_executable(t2cnometisusecmake text2code.cpp)
  add_o2(t2cnometisusecmake)
  target_compile_definitions(t2cnometisusecmake PRIVATE USE_CMAKE)
  link_blas_lapack(t2cnometisusecmake)
endif()

# ======================================================================
#  Summary of text2code configuration
# ======================================================================
message(STATUS "==================== text2code CONFIGURATION ==================")

message(STATUS "EXASIM_ROOT              = ${EXASIM_ROOT}")
message(STATUS "WITH_METIS               = ${WITH_METIS}")

message(STATUS "METIS_ROOT               = ${METIS_ROOT}")
message(STATUS "GKLIB_ROOT               = ${GKLIB_ROOT}")

message(STATUS "METIS_INCLUDE_DIR        = ${METIS_INCLUDE_DIR}")
message(STATUS "METIS_LIBRARY_DIR        = ${METIS_LIBRARY_DIR}")
message(STATUS "GKLIB_INCLUDE_DIR        = ${GKLIB_INCLUDE_DIR}")
message(STATUS "GKLIB_LIBRARY_DIR        = ${GKLIB_LIBRARY_DIR}")

message(STATUS "MODEL_FILE               = ${MODEL_FILE}")

# BLAS / LAPACK
message(STATUS "BLAS_FOUND               = ${BLAS_FOUND}")
message(STATUS "LAPACK_FOUND             = ${LAPACK_FOUND}")
message(STATUS "BLAS_LIBRARIES           = ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LIBRARIES         = ${LAPACK_LIBRARIES}")

# C++ compiler
message(STATUS "C++ COMPILER             = ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ COMPILER ID          = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ COMPILER VERSION     = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE         = ${CMAKE_BUILD_TYPE}")

# Targets summary
if (WITH_METIS)
  if (TARGET text2code)
    message(STATUS "Target 'text2code'       = ENABLED (HAVE_METIS)")
  else()
    message(STATUS "Target 'text2code'       = SKIPPED (missing METIS/GKlib dirs)")
  endif()

  if (TARGET t2cusecmake)
    message(STATUS "Target 't2cusecmake'     = ENABLED (HAVE_METIS + USE_CMAKE)")
  else()
    message(STATUS "Target 't2cusecmake'     = NOT CREATED")
  endif()
else()
  if (TARGET t2cnometis)
    message(STATUS "Target 't2cnometis'      = ENABLED (no METIS)")
  else()
    message(STATUS "Target 't2cnometis'      = NOT CREATED")
  endif()

  if (TARGET t2cnometisusecmake)
    message(STATUS "Target 't2cnometisusecmake' = ENABLED (no METIS + USE_CMAKE)")
  else()
    message(STATUS "Target 't2cnometisusecmake' = NOT CREATED")
  endif()
endif()

message(STATUS "===============================================================")

# ----------------------
# Usage examples:
#   cmake -S . -B build
#   cmake --build build -j
#
# To point at custom installs:
#   cmake -S . -B build -DWITH_METIS=ON \
#         -DMETIS_ROOT=/path/to/METIS \
#         -DGKLIB_ROOT=/path/to/GKlib
#
# If CMake cannot find BLAS/LAPACK but you know they’re on link path, it will fall back
# to linking '-lblas -llapack' as in your original commands.
