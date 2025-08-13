cmake_minimum_required(VERSION 3.16)
project(text2code_superbuild LANGUAGES NONE)  # superbuild only orchestrates

#cmake -S . -B superbuild
#cmake --build superbuild

include(ExternalProject)

# --- tools we need ---
find_program(MAKE_EXE NAMES gmake make REQUIRED)

# Root layout (everything is already checked out locally)
set(ROOT        "${CMAKE_SOURCE_DIR}")
set(GKLIB_SRC   "${ROOT}/GKlib")
set(METIS_SRC   "${ROOT}/METIS")
set(SYM_SRC     "${ROOT}/symengine")
set(T2C_SRC     "${ROOT}/text2code")

# Where you want the final text2code binaries to land (same as before)
# Change if you prefer a different staging dir.
set(EXASIM_OUT  "${ROOT}/../build")

# macOS SDK (only if on Apple)
if(APPLE)
  execute_process(COMMAND xcrun --show-sdk-path
                  OUTPUT_VARIABLE MACOS_SDK
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

# ------------------------------------------------------------
# GKlib  (equivalent to: make config prefix=$(pwd) && make install)
# ------------------------------------------------------------
ExternalProject_Add(GKlib_ep
  SOURCE_DIR   "${GKLIB_SRC}"
  BINARY_DIR   "${GKLIB_SRC}"        # build in source, like your commands
  CONFIGURE_COMMAND ${MAKE_EXE} config prefix=<SOURCE_DIR>
  BUILD_COMMAND     ${MAKE_EXE}
  INSTALL_COMMAND   ${MAKE_EXE} install
  LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
)

# Pick a job count
if(DEFINED ENV{CMAKE_BUILD_PARALLEL_LEVEL})
  set(EP_JOBS $ENV{CMAKE_BUILD_PARALLEL_LEVEL})
else()
  set(EP_JOBS 4)
endif()

# ------------------------------------------------------------
# METIS (needs GKlib path)
#   make config prefix=$(pwd) gklib_path="$(pwd)/../GKlib" && make install
# ------------------------------------------------------------
ExternalProject_Add(METIS_ep
  SOURCE_DIR   "${METIS_SRC}"
  BINARY_DIR   "${METIS_SRC}"   # or better: an out-of-source dir
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env --unset=MAKEFLAGS
                    ${MAKE_EXE} config prefix=<SOURCE_DIR> gklib_path=${GKLIB_SRC}
  BUILD_COMMAND     ${CMAKE_COMMAND} -E env --unset=MAKEFLAGS
                    ${MAKE_EXE} -j ${EP_JOBS}
  INSTALL_COMMAND   ${CMAKE_COMMAND} -E env --unset=MAKEFLAGS
                    ${MAKE_EXE} -j ${EP_JOBS} install
  DEPENDS GKlib_ep
  LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
  USES_TERMINAL_CONFIGURE 1
  USES_TERMINAL_BUILD 1
  USES_TERMINAL_INSTALL 1
)

# ------------------------------------------------------------
# SymEngine (CMake build with your flags; installs into its source dir)
# ------------------------------------------------------------
set(SYM_BUILD "${SYM_SRC}/build")
ExternalProject_Add(SymEngine_ep
  SOURCE_DIR  "${SYM_SRC}"
  BINARY_DIR  "${SYM_BUILD}"
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${SYM_SRC}
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_SHARED_LIBS=OFF
    -DWITH_GMP=OFF
    -DWITH_MPFR=OFF
    -DWITH_LLVM=OFF
    -DWITH_SYMENGINE_THREAD_SAFE=ON
    -DINTEGER_CLASS=boostmp
    $<$<BOOL:${APPLE}>:-DCMAKE_OSX_SYSROOT=${MACOS_SDK}>
    # Your original CXX flags (macOS specific setup)
    -DCMAKE_CXX_FLAGS=-std=gnu++17\ -fPIC\ $<$<BOOL:${APPLE}>:-isysroot\ ${MACOS_SDK}\ -nostdinc++\ -isystem\ ${MACOS_SDK}/usr/include/c++/v1\ -stdlib=libc++>
  BUILD_COMMAND     ${CMAKE_COMMAND} --build .
  INSTALL_COMMAND   ${CMAKE_COMMAND} --build . --target install
  LOG_CONFIGURE 1 LOG_BUILD 1 LOG_INSTALL 1
)

# ------------------------------------------------------------
# text2code (your CMake project that can use METIS/GKlib)
# We build it and then stage the two executables into EXASIM_OUT.
# ------------------------------------------------------------
set(T2C_BUILD "${T2C_SRC}/build")
ExternalProject_Add(text2code_ep
  SOURCE_DIR  "${T2C_SRC}"
  BINARY_DIR  "${T2C_BUILD}"
  CMAKE_ARGS
    -DWITH_METIS=ON
    #-DMETIS_ROOT=${METIS_SRC}
    #-DGKLIB_ROOT=${GKLIB_SRC}
    # If your text2code CMakeLists installs the binaries, you can also pass:
    # -DCMAKE_INSTALL_PREFIX=${EXASIM_OUT}/..
  BUILD_COMMAND   ${CMAKE_COMMAND} --build . --config Release -j
  INSTALL_COMMAND ""  # we'll stage manually below
  DEPENDS METIS_ep SymEngine_ep
  LOG_CONFIGURE 1 LOG_BUILD 1
)

# After text2code builds, copy the binaries where you want them (Exasim/build).
# Adjust names if your targets differ.
# ExternalProject_Add_Step(text2code_ep stage_binaries
#  COMMAND ${CMAKE_COMMAND} -E make_directory "${EXASIM_OUT}"
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${T2C_BUILD}/text2code"  "${EXASIM_OUT}/"
#  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${T2C_BUILD}/t2cnometis" "${EXASIM_OUT}/"
#  DEPENDEES build
#  ALWAYS 1
#  COMMENT "Staging text2code and t2cnometis into ${EXASIM_OUT}"
# )

# A convenient top-level target name:
add_custom_target(all_text2code_build
  DEPENDS GKlib_ep METIS_ep SymEngine_ep text2code_ep)
