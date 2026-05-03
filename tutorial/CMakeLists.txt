# Tutorial subprojects.
#
# Each subdirectory is a small self-contained consumer of the
# Exasim header-only library. Including this file from the
# top-level CMakeLists registers the tutorial binaries as standard
# build targets. The actual run logic (regenerating codegen
# artifacts, invoking the binaries, checking output) lives in
# `tutorial/run-all.sh`, which is wrapped by the
# `tutorial:run-all:<variant>` ctest test.
#
# The codegen tutorials (paths 03 and 04) require text2code-emitted
# files to exist in their source directories before `cmake --build`
# can link them. The build script in `run-all.sh` runs `text2code
# --out-dir <path>` before invoking the build.

cmake_minimum_required(VERSION 3.16)

# Path 01 — hand-written + legacy CLI. Available on every backend
# variant.
add_subdirectory(01-handwritten-cli)

# Path 02 — hand-written + ExasimSolver<M>. Available on every
# backend variant.
add_subdirectory(02-handwritten-facade)

# Paths 03 and 04 — codegen + (CLI | facade). Both require the
# codegen artifacts (my_model.hpp + libpdemodelserial.{so,dylib})
# to exist in the source directory at link time. We add them as
# subprojects only when the codegen artifacts are present, which
# tutorial/run-all.sh ensures by running text2code first.
foreach(_dir IN ITEMS 03-codegen-cli 04-codegen-facade)
    set(_codegen_hdr "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/my_model.hpp")
    if(EXISTS "${_codegen_hdr}")
        add_subdirectory("${_dir}")
    else()
        message(STATUS "tutorial: ${_dir} skipped (run text2code --out-dir <path> first)")
    endif()
endforeach()

# Path 05 — codegen + AbiAdapter. The binary is `cput2cEXASIM` (and
# friends), built by the main install/CMakeLists.txt under
# WITH_TEXT2CODE=ON. The tutorial directory holds only inputs;
# nothing to build here.
