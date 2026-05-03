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
# The sections are ordered by decreasing automation:
#   00 — overview and naming convention (no binary)
#   01 — generated model + prebuilt solver (drop two text files)
#   02 — generated model + CLI wrapper (3-line main, header-only)
#   03 — generated model + embedded solver (load_pdeapp from C++)
#   04 — handwritten model + CLI wrapper
#   05 — handwritten model + embedded solver (set_mesh, in-memory)
#   06 — handwritten model + embedded solver (set_mesh_distributed)
#
# The generated-model sections (02 and 03) require text2code-emitted
# files to exist in their source directories before `cmake --build`
# can link them; `run-all.sh` runs `text2code` before invoking the
# build.

cmake_minimum_required(VERSION 3.16)

# Section 00 — overview only; no binary to build.

# Section 01 — generated model with prebuilt solver. The binary is
# `cput2cEXASIM` (and friends), built by the main install/CMakeLists
# under WITH_TEXT2CODE=ON. The tutorial directory holds only inputs.

# Sections 02 and 03 — generated model + (CLI | embedded). Both
# require the codegen artifacts (my_model.hpp +
# libpdemodelserial.{so,dylib}) to exist in the source directory at
# link time. We add them as subprojects only when the codegen
# artifacts are present, which tutorial/run-all.sh ensures by
# running text2code first.
foreach(_dir IN ITEMS 02-generated-cli 03-generated-embedded)
    set(_codegen_hdr "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/my_model.hpp")
    if(EXISTS "${_codegen_hdr}")
        add_subdirectory("${_dir}")
    else()
        message(STATUS "tutorial: ${_dir} skipped (run text2code --out-dir <path> first)")
    endif()
endforeach()

# Section 04 — handwritten model + CLI wrapper.
add_subdirectory(04-handwritten-cli)

# Section 05 — handwritten model + embedded solver (in-memory mesh).
add_subdirectory(05-handwritten-embedded)

# Section 06 — handwritten model + embedded solver (distributed
# mesh). Builds and links on every backend variant; meaningful
# primarily on MPI-enabled builds.
add_subdirectory(06-handwritten-distributed)
