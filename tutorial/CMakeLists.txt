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
#   01 — codegen + AbiAdapter (drop two text files, run a prebuilt binary)
#   02 — codegen + legacy CLI (write a 3-line main, link header-only)
#   03 — codegen + ExasimSolver facade (load_pdeapp from C++)
#   04 — hand-written Model + legacy CLI
#   05 — hand-written Model + ExasimSolver facade (set_mesh)
#   06 — hand-written Model + ExasimSolver facade (set_mesh_distributed)
#
# The codegen sections (02 and 03) require text2code-emitted files
# to exist in their source directories before `cmake --build` can
# link them; `run-all.sh` runs `text2code` before invoking the
# build.

cmake_minimum_required(VERSION 3.16)

# Section 01 — codegen + AbiAdapter. The binary is `cput2cEXASIM`
# (and friends), built by the main install/CMakeLists.txt under
# WITH_TEXT2CODE=ON. The tutorial directory holds only inputs;
# nothing to build here.

# Sections 02 and 03 — codegen + (CLI | facade). Both require the
# codegen artifacts (my_model.hpp + libpdemodelserial.{so,dylib})
# to exist in the source directory at link time. We add them as
# subprojects only when the codegen artifacts are present, which
# tutorial/run-all.sh ensures by running text2code first.
foreach(_dir IN ITEMS 02-codegen-cli 03-codegen-facade)
    set(_codegen_hdr "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/my_model.hpp")
    if(EXISTS "${_codegen_hdr}")
        add_subdirectory("${_dir}")
    else()
        message(STATUS "tutorial: ${_dir} skipped (run text2code --out-dir <path> first)")
    endif()
endforeach()

# Section 04 — hand-written + legacy CLI. Available on every backend
# variant.
add_subdirectory(04-handwritten-cli)

# Section 05 — hand-written + ExasimSolver<M>::set_mesh. Available on
# every backend variant.
add_subdirectory(05-handwritten-facade)

# Section 06 — hand-written + ExasimSolver<M>::set_mesh_distributed.
# Builds and links on every backend variant; meaningful primarily on
# MPI-enabled builds.
add_subdirectory(06-handwritten-facade-mpi)
