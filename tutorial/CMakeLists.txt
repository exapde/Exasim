# Tutorial in-tree wiring.
#
# Each section subdirectory has its own `CMakeLists.txt` that is
# the consumer-style example a real out-of-tree user would write
# (find_package(Exasim REQUIRED) + add_executable +
# target_link_libraries). Those files are NOT consumed by the
# in-tree build — they exist as documentation and as standalone
# build entry points for external users.
#
# This file is the in-tree wiring. It registers each tutorial
# target directly using the same per-target pattern that
# install/CMakeLists.txt uses for the in-tree apps
# (poisson2d_template_gpu, solve_square_facade_gpu, etc.). It
# defines _CUDA / _HIP / _MPI for the active build variant and
# links the appropriate pdemodel{serial,cuda,hip} library — none
# of which a real external consumer needs to deal with on their
# own (the install ships Exasim::headers and the consumer's
# CMakeLists picks up the rest via find_package).
#
# Sections are ordered by decreasing automation:
#   00 — overview and naming convention (no binary)
#   01 — generated model + prebuilt solver (drop two text files)
#   02 — generated model + CLI wrapper (3-line main, header-only)
#   03 — generated model + embedded solver (load_pdeapp from C++)
#   04 — handwritten model + CLI wrapper
#   05 — handwritten model + embedded solver (in-memory mesh,
#                                             single rank only)
#   06 — handwritten model + embedded solver (set_mesh_distributed)
#
# The generated-model sections (02 and 03) require text2code-emitted
# files to exist in their source directories before `cmake --build`
# can link them; `tutorial/run-all.sh` runs `text2code` before
# invoking the build.

cmake_minimum_required(VERSION 3.16)

# Helper: configure a tutorial target with the right backend
# defines and libraries for the active in-tree build variant.
#
# `<exasim/run.hpp>` translates `_CUDA` / `_HIP` / `_MPI` into
# `HAVE_CUDA` / `HAVE_HIP` / `HAVE_MPI`. Without these defines the
# code falls back to CPU paths even on a GPU build dir, which
# causes the host BLAS to read host pointers that alias unmapped
# device memory.
function(tutorial_add_target target source_path)
    add_executable(${target} ${source_path})
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_compile_definitions(${target} PRIVATE _TEXT2CODE)
    get_filename_component(_section_dir "${source_path}" DIRECTORY)
    target_include_directories(${target} PRIVATE
        ${EXASIM_DIR}/include
        "${_section_dir}")

    if(EXASIM_MPI AND MPI_FOUND)
        target_compile_definitions(${target} PRIVATE _MPI)
        target_link_libraries(${target} PRIVATE ${MPI_LIBRARIES})
    endif()

    if(EXASIM_CUDA)
        target_compile_definitions(${target} PRIVATE _CUDA)
        target_link_libraries(${target} PRIVATE
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${T2C_CUDA_LIB})
    elseif(EXASIM_HIP)
        target_compile_definitions(${target} PRIVATE _HIP)
        target_link_libraries(${target} PRIVATE ${T2C_HIP_LIB})
    else()
        target_link_libraries(${target} PRIVATE ${T2C_CPU_LIB})
    endif()

    target_link_libraries(${target} PRIVATE
        Kokkos::kokkos ${lapackblas_libraries})
    link_metis(${target})
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
endfunction()

# Section 00 — overview only; no binary to build.

# Section 01 — generated model with prebuilt solver. The binary is
# `cput2cEXASIM` (and friends), built by install/CMakeLists.txt
# under WITH_TEXT2CODE=ON. The tutorial directory holds only inputs.

# --- Sections 02 and 03 — generated + (CLI | embedded) ----------
# Both require text2code-emitted files (my_model.hpp +
# libpdemodelserial.{so,dylib}) to exist in the section dir at
# link time. We add the targets only when those artifacts are
# present, which tutorial/run-all.sh ensures by running text2code
# first.
foreach(_pair IN ITEMS
    "02-generated-cli|tutorial_02_generated_cli"
    "03-generated-embedded|tutorial_03_generated_embedded")
    string(REPLACE "|" ";" _pair_list "${_pair}")
    list(GET _pair_list 0 _dir)
    list(GET _pair_list 1 _target)
    set(_section "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}")
    if(EXISTS "${_section}/my_model.hpp")
        tutorial_add_target(${_target} "${_section}/main.cpp")
        # The codegen library lives in the section's source dir.
        target_link_directories(${_target} PRIVATE "${_section}")
        set_target_properties(${_target} PROPERTIES
            BUILD_RPATH "${_section}")
    else()
        message(STATUS "tutorial: ${_dir} skipped (run text2code --out-dir <path> first)")
    endif()
endforeach()

# --- Section 04 — handwritten + CLI -----------------------------
tutorial_add_target(tutorial_04_handwritten_cli
    "${CMAKE_CURRENT_SOURCE_DIR}/04-handwritten-cli/main.cpp")
target_link_directories(tutorial_04_handwritten_cli PRIVATE ${Model_LIB_DIR})
set_target_properties(tutorial_04_handwritten_cli PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")

# --- Section 05 — handwritten + embedded (single rank only) -----
# Skipped on MPI builds; section 06 is the distributed counterpart.
if(NOT EXASIM_MPI)
    tutorial_add_target(tutorial_05_handwritten_embedded
        "${CMAKE_CURRENT_SOURCE_DIR}/05-handwritten-embedded/main.cpp")
    target_link_directories(tutorial_05_handwritten_embedded PRIVATE ${Model_LIB_DIR})
    set_target_properties(tutorial_05_handwritten_embedded PROPERTIES
        BUILD_RPATH "${Model_LIB_DIR}")
endif()

# --- Section 06 — handwritten + distributed ---------------------
# Builds and links on every backend variant; meaningful primarily
# under `mpirun -np N` on MPI-enabled builds.
tutorial_add_target(tutorial_06_handwritten_distributed
    "${CMAKE_CURRENT_SOURCE_DIR}/06-handwritten-distributed/main.cpp")
target_link_directories(tutorial_06_handwritten_distributed PRIVATE ${Model_LIB_DIR})
set_target_properties(tutorial_06_handwritten_distributed PROPERTIES
    BUILD_RPATH "${Model_LIB_DIR}")
