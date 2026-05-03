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
    # Compile the templated `libpdemodel.cpp` directly into each binary
    # rather than linking it as a separately-loaded shared library. The
    # in-tree harness historically links `libpdemodel{serial,cuda,hip}.so`
    # (built by `text2code` under `backend/Model/`); on GPU that ends up
    # with Kokkos statically present in BOTH the main exe and
    # `libpdemodel*.so`, and Kokkos's static class members
    # (e.g. `SerialInternal::all_instances`,
    # `CudaInternal::constantMemHostStagingPerDevice`) get duplicate
    # `__cxa_atexit` destructor registrations — one per DSO. At program
    # exit both fire and double-free the same heap storage, aborting the
    # process with `free(): double free detected in tcache 2` (or
    # `Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()`).
    #
    # The compute output is correct — the crash is purely on shutdown —
    # but it makes every tutorial binary exit non-zero on GPU.
    # Compiling `libpdemodel.cpp` as a TU of the main binary keeps the
    # process single-DSO with respect to Kokkos, so the destructor is
    # registered exactly once. The math, the per-rank `outudg_np*.bin`
    # files, and the `max|udg|` value are all unchanged; the binary
    # simply exits cleanly.
    set(_pdemodel_src "${EXASIM_DIR}/backend/Model/libpdemodel.cpp")
    add_executable(${target} ${source_path} ${_pdemodel_src})
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_compile_definitions(${target} PRIVATE _TEXT2CODE)
    get_filename_component(_section_dir "${source_path}" DIRECTORY)
    target_include_directories(${target} PRIVATE
        ${EXASIM_DIR}/include
        "${EXASIM_DIR}/backend/Model"
        "${_section_dir}")

    if(EXASIM_MPI AND MPI_FOUND)
        target_compile_definitions(${target} PRIVATE _MPI)
        target_link_libraries(${target} PRIVATE ${MPI_LIBRARIES})
    endif()

    if(EXASIM_CUDA)
        target_compile_definitions(${target} PRIVATE _CUDA)
        target_link_libraries(${target} PRIVATE
            ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
        # Compile libpdemodel.cpp as a CUDA TU so its Kokkos kernels
        # use the same backend the rest of the binary does. Match
        # the flags text2code uses when it builds libpdemodelcuda.so
        # standalone — Kokkos `KOKKOS_INLINE_FUNCTION` lambdas need
        # nvcc's extended-lambda + relaxed-constexpr.
        set_source_files_properties(${_pdemodel_src} PROPERTIES
            LANGUAGE CUDA
            COMPILE_OPTIONS "--expt-extended-lambda;--expt-relaxed-constexpr")
    elseif(EXASIM_HIP)
        target_compile_definitions(${target} PRIVATE _HIP)
        set_source_files_properties(${_pdemodel_src} PROPERTIES
            LANGUAGE HIP)
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
# Both need text2code to have emitted `my_model.hpp` into the
# section dir before configure, since main.cpp #includes it. We
# add the targets only when that artifact is present, which
# `tutorial/run-all.sh` ensures by running text2code first.
foreach(_pair IN ITEMS
    "02-generated-cli|tutorial_02_generated_cli"
    "03-generated-embedded|tutorial_03_generated_embedded")
    string(REPLACE "|" ";" _pair_list "${_pair}")
    list(GET _pair_list 0 _dir)
    list(GET _pair_list 1 _target)
    set(_section "${CMAKE_CURRENT_SOURCE_DIR}/${_dir}")
    if(EXISTS "${_section}/my_model.hpp")
        tutorial_add_target(${_target} "${_section}/main.cpp")
    else()
        message(STATUS "tutorial: ${_dir} skipped (run text2code --out-dir <path> first)")
    endif()
endforeach()

# --- Section 04 — handwritten + CLI -----------------------------
tutorial_add_target(tutorial_04_handwritten_cli
    "${CMAKE_CURRENT_SOURCE_DIR}/04-handwritten-cli/main.cpp")

# --- Section 05 — handwritten + embedded (single rank only) -----
# Skipped on MPI builds; section 06 is the distributed counterpart.
if(NOT EXASIM_MPI)
    tutorial_add_target(tutorial_05_handwritten_embedded
        "${CMAKE_CURRENT_SOURCE_DIR}/05-handwritten-embedded/main.cpp")
endif()

# --- Section 06 — handwritten + distributed ---------------------
# Builds and links on every backend variant; meaningful primarily
# under `mpirun -np N` on MPI-enabled builds.
tutorial_add_target(tutorial_06_handwritten_distributed
    "${CMAKE_CURRENT_SOURCE_DIR}/06-handwritten-distributed/main.cpp")
