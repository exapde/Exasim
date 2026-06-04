# ---- CPU core detection (Linux/macOS/WSL/MSYS2) ----
NPROCS := $(strip $(shell \
  (command -v nproc >/dev/null 2>&1 && nproc) || \
  (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu 2>/dev/null) || \
  (command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN 2>/dev/null) || \
  echo 1))

SHELL := /usr/bin/env sh

# Keep the sub-build environment clean without dropping variables required by
# CMake, Homebrew, and macOS compiler discovery.
CLEAN_ENV = env -i PATH="$$PATH" SHELL="$$SHELL" HOME="$$HOME" USER="$$USER" TMPDIR="$$TMPDIR" SDKROOT="$$SDKROOT" MACOSX_DEPLOYMENT_TARGET="$$MACOSX_DEPLOYMENT_TARGET"

# ---- stop leaking parent make flags into sub-makes ----
MAKEFLAGS :=
MFLAGS :=
GNUMAKEFLAGS :=
unexport MAKEFLAGS
unexport MFLAGS
unexport GNUMAKEFLAGS

.PHONY: metis

metis:
	@echo "Building with -j$(NPROCS)…"

	# GKlib (plain make, but still keep env clean for safety)
	cd GKlib && \
	  $(CLEAN_ENV) make config prefix="$$PWD" && \
	  $(CLEAN_ENV) make install

	# METIS: requires 'make config'; run in a clean environment to avoid the 'w' target
	cd METIS && \
	  $(CLEAN_ENV) make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" && \
	  $(CLEAN_ENV) make -j"$(NPROCS)" install

  # ParMETIS: 
	cd ParMETIS && \
	  $(CLEAN_ENV) make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" metis_path="$$PWD/../METIS" && \
	  $(CLEAN_ENV) make -j"$(NPROCS)" install

	# preprocess
	cd preprocessing && cmake -S . -B build && cmake --build build -j"$(NPROCS)"
