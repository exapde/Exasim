# ---- CPU core detection (Linux/macOS/WSL/MSYS2) ----
NPROCS := $(strip $(shell \
  (command -v nproc >/dev/null 2>&1 && nproc) || \
  (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu) || \
  (command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN) || \
  echo 1))

SHELL := /usr/bin/env sh

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
	  env -i PATH="$$PATH" SHELL="$$SHELL" make config prefix="$$PWD" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make install

	# METIS: requires 'make config'; run in a clean environment to avoid the 'w' target
	cd METIS && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make -j"$(NPROCS)" install

  # ParMETIS: 
	cd ParMETIS && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" HOME="$$HOME" make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" metis_path="$$PWD/../METIS" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" HOME="$$HOME" make -j"$(NPROCS)" install

	# preprocess
	cd preprocessing && cmake -S . -B build && cmake --build build -j"$(NPROCS)"