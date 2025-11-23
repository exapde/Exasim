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

.PHONY: text2code

text2code:
	@echo "Building with -j$(NPROCS)…"

	# GKlib (plain make, but still keep env clean for safety)
	cd GKlib && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make config prefix="$$PWD" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make install

	# METIS: requires 'make config'; run in a clean environment to avoid the 'w' target
	cd METIS && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make -j"$(NPROCS)" install

	# ParMETIS: requires 'make config'; run in a clean environment to avoid the 'w' target
	cd ParMETIS && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make config prefix="$$PWD" gklib_path="$$PWD/../GKlib" metis_path="$$PWD/../METIS" && \
	  env -i PATH="$$PATH" SHELL="$$SHELL" make -j"$(NPROCS)" install

	# SymEngine (CMake is fine)
	cd symengine && cmake -S . -B build \
	  -DCMAKE_INSTALL_PREFIX="$$PWD" \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_SHARED_LIBS=OFF \
	  -DWITH_GMP=OFF -DWITH_MPFR=OFF -DWITH_LLVM=OFF \
	  -DWITH_SYMENGINE_THREAD_SAFE=ON \
	  -DINTEGER_CLASS=boostmp \
	  -DCMAKE_CXX_FLAGS="-std=gnu++17 -fPIC" && \
	  cmake --build build -j"$(NPROCS)" && \
	  cmake --install build

	# text2code
	cd text2code && cmake -S . -B build -DWITH_METIS=ON && \
	  cmake --build build -j"$(NPROCS)"