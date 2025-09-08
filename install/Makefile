# ---- Portable CPU core detection (Linux / macOS / Git Bash/MSYS2/WSL) ----
NPROCS := $(strip $(shell \
  (command -v nproc >/dev/null 2>&1 && nproc) || \
  (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu) || \
  (command -v getconf >/dev/null 2>&1 && getconf _NPROCESSORS_ONLN) || \
  (command -v cmd.exe >/dev/null 2>&1 && cmd.exe /C echo %NUMBER_OF_PROCESSORS%) || \
  echo 1))

SHELL := /usr/bin/env sh

.PHONY: text2code

text2code:
	@echo "Building with -j$(NPROCS)…"
	cd GKlib && make config prefix=$$(pwd) && make install
	cd METIS && make config prefix=$$(pwd) gklib_path="$$(pwd)/../GKlib" && make -j$(NPROCS) install
	cd symengine && cmake -S . -B build \
	  -DCMAKE_INSTALL_PREFIX=$$(pwd) \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DBUILD_SHARED_LIBS=OFF \
	  -DWITH_GMP=OFF \
	  -DWITH_MPFR=OFF \
	  -DWITH_LLVM=OFF \
	  -DWITH_SYMENGINE_THREAD_SAFE=ON \
	  -DINTEGER_CLASS=boostmp \
	  -DCMAKE_CXX_FLAGS="-std=gnu++17 -fPIC" && \
	  cmake --build build -j$(NPROCS) && \
	  cmake --install build
	cd text2code && cmake -S . -B build -DWITH_METIS=ON && cmake --build build -j$(NPROCS)