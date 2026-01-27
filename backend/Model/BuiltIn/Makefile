# Makefile (Linux + macOS + Windows via MinGW/MSYS2/Cygwin)
# Builds:
#   - libbuiltinmodelserial.(so|dylib|dll)
#   - libbuiltinmodelcuda.(so|dylib|dll)
#   - libbuiltinmodelhip.(so|dylib|dll)

UNAME_S := $(shell uname -s)

# --- platform-specific shared-lib flags / naming ---
ifeq ($(UNAME_S),Darwin)
  SHLIB_EXT := dylib
  SHLIB_PREFIX := lib
  LDSHARED   ?= -dynamiclib
else ifneq (,$(findstring MINGW,$(UNAME_S)))
  SHLIB_EXT := dll
  SHLIB_PREFIX :=
  LDSHARED   ?= -shared
  IMPORT_LIB_FLAGS := -Wl,--out-implib,
else ifneq (,$(findstring MSYS,$(UNAME_S)))
  SHLIB_EXT := dll
  SHLIB_PREFIX :=
  LDSHARED   ?= -shared
  IMPORT_LIB_FLAGS := -Wl,--out-implib,
else ifneq (,$(findstring CYGWIN,$(UNAME_S)))
  SHLIB_EXT := dll
  SHLIB_PREFIX :=
  LDSHARED   ?= -shared
  IMPORT_LIB_FLAGS := -Wl,--out-implib,
else
  SHLIB_EXT := so
  SHLIB_PREFIX := lib
  LDSHARED   ?= -shared
endif

# --- toolchain ---
CXX      ?= g++
CXX_CUDA ?= $(CXX)
CXX_HIP  ?= $(CXX)

STD      ?= -std=c++17
PIC      ?= -fPIC
OPT      ?= -O3

# On Windows, -fPIC is typically unnecessary; keep it but allow disabling.
ifneq (,$(findstring dll,$(SHLIB_EXT)))
  PIC :=
endif

# --- Kokkos build roots (adjust to your tree) ---
KOKKOS_SERIAL ?= ../../../kokkos/buildserial
KOKKOS_CUDA   ?= ../../../kokkos/buildcuda
KOKKOS_HIP    ?= ../../../kokkos/buildhip

# --- sources ---
SRC := libbuiltinmodel.cpp model1/model.cpp model2/model.cpp 

# --- outputs ---
TARGET_SERIAL := $(SHLIB_PREFIX)builtinmodelserial.$(SHLIB_EXT)
TARGET_CUDA   := $(SHLIB_PREFIX)builtinmodelcuda.$(SHLIB_EXT)
TARGET_HIP    := $(SHLIB_PREFIX)builtinmodelhip.$(SHLIB_EXT)

# MinGW import libs (only meaningful on Windows)
IMPLIB_SERIAL := $(SHLIB_PREFIX)builtinmodelserial.dll.a
IMPLIB_CUDA   := $(SHLIB_PREFIX)builtinmodelcuda.dll.a
IMPLIB_HIP    := $(SHLIB_PREFIX)builtinmodelhip.dll.a

# --- object dirs (separate to avoid flag collisions) ---
OBJDIR_SERIAL := build/serial
OBJDIR_CUDA   := build/cuda
OBJDIR_HIP    := build/hip

OBJ_SERIAL := $(patsubst %.cpp,$(OBJDIR_SERIAL)/%.o,$(SRC))
OBJ_CUDA   := $(patsubst %.cpp,$(OBJDIR_CUDA)/%.o,$(SRC))
OBJ_HIP    := $(patsubst %.cpp,$(OBJDIR_HIP)/%.o,$(SRC))

# --- flags per backend ---
CXXFLAGS_COMMON ?= $(STD) $(PIC) $(OPT)

CXXFLAGS_SERIAL ?= $(CXXFLAGS_COMMON) -I$(KOKKOS_SERIAL)/include
CXXFLAGS_CUDA   ?= $(CXXFLAGS_COMMON) -I$(KOKKOS_CUDA)/include
CXXFLAGS_HIP    ?= $(CXXFLAGS_COMMON) -I$(KOKKOS_HIP)/include

# --- libraries ---
LIB_SERIAL ?= $(KOKKOS_SERIAL)/lib/libkokkoscore.a
LIB_CUDA   ?= $(KOKKOS_CUDA)/lib/libkokkoscore.a
LIB_HIP    ?= $(KOKKOS_HIP)/lib/libkokkoscore.a

# --- macOS install_name (dylib id) flags ---
ifeq ($(UNAME_S),Darwin)
  INSTALL_NAME_SERIAL := -Wl,-install_name,@rpath/$(TARGET_SERIAL)
  INSTALL_NAME_CUDA   := -Wl,-install_name,@rpath/$(TARGET_CUDA)
  INSTALL_NAME_HIP    := -Wl,-install_name,@rpath/$(TARGET_HIP)
else
  INSTALL_NAME_SERIAL :=
  INSTALL_NAME_CUDA   :=
  INSTALL_NAME_HIP    :=
endif

.PHONY: all serial cuda hip clean

all: serial cuda hip

serial: $(TARGET_SERIAL)
cuda:   $(TARGET_CUDA)
hip:    $(TARGET_HIP)

# --- link ---
ifeq ($(SHLIB_EXT),dll)

$(TARGET_SERIAL): $(OBJ_SERIAL)
	$(CXX) $(LDSHARED) $^ $(LIB_SERIAL) $(IMPORT_LIB_FLAGS)$(IMPLIB_SERIAL) -o $@

$(TARGET_CUDA): $(OBJ_CUDA)
	$(CXX_CUDA) $(LDSHARED) $^ $(LIB_CUDA) $(IMPORT_LIB_FLAGS)$(IMPLIB_CUDA) -o $@

$(TARGET_HIP): $(OBJ_HIP)
	$(CXX_HIP) $(LDSHARED) $^ $(LIB_HIP) $(IMPORT_LIB_FLAGS)$(IMPLIB_HIP) -o $@

else

$(TARGET_SERIAL): $(OBJ_SERIAL)
	$(CXX) $(LDSHARED) $(INSTALL_NAME_SERIAL) $^ $(LIB_SERIAL) -o $@

$(TARGET_CUDA): $(OBJ_CUDA)
	$(CXX_CUDA) $(LDSHARED) $(INSTALL_NAME_CUDA) $^ $(LIB_CUDA) -o $@

$(TARGET_HIP): $(OBJ_HIP)
	$(CXX_HIP) $(LDSHARED) $(INSTALL_NAME_HIP) $^ $(LIB_HIP) -o $@

endif

# --- compile ---
$(OBJDIR_SERIAL)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_SERIAL) -c $< -o $@

$(OBJDIR_CUDA)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX_CUDA) $(CXXFLAGS_CUDA) -c $< -o $@

$(OBJDIR_HIP)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX_HIP) $(CXXFLAGS_HIP) -c $< -o $@

clean:
	rm -rf build \
	  $(TARGET_SERIAL) $(TARGET_CUDA) $(TARGET_HIP) \
	  $(IMPLIB_SERIAL) $(IMPLIB_CUDA) $(IMPLIB_HIP)