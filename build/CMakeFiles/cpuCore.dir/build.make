# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.27.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.27.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ngoccuongnguyen/Exasim/install

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ngoccuongnguyen/Exasim/build

# Include any dependencies generated for this target.
include CMakeFiles/cpuCore.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cpuCore.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cpuCore.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cpuCore.dir/flags.make

CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o: CMakeFiles/cpuCore.dir/flags.make
CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o: /Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp
CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o: CMakeFiles/cpuCore.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/ngoccuongnguyen/Exasim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o -MF CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o.d -o CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o -c /Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp

CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp > CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.i

CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp -o CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.s

# Object files for target cpuCore
cpuCore_OBJECTS = \
"CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o"

# External object files for target cpuCore
cpuCore_EXTERNAL_OBJECTS =

libcpuCore.a: CMakeFiles/cpuCore.dir/Users/ngoccuongnguyen/Exasim/lib/opuCore.cpp.o
libcpuCore.a: CMakeFiles/cpuCore.dir/build.make
libcpuCore.a: CMakeFiles/cpuCore.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/ngoccuongnguyen/Exasim/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libcpuCore.a"
	$(CMAKE_COMMAND) -P CMakeFiles/cpuCore.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cpuCore.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cpuCore.dir/build: libcpuCore.a
.PHONY : CMakeFiles/cpuCore.dir/build

CMakeFiles/cpuCore.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cpuCore.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cpuCore.dir/clean

CMakeFiles/cpuCore.dir/depend:
	cd /Users/ngoccuongnguyen/Exasim/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ngoccuongnguyen/Exasim/install /Users/ngoccuongnguyen/Exasim/install /Users/ngoccuongnguyen/Exasim/build /Users/ngoccuongnguyen/Exasim/build /Users/ngoccuongnguyen/Exasim/build/CMakeFiles/cpuCore.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/cpuCore.dir/depend
