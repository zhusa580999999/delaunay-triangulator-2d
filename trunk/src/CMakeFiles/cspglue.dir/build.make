# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.4

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jack/triangle

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jack/triangle

# Include any dependencies generated for this target.
include src/CMakeFiles/cspglue.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/cspglue.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/cspglue.dir/flags.make

src/CMakeFiles/cspglue.dir/depend.make.mark: src/CMakeFiles/cspglue.dir/flags.make
src/CMakeFiles/cspglue.dir/depend.make.mark: src/csp_glue.cpp

src/CMakeFiles/cspglue.dir/csp_glue.o: src/CMakeFiles/cspglue.dir/flags.make
src/CMakeFiles/cspglue.dir/csp_glue.o: src/csp_glue.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jack/triangle/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/cspglue.dir/csp_glue.o"
	/usr/bin/c++   $(CXX_FLAGS) -o src/CMakeFiles/cspglue.dir/csp_glue.o -c /home/jack/triangle/src/csp_glue.cpp

src/CMakeFiles/cspglue.dir/csp_glue.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to src/CMakeFiles/cspglue.dir/csp_glue.i"
	/usr/bin/c++  $(CXX_FLAGS) -E /home/jack/triangle/src/csp_glue.cpp > src/CMakeFiles/cspglue.dir/csp_glue.i

src/CMakeFiles/cspglue.dir/csp_glue.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly src/CMakeFiles/cspglue.dir/csp_glue.s"
	/usr/bin/c++  $(CXX_FLAGS) -S /home/jack/triangle/src/csp_glue.cpp -o src/CMakeFiles/cspglue.dir/csp_glue.s

src/CMakeFiles/cspglue.dir/csp_glue.o.requires:

src/CMakeFiles/cspglue.dir/csp_glue.o.provides: src/CMakeFiles/cspglue.dir/csp_glue.o.requires
	$(MAKE) -f src/CMakeFiles/cspglue.dir/build.make src/CMakeFiles/cspglue.dir/csp_glue.o.provides.build

src/CMakeFiles/cspglue.dir/csp_glue.o.provides.build: src/CMakeFiles/cspglue.dir/csp_glue.o

src/CMakeFiles/cspglue.dir/depend: src/CMakeFiles/cspglue.dir/depend.make.mark

src/CMakeFiles/cspglue.dir/depend.make.mark:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --magenta --bold "Scanning dependencies of target cspglue"
	cd /home/jack/triangle && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jack/triangle /home/jack/triangle/src /home/jack/triangle /home/jack/triangle/src /home/jack/triangle/src/CMakeFiles/cspglue.dir/DependInfo.cmake

# Object files for target cspglue
cspglue_OBJECTS = \
"CMakeFiles/cspglue.dir/csp_glue.o"

# External object files for target cspglue
cspglue_EXTERNAL_OBJECTS =

lib/libcspglue.so: src/CMakeFiles/cspglue.dir/csp_glue.o
lib/libcspglue.so: lib/libreader.so
lib/libcspglue.so: lib/libwriter.so
lib/libcspglue.so: src/CMakeFiles/cspglue.dir/build.make
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../lib/libcspglue.so"
	cd /home/jack/triangle/src && $(CMAKE_COMMAND) -P CMakeFiles/cspglue.dir/cmake_clean_target.cmake
	cd /home/jack/triangle/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cspglue.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/cspglue.dir/build: lib/libcspglue.so

src/CMakeFiles/cspglue.dir/requires: src/CMakeFiles/cspglue.dir/csp_glue.o.requires

src/CMakeFiles/cspglue.dir/clean:
	cd /home/jack/triangle/src && $(CMAKE_COMMAND) -P CMakeFiles/cspglue.dir/cmake_clean.cmake

