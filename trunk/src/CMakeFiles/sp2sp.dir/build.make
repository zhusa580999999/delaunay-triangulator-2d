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
include src/CMakeFiles/sp2sp.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/sp2sp.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/sp2sp.dir/flags.make

src/CMakeFiles/sp2sp.dir/depend.make.mark: src/CMakeFiles/sp2sp.dir/flags.make
src/CMakeFiles/sp2sp.dir/depend.make.mark: src/sp2sp.cpp

src/CMakeFiles/sp2sp.dir/sp2sp.o: src/CMakeFiles/sp2sp.dir/flags.make
src/CMakeFiles/sp2sp.dir/sp2sp.o: src/sp2sp.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/jack/triangle/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/sp2sp.dir/sp2sp.o"
	/usr/bin/c++   $(CXX_FLAGS) -o src/CMakeFiles/sp2sp.dir/sp2sp.o -c /home/jack/triangle/src/sp2sp.cpp

src/CMakeFiles/sp2sp.dir/sp2sp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to src/CMakeFiles/sp2sp.dir/sp2sp.i"
	/usr/bin/c++  $(CXX_FLAGS) -E /home/jack/triangle/src/sp2sp.cpp > src/CMakeFiles/sp2sp.dir/sp2sp.i

src/CMakeFiles/sp2sp.dir/sp2sp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly src/CMakeFiles/sp2sp.dir/sp2sp.s"
	/usr/bin/c++  $(CXX_FLAGS) -S /home/jack/triangle/src/sp2sp.cpp -o src/CMakeFiles/sp2sp.dir/sp2sp.s

src/CMakeFiles/sp2sp.dir/sp2sp.o.requires:

src/CMakeFiles/sp2sp.dir/sp2sp.o.provides: src/CMakeFiles/sp2sp.dir/sp2sp.o.requires
	$(MAKE) -f src/CMakeFiles/sp2sp.dir/build.make src/CMakeFiles/sp2sp.dir/sp2sp.o.provides.build

src/CMakeFiles/sp2sp.dir/sp2sp.o.provides.build: src/CMakeFiles/sp2sp.dir/sp2sp.o

src/CMakeFiles/sp2sp.dir/depend: src/CMakeFiles/sp2sp.dir/depend.make.mark

src/CMakeFiles/sp2sp.dir/depend.make.mark:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --magenta --bold "Scanning dependencies of target sp2sp"
	cd /home/jack/triangle && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jack/triangle /home/jack/triangle/src /home/jack/triangle /home/jack/triangle/src /home/jack/triangle/src/CMakeFiles/sp2sp.dir/DependInfo.cmake

# Object files for target sp2sp
sp2sp_OBJECTS = \
"CMakeFiles/sp2sp.dir/sp2sp.o"

# External object files for target sp2sp
sp2sp_EXTERNAL_OBJECTS =

bin/sp2sp: src/CMakeFiles/sp2sp.dir/sp2sp.o
bin/sp2sp: lib/libreader.so
bin/sp2sp: lib/libwriter.so
bin/sp2sp: src/CMakeFiles/sp2sp.dir/build.make
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/sp2sp"
	cd /home/jack/triangle/src && $(CMAKE_COMMAND) -P CMakeFiles/sp2sp.dir/cmake_clean_target.cmake
	cd /home/jack/triangle/src && /usr/bin/c++     -g  -fPIC $(sp2sp_OBJECTS) $(sp2sp_EXTERNAL_OBJECTS)  -o ../bin/sp2sp -rdynamic -L/home/jack/triangle/lib -lreader -lwriter -Wl,-rpath,/home/jack/triangle/lib 

# Rule to build all files generated by this target.
src/CMakeFiles/sp2sp.dir/build: bin/sp2sp

src/CMakeFiles/sp2sp.dir/requires: src/CMakeFiles/sp2sp.dir/sp2sp.o.requires

src/CMakeFiles/sp2sp.dir/clean:
	cd /home/jack/triangle/src && $(CMAKE_COMMAND) -P CMakeFiles/sp2sp.dir/cmake_clean.cmake

