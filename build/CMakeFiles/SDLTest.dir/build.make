# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.23.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.23.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/sam/Desktop/arwork/OrchestrAR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/sam/Desktop/arwork/OrchestrAR/build

# Include any dependencies generated for this target.
include CMakeFiles/SDLTest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/SDLTest.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/SDLTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SDLTest.dir/flags.make

CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o: CMakeFiles/SDLTest.dir/flags.make
CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o: ../OrchestrAR.cpp
CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o: CMakeFiles/SDLTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/sam/Desktop/arwork/OrchestrAR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o -MF CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o.d -o CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o -c /Users/sam/Desktop/arwork/OrchestrAR/OrchestrAR.cpp

CMakeFiles/SDLTest.dir/OrchestrAR.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SDLTest.dir/OrchestrAR.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/sam/Desktop/arwork/OrchestrAR/OrchestrAR.cpp > CMakeFiles/SDLTest.dir/OrchestrAR.cpp.i

CMakeFiles/SDLTest.dir/OrchestrAR.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SDLTest.dir/OrchestrAR.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/sam/Desktop/arwork/OrchestrAR/OrchestrAR.cpp -o CMakeFiles/SDLTest.dir/OrchestrAR.cpp.s

# Object files for target SDLTest
SDLTest_OBJECTS = \
"CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o"

# External object files for target SDLTest
SDLTest_EXTERNAL_OBJECTS =

SDLTest: CMakeFiles/SDLTest.dir/OrchestrAR.cpp.o
SDLTest: CMakeFiles/SDLTest.dir/build.make
SDLTest: /usr/local/lib/libopencv_dnn.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_highgui.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_ml.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_objdetect.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_shape.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_stitching.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_superres.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_videostab.3.4.14.dylib
SDLTest: ../SDL2.framework/SDL2
SDLTest: ../SDL2_mixer.framework/SDL2_mixer
SDLTest: /usr/local/lib/libopencv_calib3d.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_features2d.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_flann.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_photo.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_video.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_videoio.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_imgcodecs.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_imgproc.3.4.14.dylib
SDLTest: /usr/local/lib/libopencv_core.3.4.14.dylib
SDLTest: CMakeFiles/SDLTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/sam/Desktop/arwork/OrchestrAR/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SDLTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SDLTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SDLTest.dir/build: SDLTest
.PHONY : CMakeFiles/SDLTest.dir/build

CMakeFiles/SDLTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SDLTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SDLTest.dir/clean

CMakeFiles/SDLTest.dir/depend:
	cd /Users/sam/Desktop/arwork/OrchestrAR/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/sam/Desktop/arwork/OrchestrAR /Users/sam/Desktop/arwork/OrchestrAR /Users/sam/Desktop/arwork/OrchestrAR/build /Users/sam/Desktop/arwork/OrchestrAR/build /Users/sam/Desktop/arwork/OrchestrAR/build/CMakeFiles/SDLTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SDLTest.dir/depend

