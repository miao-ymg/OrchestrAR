cmake_minimum_required(VERSION 3.4)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

project(OrchestrAR)

find_package(OpenCV REQUIRED)
find_package(SDL2 REQUIRED)
find_package(SDL2_mixer REQUIRED)

# message(STATUS ${SDL2_INCLUDE_DIR})
# message(STATUS ${SDL2_MIXER_INCLUDE_DIRS})
# message(STATUS ${SDL2_MIXER_LIBRARY})
include_directories( ${OpenCV_INCLUDE_DIRS})
include_directories( ${SDL2_INCLUDE_DIR})
include_directories( ${SDL2_MIXER_INCLUDE_DIR})

set(EXECUTABLE_NAME OrchestrAR.cpp)
        
add_executable(OrchestrAR ${EXECUTABLE_NAME})

target_link_libraries (OrchestrAR ${OpenCV_LIBS} ${SDL2_LIBRARY} ${SDL2MAIN_LIBRARY} ${SDL2_MIXER_LIBRARY})
