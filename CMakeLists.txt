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
include_directories( ${SDL2_MIXER_INCLUDE_DIRS})

set(EXECUTABLE_NAME OrchestrAR.cpp)
        
add_executable(OrchestrAR ${EXECUTABLE_NAME})

target_link_libraries (OrchestrAR ${OpenCV_LIBS} ${CMAKE_CURRENT_SOURCE_DIR}/SDL2.framework/SDL2 ${CMAKE_CURRENT_SOURCE_DIR}/SDL2_mixer.framework/SDL2_mixer)
