cmake_minimum_required(VERSION 3.4)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

project(SDLTest)

find_package(OpenCV REQUIRED)
find_package(SDL2 REQUIRED)
find_package(SDL2_mixer REQUIRED)

message(STATUS ${SDL2_INCLUDE_DIR})
message(STATUS ${SDL2_MIXER_INCLUDE_DIRS})
message(STATUS ${SDL2_MIXER_LIBRARY})
include_directories( ${OpenCV_INCLUDE_DIRS})
include_directories( ${SDL2_INCLUDE_DIR})
include_directories( ${SDL2_MIXER_INCLUDE_DIRS})

# set(Exercise01_SOURCES
  #       Exercise-1.cpp)
set(SDL_SOURCES
        OrchestrAR.cpp)
        
add_executable(SDLTest ${SDL_SOURCES})
# add_executable(Exercise01 ${Exercise01_SOURCES})

# target_link_libraries (Exercise01 ${OpenCV_LIBS})
target_link_libraries (SDLTest ${OpenCV_LIBS} /Users/sam/Desktop/arwork/OrchestrAR/SDL2.framework/SDL2 /Users/sam/Desktop/arwork/OrchestrAR/SDL2_mixer.framework/SDL2_mixer)
