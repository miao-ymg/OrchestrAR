cmake_minimum_required(VERSION 3.4)


project(OrchestrAR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

set(SOURCES Instrument.cpp)
set(EXECUTABLE_NAME OrchestrAR.cpp)
set(OWN_LIBRARY ownlib)


find_package(OpenCV REQUIRED)
find_package(SDL2 REQUIRED)
find_package(SDL2_mixer REQUIRED)

add_library(${OWN_LIBRARY} ${SOURCES})
include_directories(${OpenCV_INCLUDE_DIRS} ${SDL2_INCLUDE_DIR} ${SDL2_MIXER_INCLUDE_DIR})
        
add_executable(OrchestrAR ${EXECUTABLE_NAME})
target_link_libraries (OrchestrAR ${OpenCV_LIBS} ${SDL2_LIBRARY} ${SDL2MAIN_LIBRARY} ${SDL2_MIXER_LIBRARY} ${OWN_LIBRARY})
