cmake_minimum_required(VERSION 3.4)
project(OrchestrAR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(SDL2 REQUIRED)
find_package(SDL2_mixer REQUIRED)

set(ORCHESTRAR_SOURCES src/main.cpp)


include_directories(${OpenCV_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${SDL2_INCLUDE_DIR} ${SDL2_MIXER_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(OrchestrAR ${ORCHESTRAR_SOURCES})
target_link_libraries (OrchestrAR ${OpenCV_LIBS} ${OPENGL_LIBRARIES} ${GLEW_LIBRARY} ${GLEW_LIBRARIES} ${GLFW_LIBRARIES} ${SDL2_LIBRARY} ${SDL2MAIN_LIBRARY} ${SDL2_MIXER_LIBRARY})