cmake_minimum_required(VERSION 3.4)


project(OrchestrAR)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

find_package(OpenCV REQUIRED)
find_package(SDL2 REQUIRED)
find_package(SDL2_mixer REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)

set(ORCHESTRAR_SOURCES src/Main.cpp)


# set(CPCT_LIB ${CMAKE_CURRENT_SOURCE_DIR}/include/cpct_library/cpct_dll/Debug/cpct_dll.lib 
 #    ${CMAKE_CURRENT_SOURCE_DIR}/include/cpct_library/cpct-mstftm/Debug/cpct-mstftm.lib)

include_directories(${OpenCV_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} ${SDL2_INCLUDE_DIR} ${SDL2_MIXER_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/cpct_library/cpct_dll/cpct_dll)   
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/cpct_library/cpct-mstftm/cpct-mstftm)
include_directories( ${GLEW_INCLUDE_DIRS} "/usr/local/Cellar/glfw/3.3.7/include") 


add_executable(OrchestrAR ${ORCHESTRAR_SOURCES})
target_link_libraries (OrchestrAR ${OpenCV_LIBS} ${SDL2_LIBRARY} ${SDL2MAIN_LIBRARY} ${SDL2_MIXER_LIBRARY} ${CPCT_LIB} ${GLEW_LIBRARIES} "/usr/local/Cellar/glfw/3.3.7/lib/libglfw.3.3.dylib" ${OPENGL_LIBRARIES})
