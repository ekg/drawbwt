# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(seqwish)
# We build using c++11
set(CMAKE_CXX_STANDARD 11)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -g")

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# taywee's C++ args library, header only
ExternalProject_Add(tayweeargs
  GIT_REPOSITORY "https://github.com/Taywee/args.git"
  GIT_TAG "3de44ec671db452cc0c4ef86399b108939768abb"
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")
ExternalProject_Get_property(tayweeargs SOURCE_DIR)
set(tayweeargs_INCLUDE "${SOURCE_DIR}")

set(CMAKE_BUILD_TYPE Debug) 

# set up our target executable and specify its dependencies and includes
add_executable(drawbwt
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  )
add_dependencies(drawbwt tayweeargs)
target_include_directories(drawbwt PUBLIC
  "${tayweeargs_INCLUDE}")
target_link_libraries(drawbwt
  cairo
  pthread
  z)
#set(CMAKE_EXE_LINKER_FLAGS "-static")
