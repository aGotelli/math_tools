cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(math_tools LANGUAGES CXX VERSION 2.0)
set(ProjectVersion "v_2_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message(STATUS "${PROJECT_NAME} version ${PROJECT_VERSION}")


#   Find required packages
find_package(Eigen3 3.4 NO_MODULE REQUIRED)
find_package(benchmark REQUIRED)


#   Process the content of src, ie, the actual code
add_subdirectory(src)

#   Install the library using the default routine
include(cmake/installation_module.cmake)
