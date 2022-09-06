cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(math_tools LANGUAGES CXX VERSION 2.0)
set(ProjectVersion "v_2_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

message(STATUS "${PROJECT_NAME} version ${PROJECT_VERSION}")


#   Find required packages
find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(benchmark QUIET)
find_package(GTest QUIET)


#   External utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})


#   Give different names for debug and Release
set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")




#   Process the content of src, ie, the actual code
add_subdirectory(src)

# testing
if(GTest_FOUND)
  message(STATUS "Building tests")
  enable_testing()

  add_subdirectory(tests)

endif(GTest_FOUND)


# benchmark
if(benchmark_FOUND)
  message(STATUS "Building benchmark")
  enable_testing()

  add_subdirectory(benchmarks)

endif(benchmark_FOUND)

#   Install the library using the default routine
include(cmake/installation_module.cmake)


#include(cmake/documentation_module.cmake)
