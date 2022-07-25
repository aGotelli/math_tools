cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(math_tools LANGUAGES CXX VERSION 1.3)
set(ProjectVersion "v_1_3")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# external utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})



find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(benchmark REQUIRED)


include_directories(include)



add_library(${PROJECT_NAME} SHARED    include/chebyshev_differentiation.hpp
                                      src/chebyshev_differentiation.cpp
                                      include/lie_algebra_utilities.hpp
                                      src/lie_algebra_utilities.cpp)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)



add_executable(tests tests/tests.cpp)
target_link_libraries(tests math_tools)




install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})


# generates Config.cmake file to use our package in other projects
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/Config.cmake.in # input template
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake # output config file
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake # where to put the config file during install
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR # paths to be used
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# generates a config file to ensure that URL's version is checked when importing it
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# when running make install, config files should be copied as well
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# add the possibility tu run 'make uninstall' to remove files added via 'make install'
configure_file(
  cmake/cmake_uninstall.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  IMMEDIATE @ONLY
)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
