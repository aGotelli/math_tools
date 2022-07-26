cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

project(math_tools LANGUAGES CXX VERSION 2.0)
set(ProjectVersion "v_2_0")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RELEASE_POSTFIX "")
set(CMAKE_DEBUG_POSTFIX "-debug")

message(STATUS "${PROJECT_NAME} version ${PROJECT_VERSION}")

# external utilities to configure the package
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})





find_package(Eigen3 3.4 NO_MODULE REQUIRED)

find_package(benchmark REQUIRED)


# process the content of src, ie, the actual code
add_subdirectory(src)






# generates URLConfig.cmake file to use our package in other projects
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
