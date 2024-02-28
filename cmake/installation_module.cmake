#   Generates ${PROJECT_NAME}Config.cmake file to use our package in other projects
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/Config.cmake.in # input template
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake # output config file
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake # where to put the config file during install
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR # paths to be used
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

#   Generates a config file to ensure that URL's version is checked when importing it
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION
    ${PROJECT_VERSION}
  COMPATIBILITY
    SameMajorVersion
)



#   When running make install, config files should be copied as well
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

#   Add the possibility tu run 'make uninstall' to remove files added via 'make install'
configure_file(
  cmake/cmake_uninstall.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
  IMMEDIATE @ONLY
)
add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
