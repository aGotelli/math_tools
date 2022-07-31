find_package(Doxygen REQUIRED)

set(DOXYGEN_QUIET YES)
set(DOXYGEN_RECURSIVE YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_MAN NO)
set(DOXYGEN_MARKDOWN_SUPPORT YES)
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
set(DOXYGEN_EXTRACT_PACKAGE YES)
set(DOXYGEN_EXTRACT_STATIC YES)
set(DOXYGEN_SHOW_INCLUDE_FILES YES)
set(DOXYGEN_BINARY_TOC YES)
set(DOXYGEN_TOC_EXPAND YES)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "index.md")

#doxygen_add_docs(benchmark_doxygen
#    docs
#    include
#    src
#    #ALL
#    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#    COMMENT "Building documentation with Doxygen."
#)

#install(
#    DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html/"
#    DESTINATION ${CMAKE_INSTALL_DOCDIR})
