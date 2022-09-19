

set(SUBPROJECT_NAME "LieAlgebra")


add_library(${SUBPROJECT_NAME}
    SHARED
        lie_algebra_utilities.cpp
        include/lie_algebra_utilities.hpp
)
add_library(${PROJECT_NAME}::${SUBPROJECT_NAME} ALIAS ${SUBPROJECT_NAME})
target_link_libraries(${SUBPROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
)
target_compile_options(${SUBPROJECT_NAME}
    PRIVATE
        -Wall
        -Wextra
)



install(
    TARGETS
        ${SUBPROJECT_NAME}
    DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/${SUBPROJECT_NAME}
)
install(
    DIRECTORY
        include
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${SUBPROJECT_NAME}
)

