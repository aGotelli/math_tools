

set(SUBPROJECT_NAME "Chebyshev")

include_directories(include)


add_library(${SUBPROJECT_NAME}
    SHARED
        src/chebyshev_differentiation.cpp
        include/${SUBPROJECT_NAME}/chebyshev_differentiation.hpp
)
target_link_libraries(${SUBPROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
)
target_compile_options(${SUBPROJECT_NAME}
    PRIVATE
        -Wall
        -Wextra
)

add_library(${PROJECT_NAME}::${SUBPROJECT_NAME} ALIAS ${SUBPROJECT_NAME})



install(
    TARGETS
        ${SUBPROJECT_NAME}
    DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/${SUBPROJECT_NAME}
)
install(
    DIRECTORY
        include/${SUBPROJECT_NAME}/
    DESTINATION
        ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${SUBPROJECT_NAME}
)

