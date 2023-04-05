include_directories(include)



add_library(${PROJECT_NAME}
    SHARED
        Chebyshev/chebyshev_differentiation.cpp
        include/math_tools/Chebyshev/chebyshev_differentiation.hpp
        LieAlgebra/lie_algebra_utilities.cpp
        include/math_tools/LieAlgebra/lie_algebra_utilities.hpp
)
target_link_libraries(${PROJECT_NAME}
    PUBLIC
        Eigen3::Eigen
)
target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wall
        -Wextra
)



install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Exports
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT ${PROJECT_NAME}Exports
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION lib/cmake/${PROJECT_NAME}
)

