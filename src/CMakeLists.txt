include_directories(include)



add_library(${PROJECT_NAME}
    SHARED
        Chebyshev/chebyshev_differentiation.cpp
        include/math_tools/Chebyshev/chebyshev_differentiation.hpp
        LieAlgebra/lie_algebra_utilities.cpp
        include/math_tools/LieAlgebra/lie_algebra_utilities.hpp
)
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        Eigen3::Eigen
)
install(TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    EXPORT ${PROJECT_NAME}Export
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

