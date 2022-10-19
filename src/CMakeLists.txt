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
target_compile_options(${PROJECT_NAME}
    PRIVATE
        #-Werror
        -Wall
        -Wextra
)



install(TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

