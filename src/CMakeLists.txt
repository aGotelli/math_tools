include_directories(include)



add_library(Chebyshev Chebyshev/chebyshev_differentiation.cpp
                      include/math_tools/Chebyshev/chebyshev_differentiation.hpp)
target_link_libraries(Chebyshev PRIVATE Eigen3::Eigen)
install(TARGETS Chebyshev DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(LieAlgebra LieAlgebra/lie_algebra_utilities.cpp
                      include/math_tools/LieAlgebra/lie_algebra_utilities.hpp)
target_link_libraries(LieAlgebra PRIVATE Eigen3::Eigen)
install(TARGETS LieAlgebra DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
