include_directories(include)


add_subdirectory(Chebyshev)

add_subdirectory(LieAlgebra)



add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}
  INTERFACE
    ${PROJECT_NAME}::Chebyshev
    ${PROJECT_NAME}::LieAlgebra
)


