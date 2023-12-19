




include_directories(../src/include)
include_directories(include)






add_executable(benchmark_error_SO3 benchmark_error_SO3.cpp)
target_link_libraries(benchmark_error_SO3
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        benchmark::benchmark
)


add_executable(benchmark_T_SO3 benchmark_T_SO3.cpp)
target_link_libraries(benchmark_T_SO3
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        benchmark::benchmark
)



