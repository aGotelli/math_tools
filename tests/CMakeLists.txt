




include_directories(../src/include)
include_directories(include)




add_executable(test_log test_log.cpp)
target_link_libraries(test_log
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
)

add_executable(test_SE3Pose test_SE3Pose.cpp)
target_link_libraries(test_SE3Pose
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        GTest::GTest
)
add_test(test_SE3Pose  test_SE3Pose)

add_executable(test_SO3_difference test_SO3_difference.cpp)
target_link_libraries(test_SO3_difference
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        GTest::GTest
)
add_test(test_SO3_difference  test_SO3_difference)


add_executable(test_Adg test_Adg.cpp)
target_link_libraries(test_Adg
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
)

add_executable(test_T test_T.cpp)
target_link_libraries(test_T
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
)

add_executable(test_log_SO3 test_log_SO3.cpp)
target_link_libraries(test_log_SO3
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
)


find_package(utilities REQUIRED)

add_executable(test_chebyshev_reconstructor test_chebyshev_reconstructor.cpp)
target_link_libraries(test_chebyshev_reconstructor
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        utilities
)


add_executable(test_expSO3 test_expSO3.cpp)
target_link_libraries(test_expSO3
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
)
