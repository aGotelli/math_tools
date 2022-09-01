




include_directories(../src/include)
include_directories(include)






add_executable(test_SE3Pose test_SE3Pose.cpp)
target_link_libraries(test_SE3Pose
    PRIVATE
        Eigen3::Eigen
        ${PROJECT_NAME}
        GTest::GTest
)



