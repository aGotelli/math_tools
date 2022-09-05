
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <gtest/gtest.h>

using namespace ::LieAlgebra;


TEST(test_SE3Pose, translation)
{
    SE3Pose pose1(Eigen::Vector3d(1, 2, 1));
    SE3Pose pose2(Eigen::Vector3d(1, 0.5, -1.8));

    const SE3Pose pose_result = pose1 * pose2;

    EXPECT_NEAR(pose_result.m_position.x(), 2.0, 1e-14);
    EXPECT_NEAR(pose_result.m_position.y(), 2.5, 1e-14);
    EXPECT_NEAR(pose_result.m_position.z(), -0.8, 1e-14);
}


TEST(test_SE3Pose, null_rotation)
{
    double roll = 0.0;
    double pitch = 0.0;
    double yaw = 0.0;

    Eigen::Quaterniond q1;
    q1 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());


    roll = 0.0;
    pitch = 0.0;
    yaw = 0.0;


    Eigen::Quaterniond q2;
    q2 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());


    SE3Pose pose1(q1, Eigen::Vector3d::Zero());
    SE3Pose pose2(q2, Eigen::Vector3d::Zero());

    const SE3Pose pose_result = pose1 * pose2;


    EXPECT_EQ(pose_result.m_quaternion.w(), 1.0);
    EXPECT_EQ(pose_result.m_quaternion.x(), 0.0);
    EXPECT_EQ(pose_result.m_quaternion.y(), 0.0);
    EXPECT_EQ(pose_result.m_quaternion.z(), 0.0);
}


TEST(test_SE3Pose, rotation_along_x)
{
    double roll = M_PI/2;
    double pitch = 0.0;
    double yaw = 0.0;

    Eigen::Quaterniond q1;
    q1 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());

    const auto rotated_quaternion = rotateAlongAxis(roll, Eigen::Vector3d::UnitX());


    roll = 0.0;
    pitch = 0.0;
    yaw = 0.0;


    Eigen::Quaterniond q2;
    q2 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());


    SE3Pose pose1(q1, Eigen::Vector3d::Zero());
    SE3Pose pose2(q2, Eigen::Vector3d::Zero());

    const SE3Pose pose_result = pose1 * pose2;




    EXPECT_EQ(pose_result.m_quaternion.w(), rotated_quaternion.w());
    EXPECT_EQ(pose_result.m_quaternion.x(), rotated_quaternion.x());
    EXPECT_EQ(pose_result.m_quaternion.y(), rotated_quaternion.y());
    EXPECT_EQ(pose_result.m_quaternion.z(), rotated_quaternion.z());
}


TEST(test_SE3Pose, rotation_along_xy)
{
    double roll = M_PI/2;
    double pitch = 0.0;
    double yaw = 0.0;

    Eigen::Quaterniond q1;
    q1 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());

    const auto rotated_quaternion1 = rotateAlongAxis(roll, Eigen::Vector3d::UnitX());


    roll = 0.0;
    pitch = -M_PI/4;
    yaw = 0.0;


    Eigen::Quaterniond q2;
    q2 = Eigen::AngleAxisd(roll  , Eigen::Vector3d::UnitX())
        * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(yaw  , Eigen::Vector3d::UnitZ());

    const auto rotated_quaternion2 = rotateAlongAxis(pitch, Eigen::Vector3d::UnitY());

    const auto result_quaternion = rotated_quaternion1 * rotated_quaternion2;


    SE3Pose pose1(q1, Eigen::Vector3d::Zero());
    SE3Pose pose2(q2, Eigen::Vector3d::Zero());

    const SE3Pose pose_result = pose1 * pose2;




    EXPECT_EQ(pose_result.m_quaternion.w(), result_quaternion.w());
    EXPECT_EQ(pose_result.m_quaternion.x(), result_quaternion.x());
    EXPECT_EQ(pose_result.m_quaternion.y(), result_quaternion.y());
    EXPECT_EQ(pose_result.m_quaternion.z(), result_quaternion.z());
}


TEST(test_SE3Pose, rotation_along_xyz)
{
    const Eigen::Quaterniond q1 = Eigen::Quaterniond::UnitRandom();



    const Eigen::Quaterniond q2 = Eigen::Quaterniond::UnitRandom();


    SE3Pose pose1(q1, Eigen::Vector3d::Zero());
    SE3Pose pose2(q2, Eigen::Vector3d::Zero());

    const SE3Pose pose_result = pose1 * pose2;

    const Eigen::Quaterniond q_result = q1 * q2;



    EXPECT_EQ(pose_result.m_quaternion.w(), q_result.w());
    EXPECT_EQ(pose_result.m_quaternion.x(), q_result.x());
    EXPECT_EQ(pose_result.m_quaternion.y(), q_result.y());
    EXPECT_EQ(pose_result.m_quaternion.z(), q_result.z());
}



TEST(test_SE3Pose, full)
{
    const Eigen::Vector3d r1 = Eigen::Vector3d::Random();
    const Eigen::Quaterniond q1 = Eigen::Quaterniond::UnitRandom();


    const Eigen::Vector3d r2 = Eigen::Vector3d::Random();
    const Eigen::Quaterniond q2 = Eigen::Quaterniond::UnitRandom();


    SE3Pose pose1(q1, r1);
    SE3Pose pose2(q2, r2);

    const SE3Pose pose_result = pose1 * pose2;

    const Eigen::Quaterniond q_result = q1 * q2;
    const Eigen::Vector3d r_result = q1.toRotationMatrix()*r2 + r1;



    EXPECT_EQ(pose_result.m_quaternion.w(), q_result.w());
    EXPECT_EQ(pose_result.m_quaternion.x(), q_result.x());
    EXPECT_EQ(pose_result.m_quaternion.y(), q_result.y());
    EXPECT_EQ(pose_result.m_quaternion.z(), q_result.z());

    EXPECT_EQ(pose_result.m_position.x(), r_result.x());
    EXPECT_EQ(pose_result.m_position.y(), r_result.y());
    EXPECT_EQ(pose_result.m_position.z(), r_result.z());
}


TEST(test_SE3Pose, inverse)
{
    //  Initialise the transformation
    const Eigen::Vector3d r = Eigen::Vector3d::Ones();
    const Eigen::Quaterniond q = rotateAlongAxis(M_PI/3, Eigen::Vector3d::UnitZ());
    SE3Pose pose(q, r);

    //  Compute the inverse internally
    const auto inverse_pose = pose.inverse();

    //  Compute the inverse analytically
    const Eigen::Matrix3d analytic_R_inv = q.toRotationMatrix().transpose();
    const auto analytic_r_inv = -analytic_R_inv*r;

    //  Compare results in rotation
    const auto so3_error = differenceInSO3(inverse_pose.getRotationMatrix(),analytic_R_inv);
    EXPECT_NEAR(so3_error.x(), 0.0, 1e-14);
    EXPECT_NEAR(so3_error.y(), 0.0, 1e-14);
    EXPECT_NEAR(so3_error.z(), 0.0, 1e-14);


    //  Compare result in position
    EXPECT_EQ(inverse_pose.m_position.x(), analytic_r_inv.x());
    EXPECT_EQ(inverse_pose.m_position.y(), analytic_r_inv.y());
    EXPECT_EQ(inverse_pose.m_position.z(), analytic_r_inv.z());


}



int main(int argc, char *argv[])
{

    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
