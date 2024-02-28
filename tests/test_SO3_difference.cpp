#include <gtest/gtest.h>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"



const static auto alpha_1 = M_PI/6;
const static auto alpha_2 = -M_PI/6;

const static Eigen::Matrix3d Ra = ::LieAlgebra::getRx(alpha_1);
const static Eigen::Matrix3d Rb = ::LieAlgebra::getRx(alpha_2);


TEST(test_SO3_difference, rotation_along_x)
{
    const auto so3_difference = ::LieAlgebra::differenceInSO3(Ra, Rb);

    EXPECT_NEAR(so3_difference.x(), 0.0, 1e-14);
    EXPECT_NEAR(so3_difference.y(), 0.0, 1e-14);
    EXPECT_NEAR(so3_difference.z(), 0.0, 1e-14);
}



TEST(test_SO3_difference, log_rotation_along_x)
{
    const Eigen::Vector3d log_SO3 = ::LieAlgebra::logSO3(Ra, Rb);


    EXPECT_NEAR(log_SO3.x(), 0.0, 1e-14);
    EXPECT_NEAR(log_SO3.y(), 0.0, 1e-14);
    EXPECT_NEAR(log_SO3.z(), 0.0, 1e-14);
}





int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
