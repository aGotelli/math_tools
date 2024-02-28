#include <iostream>
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"











int main(int argc, char *argv[])
{
    ::LieAlgebra::SE3Pose pose;
    pose.m_position = Eigen::Vector3d(1, 2, 3);
    pose.m_quaternion = Eigen::Quaterniond::UnitRandom();

    std::cout << "r : \n\n" << pose.m_position << "\n";
    std::cout << "R : \n\n" << pose.m_quaternion.toRotationMatrix() << "\n";

    std::cout << "\n\n";

    std::cout << "g : \n\n" << pose.getSE3PoseAsMatrix() << "\n";

    return 0;
}
