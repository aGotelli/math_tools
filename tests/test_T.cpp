#include <iostream>

#include <Eigen/Dense>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

int main(int argc, char *argv[])
{
    const Eigen::Quaterniond q = Eigen::Quaterniond::UnitRandom();
    std::cout << "q : \n" << q << "\n\n" << std::endl;

    const Eigen::Matrix3d R = q.toRotationMatrix();
    std::cout << "R : \n" << R << "\n\n" << std::endl;

    const Eigen::Vector3d Theta = Eigen::Vector3d(-0.85911, -2.44664, 0.51629); //::LieAlgebra::logSO3(R);
    std::cout << "Theta : \n" << Theta << "\n\n" << std::endl;

//    const Eigen::Matrix3d T_SO3 = ::LieAlgebra::TSO3(Theta);
//    std::cout << "T_SO3 : \n" << T_SO3 << "\n\n" << std::endl;

//    const Eigen::Matrix3d T_SO3_inv = ::LieAlgebra::TSO3_inverse(Theta);
//    std::cout << "T_SO3_inv : \n" << T_SO3_inv << "\n\n" << std::endl;

//    const Eigen::Matrix3d T_SO3_inv_numeric = ::LieAlgebra::TSO3(Theta).inverse();
//    std::cout << "T_SO3_inv_numeric : \n" << T_SO3_inv_numeric << "\n\n" << std::endl;



    const Eigen::Vector3d Delta_Theta = Eigen::Vector3d(0.59688, 0.823295, -0.604897); //Eigen::Vector3d::Random();
    std::cout << "Delta_Theta : \n" << Delta_Theta << "\n\n" << std::endl;

    const Eigen::Matrix3d Delta_T_SO3 = ::LieAlgebra::DeltaTSO3(Theta, Delta_Theta);
//    std::cout << "Delta_T_SO3 : \n" << Delta_T_SO3 << "\n\n" << std::endl;
    std::cout << "Delta_T_SO3_cpp = [\n" << Delta_T_SO3.row(0) << ";\n"
                                     << Delta_T_SO3.row(1) << ";\n"
                                     << Delta_T_SO3.row(2) << ";\n"


              << "]\n\n" << std::endl;


    return 0;
}
