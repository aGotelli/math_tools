#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <iostream>


using namespace ::LieAlgebra;
using namespace Eigen;



Matrix6d Ad_invg(Matrix3d R, Vector3d v)
{
    Matrix6d Ad;
    Ad << R.transpose(), - R.transpose() * skew(v),
          Matrix3d::Zero(), R.transpose();

    return Ad;
}

Matrix6d Adg(Matrix3d R, Vector3d v)
{
    Matrix6d Ad;
    Ad << R               , skew(v) *R,
          Matrix3d::Zero(),     R;

    return Ad;
}

int main(int argc, char *argv[])
{

    Vector3d v = Vector3d::Random();


    Quaterniond q = Quaterniond::UnitRandom();
    Matrix3d R = q.toRotationMatrix();





    std::cout << "Ad_invg : \n" << Ad_invg(R, v) << "\n\n";


    std::cout << "inv_Adg : \n" << Adg(R, v).inverse() << "\n\n";


    auto R_inv = R.transpose();
    auto v_inv = -R.transpose() * v;


    std::cout << "Ad_ginv : \n" << Adg(R_inv, v_inv) << "\n\n";
    std::cout << "Ad_ginv : \n" << Adg(R_inv, v_inv).inverse() << "\n\n";

    return 0;
}
