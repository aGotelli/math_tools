#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <iostream>







int main(int argc, char *argv[])
{


    const Eigen::Matrix3d R1 = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d R2;

    const unsigned int iterations = 1001;
    const double max_angle = 2*M_PI;
    const double step = max_angle / (iterations-1);
    double angle;
    for(unsigned int i=0; i< iterations; i++){
        angle = i*step;

        R2 = ::LieAlgebra::getRz( angle );

        const auto error = ::LieAlgebra::logSO3(R1, R2);

        std::cout << "angle : " << angle << " (" << angle*180/M_PI << "°); error : " << error.z() << " (" << error.z()*180/M_PI << "°)\n" << std::endl;

    }




    return 0;
}
