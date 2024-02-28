#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <iostream>


using namespace ::LieAlgebra;
using namespace Eigen;



int main(int argc, char *argv[])
{

    Quaterniond q = Quaterniond::UnitRandom();

    const auto R = q.toRotationMatrix();


    const auto Theta = ::LieAlgebra::logSO3(R);


    std::cout << "exp Rod : \n" << expRodigues(Theta) << "\n\n";
    std::cout << "exp SO(3) : \n" << expSO3(Theta) << "\n\n";

    return 0;
}
