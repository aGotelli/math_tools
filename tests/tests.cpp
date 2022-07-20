

#include "lie_algebra_utilities.hpp"

#include "chebyshev_differentiation.hpp"

int main(int argc, char *argv[])
{
    std::cout << getDn(5) << std::endl;
    std::cout << "Forward :\n   D_NN :\n" << getD_NN(5, 1, INTEGRATION_DIRECTION::FORWARD) << "\n   D_IN :\n" << getD_IN(5,1, INTEGRATION_DIRECTION::FORWARD) << std::endl;
    std::cout << "Backward :\n   D_NN :\n" << getD_NN(5, 1, INTEGRATION_DIRECTION::BACKWARD) << "\n   D_IN :\n" << getD_IN(5,1, INTEGRATION_DIRECTION::BACKWARD) << std::endl;
    return 0;
}
