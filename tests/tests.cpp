

#include "lie_algebra_utilities.hpp"

#include "chebyshev_differentiation.hpp"

int main(int argc, char *argv[])
{
    std::cout << ::Chebyshev::getDn(5) << std::endl;
    std::cout << "Forward :\n   ::Chebyshev::D_NN :\n" << ::Chebyshev::getD_NN(5, 1, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD) << "\n   ::Chebyshev::D_IN :\n" << ::Chebyshev::getD_IN(5,1, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD) << std::endl;
    std::cout << "Backward :\n   ::Chebyshev::D_NN :\n" << ::Chebyshev::getD_NN(5, 1, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD) << "\n   ::Chebyshev::D_IN :\n" << ::Chebyshev::getD_IN(5,1, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD) << std::endl;

    std::cout << "\n\n\n";
    std::cout << "Forward integration points : ";
    for(const auto point : ::Chebyshev::defineIntegrationPoints(11, ::Chebyshev::INTEGRATION_DIRECTION::FORWARD))
        std::cout << point << ", ";
    std::cout << "\n\n\n";

    std::cout << "Backward integration points : ";
    for(const auto point : ::Chebyshev::defineIntegrationPoints(11, ::Chebyshev::INTEGRATION_DIRECTION::BACKWARD))
        std::cout << point << ", ";
    std::cout << "\n\n\n";

    return 0;
}
