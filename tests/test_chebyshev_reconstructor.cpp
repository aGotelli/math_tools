#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"


static const double center = 0.1;
static const double height = 1;
static const double base = 0.5;

double gaussian(double x)
{
    const double y = x - base;
    const double exp_mem = pow(y, 2)/(2 * pow(center, 2));
    const double gaussian_exp = exp( -exp_mem );

    return height * gaussian_exp;
}


#include "utilities/Eigen/eigen_io.hpp"

int main(int argc, char *argv[])
{
    //  Define the dimension of the state
    const unsigned int state_dimension = 7;
    //  Define number of Chebyshev points
    const unsigned int Nc = 21;

    //  translate into the same variable of book
    const unsigned int N  = Nc - 1;

    //  Define the number of interpolation points
    const unsigned int number_of_interpolation_points = N*3;


    //  Define the domain
    const double min = 0;
    const double max = 1.5;

    //  Define a function to interpolate
    auto f = [&](const double x) {
        Eigen::VectorXd state(state_dimension);
        state << x,
                 x*x,
                 x*x *(1 + sin(x)),
                 exp(x),
                 sqrt(x),
                 atan(x),
                 gaussian(x);
//        state.setConstant( x*x + sin(x) );

        return state;};

    auto Chebyshev_points = ::Chebyshev::ComputeChebyshevPoints(Nc, min, max);

//    std::cout << "Chebyshev points : \n";
//    for(const auto point : Chebyshev_points)
//        std::cout << point << "\n";
//    std::cout << "\n\n\n";



    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(state_dimension, Nc);
    for(unsigned int i=0; i<Nc; i++){
        auto x = Chebyshev_points[i];
        F.col(i) = f( x );
    }

//    std::cout << "F: \n" << F << "\n\n";





    ::Chebyshev::ChebyshevInterpolator reconstructor(Chebyshev_points, number_of_interpolation_points);

    auto interpolated_function = reconstructor(F);










    //  Define the actual interpolation points
    Eigen::VectorXd interpolation_points = Eigen::VectorXd::LinSpaced(number_of_interpolation_points, min, max);

    //  Pre-allocate memory for the matrix
    Eigen::MatrixXd f_analytic = Eigen::MatrixXd::Zero(state_dimension, number_of_interpolation_points);


    //  Compute the approximation
    for(unsigned int i=0; i<number_of_interpolation_points; i++){

            auto x = interpolation_points(i);
            f_analytic(Eigen::all,i) = f(x);
    }

//    std::cout << "f_analytic: \n" << f_analytic << "\n\n";

//    std::cout << "interpolated_function: \n" << interpolated_function << "\n\n";

    Eigen::MatrixXd error = interpolated_function - f_analytic;



    std::cout << "f = x error norm:               " << error.row(0).norm() << "\n";
    std::cout << "f = x² error norm:              " << error.row(1).norm() << "\n";
    std::cout << "f = x² (1 + sin(x)) error norm: " << error.row(2).norm() << "\n";
    std::cout << "f = e^x error norm:             " << error.row(3).norm() << "\n";
    std::cout << "f = sqrt(x) error norm:         " << error.row(4).norm() << "\n";
    std::cout << "f = atan(x) error norm:         " << error.row(5).norm() << "\n";
    std::cout << "f = gaussian(x) error norm:     " << error.row(6).norm() << "\n";

    std::cout << "\n\n\n" <<
                 "global error norm : " << error.norm() << "\n\n";


    Eigen::VectorXd Chebyshev_points_Eigen(Chebyshev_points.size());
    for(unsigned int i=0; i<Chebyshev_points.size(); i++)
        Chebyshev_points_Eigen[i] = Chebyshev_points[i];


    writeToFile("Chebyshev_points", Chebyshev_points_Eigen, "../../../tests/data");
    writeToFile("interpolation_points", interpolation_points, "../../../tests/data");
    writeToFile("analytic_functions", f_analytic, "../../../tests/data");
    writeToFile("interpolated_function", interpolated_function, "../../../tests/data");
    writeToFile("observed_function", F, "../../../tests/data");


    return 0;
}
