/**
 * \file chebyshev_differentiation.cpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the base class for integrating a linear ode with the spectral numerical integration
 * \version 0.1
 * \date 06-07-2022
 *
 * \copyright Copyright (c) 2022
 *
 */

#include "chebyshev_differentiation.hpp"


std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_chebyshev_nodes,
                                                    const unsigned int t_L)
{
    std::vector<double> x(t_number_of_chebyshev_nodes);

    unsigned int j = 0;
    std::generate(x.begin(), x.end(), [&](){
        return (static_cast<double>(t_L)/2)*(1 +cos( M_PI * static_cast<double>(j++) / static_cast<double>(t_number_of_chebyshev_nodes-1) ));
    });

    return x;
}




std::vector<double> GetCoefficients_c(const unsigned int t_number_of_chebyshev_nodes)
{
    std::vector<double> c(t_number_of_chebyshev_nodes);

    unsigned int i = 0;
    std::generate(c.begin(), c.end(), [&](){
        //  gain is 2 in the edges and 1 elsewhere
        const unsigned int gain = (i==0 or i==t_number_of_chebyshev_nodes-1) ? 2 : 1;

        //  Follows the formula
        return pow(-1, i++)*gain;
    });

    return c;
}



Eigen::MatrixXd getDn(const unsigned int t_number_of_chebyshev_nodes)
{

    //  Define the Chebyshev points on the unit circle
    const auto x = ComputeChebyshevPoints(t_number_of_chebyshev_nodes);


    //  Create a matrix every row filled with a point value
    Eigen::MatrixXd X(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);
    for(unsigned int i=0; i<X.rows(); i++)
        X(i, Eigen::all) = Eigen::RowVectorXd::Constant(1, X.cols(), x[i]);


    //  Now compute the array containing the coefficients used in the definition of Dn
    const auto c = GetCoefficients_c(t_number_of_chebyshev_nodes);


    //  Create the appropriate matrix of coefficients
    Eigen::MatrixXd C(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);
    for(unsigned int i=0; i<t_number_of_chebyshev_nodes;i++) {
        for(unsigned int j=0; j<t_number_of_chebyshev_nodes;j++) {
            C(i,j) = c[i]/c[j];
        }
    }

    //  Definition of the temporary matrix Y
    const auto dX = X - X.transpose() + Eigen::MatrixXd::Identity(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);

    //  Declare the differentiation matrix
    Eigen::MatrixXd  Dn(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);


    //  Obtain off diagonal element for the differentiation matrix
    for(unsigned int i=0; i<t_number_of_chebyshev_nodes;i++) {
        for(unsigned int j=0; j<t_number_of_chebyshev_nodes;j++) {
            Dn(i,j) = C(i, j) / dX(i, j);
        }
    }


    //  Remove row sum from the diagonal of Dn
    Dn.diagonal() -= Dn.rowwise().sum();

    //  Finally return the matrix
    return Dn;
}

