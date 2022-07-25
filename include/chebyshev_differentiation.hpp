/**
 * \file chebyshev_differentiation.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This file contains the functions related to Chebyshev points and differentiation matrix
 * \date 06/07/2022
 * 
 * \copyright Copyright (c) 2022
 * 
 */


#ifndef CHEBYSHEV_DIFFERENTIATION_HPP
#define CHEBYSHEV_DIFFERENTIATION_HPP


#include <eigen3/Eigen/Dense>
#include <iostream>

namespace Chebyshev {



enum class INTEGRATION_DIRECTION {
    FORWARD,
    BACKWARD
};



/*!
 * \brief ComputeChebyshevPoints Computes the Chebyshev points in the given interval
 * \tparam t_N The number of Chebyshev points.
 * \tparam t_L The length of the interval. Default 1 for the interval [0, 1]
 * \return An std::vector containing the Chebyshev points
 */
std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_chebyshev_nodes,
                                                    const unsigned int t_L=1);


/*!
 * \brief getDn Computes the Chebyshev differentiation matrix
 * \tparam t_N The number of Chebyshev points.
 * \return The Chebyshev differentiation matrix
 */
Eigen::MatrixXd getDn(const unsigned int t_number_of_chebyshev_nodes);


/*!
 * \brief getD_NN computes the block-diagonal matrix corresponding to the influence of the unknown values onto themselfs
 * \param t_number_of_chebyshev_nodes the number of Chebyshev point to be used in the numerical integration
 * \param t_state_dimension the dimension of the system that is integrated
 * \return The block-diagonal matrix corresponding to the influence of the unknown values onto themselfs
 */
Eigen::MatrixXd getD_NN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction);

/*!
 * \brief getD_IN computes the block-diagonal matrix corresponding to the influence of the initial condition ontt the unknown values
 * \param t_number_of_chebyshev_nodes
 * \param t_state_dimension
 * \return The block-diagonal matrix corresponding to the influence of the initial condition ontt the unknown values
 */
Eigen::MatrixXd getD_IN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction);


std::vector<unsigned int> defineIntegrationPoints(unsigned int t_number_of_chebyshev_nodes,
                                                  INTEGRATION_DIRECTION t_integration_direction);


}   //  namespace Chebyshev

#endif // CHEBYSHEV_DIFFERENTIATION_HPP
