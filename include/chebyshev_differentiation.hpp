/**
 * \file chebyshev_differentiation.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This file contains the functions related to Chebyshev points and differentiation matrix
 * \version 0.1
 * \date 06/07/2022
 * 
 * \copyright Copyright (c) 2022
 * 
 */


#ifndef CHEBYSHEV_DIFFERENTIATION_HPP
#define CHEBYSHEV_DIFFERENTIATION_HPP


#include <eigen3/Eigen/Dense>
#include <iostream>

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



#endif // CHEBYSHEV_DIFFERENTIATION_HPP
