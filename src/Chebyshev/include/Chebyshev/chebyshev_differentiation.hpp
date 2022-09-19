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


#include <Eigen/Dense>
#include <iostream>

namespace Chebyshev {

constexpr unsigned int default_number_of_Chebyshev_points = 17;

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
 * \brief getDn_NN gives the submatrix of Dn for the spectral numerical integration
 * \param t_number_of_chebyshev_nodes is the number of Chebyshev points used in the integration
 * \param t_integration_direction is the direction of integration
 * \return The Dn_NN matrix giving the influeces of the unknown points onto themself
 */
Eigen::MatrixXd getDn_NN(const unsigned int t_number_of_chebyshev_nodes,
                        const INTEGRATION_DIRECTION &t_integration_direction);


/*!
 * \brief getDn_IN gives the submatrix of Dn for the spectral numerical integration
 * \param t_number_of_chebyshev_nodes is the number of Chebyshev points used in the integration
 * \param t_integration_direction is the direction of integration
 * \return The Dn_IN matrix giving the influeces of the initial conditions onto the unknown points
 */
Eigen::MatrixXd getDn_IN(const unsigned int t_number_of_chebyshev_nodes,
                        const INTEGRATION_DIRECTION &t_integration_direction);



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



/*!
 * \brief getintialConditionAddress defines the addess of the intial condition for the integration
 * \param t_number_of_chebyshev_nodes defines how many Chebyshev points are used in the discretization
 * \param t_integration_direction defines the direction of integration
 * \return The address of the initial conditions
 */
unsigned int getintialConditionAddress(unsigned int t_number_of_chebyshev_nodes,
                                       INTEGRATION_DIRECTION t_integration_direction);





class ChebyshevReconstructor {
public:
    ChebyshevReconstructor()=default;

    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points);


    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points,
                         const unsigned int t_number_of_reconstruction_points);


    Eigen::MatrixXd ReconstructRodShape(const Eigen::Matrix<double, 3, Eigen::Dynamic> &t_centerline_points);


    inline unsigned int getNumberOfReconstructionPoints()const{return m_number_of_Chebyshev_points;}

    private:

    //  Number of Chebyshev points
    const unsigned int m_number_of_Chebyshev_points { ::Chebyshev::default_number_of_Chebyshev_points };

    //  Number of points describing the rod centerline
    const unsigned int m_number_of_reconstruction_points { 101 };

    //  Discretize the rod in several points of observations
    const double m_step { 1.0/(m_number_of_reconstruction_points -1) };

    //  Direct Discrete Cosine Transform
    const Eigen::MatrixXd m_DDCT { [&](){
            //  Initialize the matrix
            Eigen::MatrixXd DDCT(m_number_of_Chebyshev_points, m_number_of_Chebyshev_points);

            const auto k = [](const unsigned int h){ return h==0 ? 1/sqrt(2) : pow((-1), h);};

            for(unsigned int h=0; h<=m_number_of_Chebyshev_points-1;h++) {
                for(unsigned int j=0; j<=m_number_of_Chebyshev_points-1;j++){
                    DDCT(h, j) = (2.0/m_number_of_Chebyshev_points)
                                * k(h)*cos( (2.0*j + 1.0)*h*M_PI/(2.0*m_number_of_Chebyshev_points) );
                }
            }

            //  Return the matrix
            return DDCT;

          }()};




    //  Cosine transform of the fucntion
    Eigen::MatrixXd m_CN { Eigen::MatrixXd(m_number_of_Chebyshev_points, 3) };


    //  Matrix of reconstructed points
    Eigen::MatrixXd m_reconstructed_points { Eigen::MatrixXd(m_number_of_reconstruction_points, 3) };


};

}   //  namespace Chebyshev

#endif // CHEBYSHEV_DIFFERENTIATION_HPP
