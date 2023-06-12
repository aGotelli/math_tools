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


/// \namespace Chebyshev contains all the function related to the Chebyshev discretisation theory
namespace Chebyshev {

/// \brief Defines a default value for the Chebyshev to be used in the depending packages
constexpr unsigned int default_number_of_Chebyshev_points = 17;


/*!
 * \brief The INTEGRATION_DIRECTION enum defines the two possible direction of integration: FORWARD and BACKWARD
 */
enum class INTEGRATION_DIRECTION {
    FORWARD,
    BACKWARD
};



/*!
 * \brief ComputeChebyshevPoints Computes the Chebyshev points in the given interval
 * \param t_number_of_chebyshev_points The number of Chebyshev points.
 * \param t_L The length of the interval. Default 1 for the interval [0, 1]
 * \return An std::vector containing the Chebyshev points
 */
std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_chebyshev_points,
                                           const double &t_L=1.0f);


/*!
 * \brief ComputeChebyshevPoints Computes the Chebyshev points in the given interval
 * \param t_number_of_chebyshev_points The number of Chebyshev points.
 * \param t_lower_bound the lower bound of the points domain
 * \param t_upper_bound the upper bound of the points domain
 * \return An std::vector containing the Chebyshev points
 */
std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_chebyshev_points,
                                           const double t_lower_bound,
                                           const double t_upper_bound);


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
 * \param t_integration_direction the direction of integration, specified by the INTEGRATION_DIRECTION enum
 * \return The block-diagonal matrix corresponding to the influence of the unknown values onto themselfs
 */
Eigen::MatrixXd getD_NN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction);

/*!
 * \brief getD_IN computes the block-diagonal matrix corresponding to the influence of the initial condition ontt the unknown values
 * \param t_number_of_chebyshev_nodes
 * \param t_state_dimension
 * \param t_integration_direction the direction of integration, specified by the INTEGRATION_DIRECTION enum
 * \return The block-diagonal matrix corresponding to the influence of the initial condition ontt the unknown values
 */
Eigen::MatrixXd getD_IN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction);


/*!
 * \brief defineIntegrationPoints gives the sequence of integration points coherent with the direction
 * \param t_number_of_chebyshev_nodes is the number of Chebyshev nodes
 * \param t_integration_direction is the direction of integration
 * \return the sequence of integration points as a vector
 *
 * This function gives the sequence of Chebyshev points coherent with the direction.
 * If the intgration if forward then we integrate from the first point up to the second last point in the grid.
 * On the other end, if backward we integrate from the second point up to the last point in the grid.
 */
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




/*!
 * \brief The ChebyshevReconstructor class an handle for the reconstruction of a shape observed with Chebyshev points
 */

class ChebyshevReconstructor {
public:

    /*!
     * \brief ChebyshevReconstructor Initialises the object ChebyshevReconstructor
     */
    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    ChebyshevReconstructor()=default;

    /*!
     * \brief ChebyshevReconstructor Initialises the object ChebyshevReconstructor
     * \param t_number_of_Chebyshev_points the number of Chebyshev point used in the observation
     */
    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points);


    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points,
                           const unsigned int t_number_of_reconstruction_points,
                           const int t_state_dimension)
        : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
          m_state_dimension(t_state_dimension),
          m_number_of_reconstruction_points(t_number_of_reconstruction_points)
    {}

    /*!
     * \brief ChebyshevReconstructor Initialises the object ChebyshevReconstructor
     * \param t_number_of_Chebyshev_points the number of Chebyshev point used in the observation
     * \param t_number_of_reconstruction_points the number of points to use in the reconstruction
     */
    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points,
                          const unsigned int t_number_of_reconstruction_points);

    /*!
     * \brief reconstructRodShape reconstruct the shape of the rod from the observed centerline positions at the Chebyshev points
     * \param t_centerline_points the observed centerline positions at the Chebyshev points
     * \return the shape of the rod from the observed centerline positions at the Chebyshev points
     */
    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    Eigen::MatrixXd reconstructRodShape(const Eigen::MatrixXd &t_centerline_points)const;


    /*!
     * \brief getNumberOfReconstructionPoints gives the number of reconstruction points used to reconstruct the shape of the rod
     * \return the number of reconstruction points used to reconstruct the shape of the rod
     */
    [[deprecated("this class is outdated, you should use the Chebyshev interpolator that is correct and unsensitive to the state dymension")]]
    inline unsigned int getNumberOfReconstructionPoints()const{return m_number_of_reconstruction_points;}

    private:




    //  Number of Chebyshev points
    const unsigned int m_number_of_Chebyshev_points { ::Chebyshev::default_number_of_Chebyshev_points };




    const unsigned int m_state_dimension { 3 };

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




    //  Cosine transform of the fucntion (marked mutable, here we want just to allocate the memory)
    mutable Eigen::MatrixXd m_CN { Eigen::MatrixXd(m_number_of_Chebyshev_points, m_state_dimension) };


    //  Matrix of reconstructed points (marked mutable, here we want just to allocate the memory)
    mutable Eigen::MatrixXd m_reconstructed_points { Eigen::MatrixXd(m_number_of_reconstruction_points, m_state_dimension) };


};

/*!
 * \brief The ChebyshevInterpolator class is a functor to interpolate a given function observed on the Chebyshev grid
 *
 * It uses the Chebyshev interpolation for vector spaces. Thus it is correct for vectors that does not lie in a group.
 *
 * The Chebyshev points can be defined in an albitrary domain.
 */
struct ChebyshevInterpolator {

    /*!
     * \brief ChebyshevInterpolator construct the objects for a given number of Chebyshev points using equispaced interpolation points
     * \param t_number_of_Chebyshev_points the number of Chebyshev points on which the function is observed
     * \param t_number_of_interpolation_points the number of equispaced points in which the function will be interpolated
     * \param t_lower_bound the lower bound of the function domain
     * \param t_upper_bound the upper bound of the function domain
     */
    ChebyshevInterpolator(const unsigned int t_number_of_Chebyshev_points,
                          const unsigned int t_number_of_interpolation_points=101,
                          const double t_lower_bound=0.0,
                          const double t_upper_bound=1.0);


    /*!
     * \brief ChebyshevInterpolator construct the object for a given number of Chebyshev using an albitrary set of interpolation points
     * \param t_Chebyshev_points the set of Chebyshev points on which the function is observed
     * \param t_number_of_interpolation_points the number of equispaced points in which the function will be interpolated
     * \param t_lower_bound the lower bound of the function domain
     * \param t_upper_bound the upper bound of the function domain
     *
     */
    ChebyshevInterpolator(const std::vector<double> &t_Chebyshev_points,
                          const unsigned int t_number_of_interpolation_points=101,
                          const double t_lower_bound=0.0,
                          const double t_upper_bound=1.0);


    /*!
     * \brief ChebyshevInterpolator construct the object for a given number of Chebyshev using an albitrary set of interpolation points
     * \param t_Chebyshev_points the set of Chebyshev points on which the function is observed
     * \param t_interpolation_points the set of interpolation points in the normalised domain [0, 1]
     *
     */
    ChebyshevInterpolator(const unsigned int t_number_of_Chebyshev_points,
                          const Eigen::VectorXd &t_interpolation_points);


    /*!
     * \brief ChebyshevInterpolator construct the object for a given number of Chebyshev using an albitrary set of interpolation points
     * \param t_Chebyshev_points the set of Chebyshev points on which the function is observed
     * \param t_interpolation_points the set of interpolation points in the normalised domain [0, 1]
     *
     */
    ChebyshevInterpolator(const unsigned int t_number_of_Chebyshev_points,
                          const std::vector<double> &t_interpolation_points);


    /*!
     * \brief ChebyshevInterpolator construct the object for a given number of Chebyshev using an albitrary set of interpolation points
     * \param t_number_of_Chebyshev_points the number of Chebyshev points on which the function is observed
     * \param t_interpolation_points the set of interpolation points in the normalised domain [0, 1]
     *
     */
    ChebyshevInterpolator(const std::vector<double> &t_Chebyshev_points,
                          const Eigen::VectorXd &t_interpolation_points);


    /*!
     * \brief operator () interpolates the given function on the defines grid of interpolation points
     * \param t_function_to_interpolate is the function to interpolate of an albitrary state dimension
     *
     *
     * \return a stack containing all the interpolated points
     *
     * The function to interpolate is a stack of column vector being the observed states. That is, the state of the system
     * is observed at every Chebyshev pont. Starting from the begin of the integration domain, the observed states are stacked
     * one after the other. As a result, being the observed state y, then \f$ y(X=0) \f$ will be the first column of the
     * matrix while \f$ y(X=1)\f$ will be the last column of the matrix.
     *
     * As a result, the matrix has as many row as the dimension of the observed state and as many columns as the
     * number of Chebyshev points on which it is observed.
     */
    Eigen::MatrixXd operator()(const Eigen::MatrixXd &t_function_to_interpolate)const;


    /*!
     * \brief getNumberOfInterpolationPoints returns the number of points at which we interpolate (the number of interpolation points)
     * \return the number of interpolation points
     */
    unsigned int getNumberOfInterpolationPoints()const;


private:

    //  The matrix used to interpolate the function on its domain
    Eigen::MatrixXd m_interpolation_matrix;

};




}   //  namespace Chebyshev

#endif // CHEBYSHEV_DIFFERENTIATION_HPP
