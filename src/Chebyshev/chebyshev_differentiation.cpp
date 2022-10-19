#include "math_tools/Chebyshev/chebyshev_differentiation.hpp"
#include <unsupported/Eigen/KroneckerProduct>
#include <numeric>
#include <vector>

#include <boost/math/special_functions/chebyshev.hpp>
namespace Chebyshev {

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



/*!
 * \brief GetCoefficients_c computes the c coefficients of the Chebyshev matrix as an std::vector<double>
 * \param t_number_of_chebyshev_nodes the number of Chebyshev nodes used in the discretisation
 * \return the c coefficients of the Chebyshev matrix as an std::vector<double>
 */
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




Eigen::MatrixXd getDn_NN(const unsigned int t_number_of_chebyshev_nodes,
                        const INTEGRATION_DIRECTION &t_integration_direction)
{
    //  Get the Chebyshev differentiation matrix
    const Eigen::MatrixXd Dn = getDn(t_number_of_chebyshev_nodes);

    //  Extract the block that define influence of initial condition into the unknown states and invert it
    if(t_integration_direction == INTEGRATION_DIRECTION::FORWARD)
        return Dn.block(0, 0, Dn.rows()-1, Dn.cols()-1);
    else
        return Dn.block(1, 1, Dn.rows()-1, Dn.cols()-1);

}



Eigen::MatrixXd getDn_IN(const unsigned int t_number_of_chebyshev_nodes,
                        const INTEGRATION_DIRECTION &t_integration_direction)
{
    //  Get the Chebyshev differentiation matrix
    const Eigen::MatrixXd Dn = getDn(t_number_of_chebyshev_nodes);


    //  Extract the block that define influence of initial condition into the unknown states
    if(t_integration_direction == INTEGRATION_DIRECTION::FORWARD)
        return Dn.block(0, t_number_of_chebyshev_nodes-1, t_number_of_chebyshev_nodes-1, 1);
    else
        return Dn.block(1, 0, t_number_of_chebyshev_nodes-1, 1);

}




Eigen::MatrixXd getD_NN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction)
{
    //  Get the Chebyshev differentiation matrix
    const Eigen::MatrixXd Dn = getDn(t_number_of_chebyshev_nodes);

    Eigen::MatrixXd Dn_NN;

    if(t_integration_direction == INTEGRATION_DIRECTION::FORWARD) {

        //  Extract the block that define mutual infualces of the unknown states
        Dn_NN = Dn.block(0, 0, t_number_of_chebyshev_nodes-1, t_number_of_chebyshev_nodes-1);

    } else { // INTEGRATION_DIRECTION::BACKWARD

        //  Extract the block that define mutual infualces of the unknown states
        Dn_NN = Dn.block(1, 1, t_number_of_chebyshev_nodes-1, t_number_of_chebyshev_nodes-1);

    }

    //  Make a block diagonal matrix with the Dn_NN matrices
    const Eigen::MatrixXd D_NN  = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(t_state_dimension, t_state_dimension), Dn_NN);

    return D_NN;

}

Eigen::MatrixXd getD_IN(const unsigned int t_number_of_chebyshev_nodes,
                        const unsigned int t_state_dimension,
                        const INTEGRATION_DIRECTION &t_integration_direction)
{

    //  Get the Chebyshev differentiation matrix
    const Eigen::MatrixXd Dn = getDn(t_number_of_chebyshev_nodes);

    Eigen::MatrixXd Dn_IN;

    if(t_integration_direction == INTEGRATION_DIRECTION::FORWARD) {

        //  Extract the block that define influence of initial condition into the unknown states
        Dn_IN = Dn.block(0, t_number_of_chebyshev_nodes-1, t_number_of_chebyshev_nodes-1, 1);

    } else { // INTEGRATION_DIRECTION::BACKWARD

        //  Extract the block that define influence of initial condition into the unknown states
        Dn_IN = Dn.block(1, 0, t_number_of_chebyshev_nodes-1, 1);

    }

    //  Make a block diagonal matrix with the Dn_NN matrices
    const Eigen::MatrixXd D_IN = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(t_state_dimension, t_state_dimension), Dn_IN);

    return D_IN;


}


std::vector<unsigned int> defineIntegrationPoints(unsigned int t_number_of_chebyshev_nodes,
                                                  INTEGRATION_DIRECTION t_integration_direction)
{
    std::vector<unsigned int> integration_points(t_number_of_chebyshev_nodes-1);

    if(t_integration_direction == INTEGRATION_DIRECTION::FORWARD)
        std::iota(integration_points.begin(), integration_points.end(), 0);
    else
        std::iota(integration_points.begin(), integration_points.end(), 1);

    return integration_points;
}


unsigned int getintialConditionAddress(unsigned int t_number_of_chebyshev_nodes,
                                       INTEGRATION_DIRECTION t_integration_direction)
{
     return t_integration_direction == INTEGRATION_DIRECTION::FORWARD ? t_number_of_chebyshev_nodes-1 : 0;
}





ChebyshevReconstructor::ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points)
      : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points)
  {}


ChebyshevReconstructor::ChebyshevReconstructor(const unsigned int t_number_of_Chebyshev_points,
                         const unsigned int t_number_of_reconstruction_points)
      : m_number_of_Chebyshev_points(t_number_of_Chebyshev_points),
        m_number_of_reconstruction_points(t_number_of_reconstruction_points)
  {}


Eigen::MatrixXd ChebyshevReconstructor::reconstructRodShape(const Eigen::Matrix<double, 3, Eigen::Dynamic> &t_centerline_points)const
{
    //  Perform the cosine transform
    m_CN = m_DDCT * t_centerline_points.transpose();

    double s;                   //  Rod centerline coordinate
    double T0 = 1.0/sqrt(2.0);  //  First Chebyshev polynomial (Normalized)
    double Th;                  //  Other Chebyshev polynomial (Normalized)
    double circ_coord;          //  Coordinate on the unit circle
    Eigen::Vector3d reconstructed_point;
    for(unsigned int i=0; i<m_number_of_reconstruction_points; i++) {

        s = static_cast<double>(i*m_step);

        //  Compute for the first point
        const Eigen::VectorXd res = m_CN.row(0) * T0;
        reconstructed_point = m_CN.row(0) * T0;

        //  Obtain the others Chebyshev points
        for(unsigned int h=1; h<=m_number_of_Chebyshev_points-1;h++) {

            //  Define the corresponding point on the unit circle
            circ_coord = 2*s - 1;
            Th = boost::math::chebyshev_t(h, circ_coord);

            reconstructed_point += m_CN.row(h) * Th;

        }

        m_reconstructed_points.row(i) = reconstructed_point;

    }

    return m_reconstructed_points;
}




}   //  namespace Chebyshev
