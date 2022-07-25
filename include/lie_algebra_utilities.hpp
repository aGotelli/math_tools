/**
 * \file lie_algebra_utilities.hpp
 * \author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * \brief This files contains the formulas for the lie algebra
 * \version 0.1
 * \date 06-07-2022
 *
 * \copyright Copyright (c) 2022
 *
 */

#ifndef LIE_ALGEBRA_UTILITIES_H
#define LIE_ALGEBRA_UTILITIES_H

#include <Eigen/Dense>

namespace LieAlgebra {



Eigen::Matrix3d skew(const Eigen::Vector3d &t_v);



Eigen::MatrixXd ad(const Eigen::VectorXd &t_strain);



Eigen::MatrixXd Ad(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r);



Eigen::MatrixXd DeltaAd(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r);



Eigen::Quaterniond rotateAlongAxis(const double &t_angle, const Eigen::Vector3d &t_axis);


}   //  namespace LieAlgebra


#endif // LIE_ALGEBRA_UTILITIES_H
