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






typedef Eigen::Matrix<double, 6, 1> Vector6d ;



struct SE3Pose {
    SE3Pose()=default;

    SE3Pose(const Eigen::VectorXd &t_quaternion,
            const Eigen::VectorXd &t_position);


    std::string operator <<(const SE3Pose &t_other);

    Eigen::Matrix3d getRotationMatrix()const;

    Eigen::Matrix4d getSE3Pose()const;


    Eigen::Quaterniond m_quaternion { Eigen::Quaterniond(1, 0, 0, 0) };
    Eigen::Vector3d m_position { Eigen::Vector3d::Zero() };
};






struct Kinematics {
    Kinematics()=default;

    SE3Pose m_pose;

    Vector6d m_twist;
    Vector6d m_accelerations;
};


struct TangentKinematics {
    Vector6d m_Delta_zeta;

    Vector6d m_Delta_twist;

    Vector6d m_Delta_acceleration;
};

template<class DataType>
struct GeneralizedCoordinates {
    DataType m_q;
    DataType m_dot_q;
    DataType m_ddot_q;
};



struct Screw {
    Screw()=default;

    Screw(const Eigen::Vector3d &t_angular,
          const Eigen::Vector3d &t_linear);

    Vector6d operator=(Eigen::VectorXd)
    {
        return ( Vector6d() << m_angular, m_linear ).finished();
    }


    Eigen::Vector3d m_angular { Eigen::Vector3d::Zero() };
    Eigen::Vector3d m_linear { Eigen::Vector3d::Zero() };
};


Eigen::Matrix3d skew(const Eigen::Vector3d &t_v);



Eigen::MatrixXd ad(const Eigen::VectorXd &t_strain);



Eigen::MatrixXd Ad(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r);



Eigen::MatrixXd DeltaAd(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r);



Eigen::Quaterniond rotateAlongAxis(const double &t_angle, const Eigen::Vector3d &t_axis);



}   //  namespace LieAlgebra


#endif // LIE_ALGEBRA_UTILITIES_H
