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


/// Define a scalar which is compatible with Eigen::VectorXd
typedef Eigen::Matrix<double, 1, 1> Vector1d ;


/// Define a vector of 6 elements
typedef Eigen::Matrix<double, 6, 1> Vector6d ;


/// Define a 6x6 matrix
typedef Eigen::Matrix<double, 6, 6> Matrix6d ;



/*!
 * \brief The SE3Pose struct embeds all the functionalities describing a configuration, or pose, in SE(3)
 */
struct SE3Pose {
    SE3Pose()=default;

    SE3Pose(const Eigen::Vector3d &t_position);


    SE3Pose(const Eigen::Quaterniond &t_quaternion,
            const Eigen::Vector3d &t_position);


    SE3Pose(const Eigen::Vector4d &t_quaternion,
            const Eigen::Vector3d &t_position);

    /*!
     * \brief operator * defines the multiplication of a pose with another one
     * \param t_other The other pose to be bultiplied
     * \return the result of the multiplication of t_other pose by this pose
     *
     * This operator is used to multiply directly two elements of type SE3Pose
     */
    SE3Pose operator*(const SE3Pose &t_other)const;

    /*!
     * \brief inverse computes the inverse of the current transformation
     * \return the inverse of this transformation
     */
    SE3Pose inverse()const;

    /*!
     * \brief toString creates a string containing all the information about the pose
     * \return the string containing all the information about the pose
     */
    std::string toString()const;

    /*!
     * \brief getRotationMatrix gives the rotation matrix associated to the SE(3) pose
     * \return the rotation matrix associated to the SE(3) pose
     */
    Eigen::Matrix3d getRotationMatrix()const;

    /*!
     * \brief getSE3PoseAsMatrix gives the 4x4 matrix corresponding to the class
     * \return the 4x4 matrix corresponding to the class
     */
    Eigen::Matrix4d getSE3PoseAsMatrix()const;

    //  The quaternion representing SO(3)
    Eigen::Quaterniond m_quaternion { Eigen::Quaterniond(1, 0, 0, 0) };

    //  The vector of position belonging to R^3
    Eigen::Vector3d m_position { Eigen::Vector3d::Zero() };
};





/*!
 * \brief The Kinematics struct contains all the elements to describe the kinematics of a body
 *
 * This class collects the SE3Pose and two Vector6 in order to describe the kinematics of a body
 * in terms of its pose, velocities and accelerations
 */
struct Kinematics {
    Kinematics()=default;


    Kinematics(const SE3Pose &t_pose);

    /*!
     * \brief toString gives a string containing all the information about the body kinematics
     * \return the string containing all the information about the body kinematics
     */
    std::string toString()const;

    SE3Pose m_pose;

    Vector6d m_twist;
    Vector6d m_accelerations;
};


/*!
 * \brief The TangentKinematics struct contains all the elements to describe the tangent kinematics of a body
 *
 * This class collects three Vector6 in order to describe tangent kinematics of a body
 * in terms of its variation of pose, the tangent velocities and tangent accelerations
 */
struct TangentKinematics {
    Vector6d m_Delta_zeta;

    Vector6d m_Delta_twist;

    Vector6d m_Delta_acceleration;
};


/*!
 * \brief The GeneralizedCoordinates struct utility class to declare generalised coordinates with a template data type
 */
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


/*!
 * \brief skew computes the skew simmetric matrix associated to a vector in R^3
 * \param t_v the vector in R^3
 * \return the skew simmetric matrix associated to a vector in R^3
 */
Eigen::Matrix3d skew(const Eigen::Vector3d &t_v);


/*!
 * \brief ad computes the adjoint transformation associated to a twist
 * \param t_twist the twist to transform
 * \return the adjoint transformation associated to a twist
 */
Matrix6d ad(const Vector6d &t_twist);


/*!
 * \brief Ad computes the Adjoint transformation associated to a homogeneous transformation
 * \param t_aR_b the orientation of frame b with respect of frame a
 * \param t_ar_b the position of frame b with respect of frame a
 * \return the adjoint transformation associated to a homogeneous transformation
 *
 * This function computes the Adjoint transformation $Ad({^a\textbf{g}_b} )$.
 * With this function, aR_b expresses the orientation of a frame $\mathcal{F}_b$
 * with respect of a frame $\mathcal{F}_a$, while ar_b expresses the position of the
 * frame $\mathcal{F}_b$ with respect to the frame $\mathcal{F}_a$.
 */
Matrix6d Ad(const Eigen::Matrix3d &t_aR_b,
            const Eigen::Vector3d &t_ar_b);


/*!
 * \brief Ad computes the Adjoint transformation associated to the homogeneous transformation
 * \param t_ag_b the homogeneous transformation expressing the frame b with respect to the frame a
 * \return
 *
 * This function is a warper for the function Ad(const Eigen::Matrix3d &,const Eigen::Vector3d &)
 */
Matrix6d Ad(const SE3Pose &t_ag_b);


/*!
 * \brief dotAd computes the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 * \param t_aR_b the orientation of frame b with respect of frame a
 * \param t_ar_b the position of frame b with respect of frame a
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \return the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the derivative with respect to time of the Adjoint transformation $Ad({^a\textbf{g}_b} )$.
 * We thus compute the map $\dot{Ad}({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} ) ad(\eta_{a/b})$.
 * With this function, aR_b expresses the orientation of a frame $\mathcal{F}_b$
 * with respect of a frame $\mathcal{F}_a$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_eta_ab expresses the relative twist between the frame $\mathcal{F}_b$ and the frame $\mathcal{F}_a$.
 * In other words: the twist of the frame $\mathcal{F}_b$ with respect to the frame $\mathcal{F}_a$.
 */
Matrix6d dotAd(const Eigen::Matrix3d &t_aR_b,
               const Eigen::Vector3d &t_ar_b,
               const Vector6d &t_eta_ab);


/*!
 * \brief dotAd computes the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 * \param t_ag_b the homogeneous transformation expressing the frame b with respect to the frame a
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \return the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function dotAd(const Eigen::Matrix3d &,const Eigen::Vector3d &,const Vector6d &)
 *
 */
Matrix6d dotAd(const SE3Pose &t_ag_b,
               const Vector6d &t_eta_ab);


/*!
 * \brief dotAd computes the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 * \param t_relative_kinematics_ab the relative kinematics of frame b with respect to frame a
 * \return the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function dotAd(const SE3Pose &,const Vector6d &)
 *
 */
Matrix6d dotAd(const Kinematics &t_relative_kinematics_ab);


/*!
 * \brief DeltaAd computes the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_aR_b the orientation of frame b with respect of frame a
 * \param t_ar_b the position of frame b with respect of frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \return the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the increment of the Adjoint transformation $Ad({^a\textbf{g}_b} )$.
 * We thus compute the map $\Delta Ad({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} ) ad(\Delta \zeta_{a/b})$.
 * With this function, aR_b expresses the orientation of a frame $\mathcal{F}_b$
 * with respect of a frame $\mathcal{F}_a$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_Delta_zeta_ab expresses the increment in the SE(3) pose associated to the homogeneous transformation
 * of frame b with respect to frame a.
 * In other words: it expresses the variation of the pose of frame $\mathcal{F}_b$ with respect to frame $\mathcal{F}_a$.
 */
Matrix6d DeltaAd(const Eigen::Matrix3d &t_aR_b,
                 const Eigen::Vector3d &t_ar_b,
                 const Vector6d &t_Delta_zeta_ab);


/*!
 * \brief DeltaAd computes the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_ag_b the homogeneous transformation expressing the frame b with respect to the frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \return the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function DeltaAd(const Eigen::Matrix3d &,const Eigen::Vector3d &,const Vector6d &t_Delta_zeta_ab)
 *
 */
Matrix6d DeltaAd(const SE3Pose &t_ag_b,
                 const Vector6d &t_Delta_zeta_ab);



/*!
 * \brief DeltaDotAd computes the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_aR_b the orientation of frame b with respect of frame a
 * \param t_ar_b the position of frame b with respect of frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \param t_Delta_eta_ab the increment in se(3) for the twist of the frame b with respect to the frame a
 * \return the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the derivative, with respect to time, of the increment of the Adjoint transformation $Ad({^a\textbf{g}_b} )$.
 * We thus compute the map $\Delta \dot{Ad}({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} )
 * \left( ad(\Delta \zeta_{a/b})ad(\eta_{a/b}) + ad(\Delta \eta_{a/b}) \right)$.
 * With this function, aR_b expresses the orientation of a frame $\mathcal{F}_b$
 * with respect of a frame $\mathcal{F}_a$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_Delta_zeta_ab expresses the increment in the SE(3) pose associated to the homogeneous transformation
 * of frame b with respect to frame a.
 * In other words: it expresses the variation of the pose of frame $\mathcal{F}_b$ with respect to frame $\mathcal{F}_a$.
 * The parameter t_eta_ab expresses the relative twist between the frame $\mathcal{F}_b$ and the frame $\mathcal{F}_a$.
 * In other words: the twist of the frame $\mathcal{F}_b$ with respect to the frame $\mathcal{F}_a$.
 * The parameter t_Delta_eta_ab expresses the increment in the relative twist between the frame $\mathcal{F}_b$ and the frame $\mathcal{F}_a$.
 * In other words: it expresses the variation of the twist of the frame $\mathcal{F}_b$ with respect to the frame $\mathcal{F}_a$.
 */
Matrix6d DeltaDotAd(const Eigen::Matrix3d &t_aR_b,
                    const Eigen::Vector3d &t_ar_b,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_eta_ab,
                    const Vector6d &t_Delta_eta_ab);


/*!
 * \brief DeltaDotAd computes the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_ag_b the homogeneous transformation expressing the frame b with respect to the frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \param t_Delta_eta_ab the increment in se(3) for the twist of the frame b with respect to the frame a
 * \return the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function DeltaDotAd(const Eigen::Matrix3d &,const Eigen::Vector3d &,const Vector6d &,const Vector6d &,const Vector6d &).
 */
Matrix6d DeltaDotAd(const SE3Pose &t_ag_b,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_eta_ab,
                    const Vector6d &t_Delta_eta_ab);


/*!
 * \brief DeltaDotAd computes the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_relative_kinematics_ab the relative kinematics of frame b with respect to frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \param t_Delta_eta_ab the increment in se(3) for the twist of the frame b with respect to the frame a
 * \return the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function DeltaDotAd(const SE3Pose &,const Vector6d &,const Vector6d &,const Vector6d &).
 */
Matrix6d DeltaDotAd(const Kinematics &t_relative_kinematics_ab,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_Delta_eta_ab);


/*!
 * \brief DeltaDotAd computes the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 * \param t_relative_kinematics_ab the relative kinematics of frame b with respect to frame a
 * \param t_relative_tangent_kinematics_ab the relative tangent, or incremental, kinematics of frame b with respect to frame a
 * \return the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function is a warper for the function DeltaDotAd(const Kinematics &,const Vector6d &,const Vector6d &).
 */
Matrix6d DeltaDotAd(const Kinematics &t_relative_kinematics_ab,
                    const TangentKinematics &t_relative_tangent_kinematics_ab);


/*!
 * \brief rotateAlongAxis rotates the quaternion along the given axis
 * \param t_angle is the angle of rotation
 * \param t_axis is the normalised rotation axis
 * \return the quaternion rotated along the given axis
 */
Eigen::Quaterniond rotateAlongAxis(const double &t_angle,
                                   const Eigen::Vector3d &t_axis);



}   //  namespace LieAlgebra


#endif // LIE_ALGEBRA_UTILITIES_H
