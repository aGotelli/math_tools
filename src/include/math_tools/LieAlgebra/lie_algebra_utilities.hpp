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

/// \namespace LieAlgebra contains all the function related to the Lie Algebra theory
namespace LieAlgebra {


/// Define a scalar which is compatible with Eigen::VectorXd
typedef Eigen::Matrix<double, 0, 1> Vector0d ;


/// Define a scalar which is compatible with Eigen::VectorXd
typedef Eigen::Matrix<double, 1, 1> Vector1d ;


/// Define a vector of 6 elements
typedef Eigen::Matrix<double, 6, 1> Vector6d ;


/// Define a 6x6 matrix
typedef Eigen::Matrix<double, 6, 6> Matrix6d ;



/*!
 * \brief The SE3Pose struct embeds all the functionalities describing a configuration, or pose, in SE(3)
 *
 * This class is basically a warper around a quaternion and a position vector.
 * It provides all the conversions to rotation matrix and rigid body transformation that we need, as
 * well as the inverse of the transformation matrix.
 */
struct SE3Pose {
    SE3Pose()=default;

    /*!
     * \brief SE3Pose constructs the object using only the position
     * \param t_position the user defined position of the SE(3) pose
     *
     * This constructor assign the position to the pose leaving the orientation as identity
     */
    SE3Pose(const Eigen::Vector3d &t_position);


    /*!
     * \brief SE3Pose constructs the object using an Eigen::Quaterniond and position
     * \param t_quaternion the user defined orientation of the SE(3) pose, using Eigen::Quaterniond
     * \param t_position the user defined position of the SE(3) pose
     */
    SE3Pose(const Eigen::Quaterniond &t_quaternion,
            const Eigen::Vector3d &t_position);

    /*!
     * \brief SE3Pose constructs the object using a quaternion and position as Eigen::Vector
     * \param t_quaternion the user defined orientation of the SE(3) pose, using Eigen::Vector4d
     * \param t_position the user defined position of the SE(3) pose
     *
     * In this constructor, the quaternion is given as a vector of 4.
     * Given the quaternion with the definition \f$ q = w + x \textbf{i} + y \textbf{j} + z \textbf{k} \f$,
     * the vector should then be filled as \f$ t\_quaternion = \begin{bmatrix} w    &   x   &   y   &   z \end{bmatrix} \f$
     */
    SE3Pose(const Eigen::Vector4d &t_quaternion,
            const Eigen::Vector3d &t_position);


    /*!
     * \brief SE3Pose constructs the object using a rotation matrix and a position vector
     * \param t_R the user defined orientation of the SE(3) pose as a rotation matrix
     * \param t_position the user defined position of the SE(3) pose
     *
     */
    SE3Pose(const Eigen::Matrix3d &t_R,
            const Eigen::Vector3d &t_position);


    /*!
     * \brief SE3Pose constructs the object using
     * \param t_position the user defined position of the SE(3) pose
     * \param t_roll the roll angle of rotation along the x asis
     * \param t_pitch the pitch angle of rotation along the y asis
     * \param t_yaw the yaw angle of rotation along the z asis
     *
     * This constructor defines an SE(3) pose from the given position and the given agles of rotations.
     * In this constructor we use the Euler angles with the convention XYZ
     */
    SE3Pose(const Eigen::Vector3d &t_position,
            const double &t_roll,
            const double &t_pitch,
            const double &t_yaw);

    /*!
     * \brief SE3Pose constructs the object using
     * \param t_position the user defined position of the SE(3) pose
     * \param t_theta the angle of rotation
     * \param t_axis the axis of rotation
     *
     * This constructor uses the angle and the given axis to define the orientation in SO(3)
     * for the pose.
     */
    SE3Pose(const Eigen::Vector3d &t_position,
            const double &t_theta,
            const Eigen::Vector3d &t_axis=Eigen::Vector3d::UnitZ());



    static SE3Pose Identity();

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
    std::string toString(const std::string &t_indentation="")const;

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

    ///  The quaternion representing SO(3)
    Eigen::Quaterniond m_quaternion { Eigen::Quaterniond(1, 0, 0, 0) };

    ///  The vector of position belonging to R^3
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

    /*!
     * \brief Kinematics defines the Kinematics object
     * \param t_pose the pose in SE(3) of the body
     */
    Kinematics(const SE3Pose &t_pose);

    /*!
     * \brief toString gives a string containing all the information about the body kinematics
     * \return the string containing all the information about the body kinematics
     */
    std::string toString(const std::string &t_indentation="")const;


    /*!
     * \brief operator * defines the multiplication of two kinematics configurations
     * \param t_other The other kinematic configuration
     * \return the result of the multiplication of t_other kinematics by this pose
     *
     * This operator is used to multiply directly two elements of type Kinematics
     */
    Kinematics operator*(const Kinematics &t_other)const;


    /// The pose in SE(3) of the body
    SE3Pose m_pose;

    /// The twist vector in se(3)
    Vector6d m_twist { Vector6d::Zero() };

    /// The accelerations vector in se(3)
    Vector6d m_accelerations { Vector6d::Zero() };
};


/*!
 * \brief The TangentKinematics struct contains all the elements to describe the tangent kinematics of a body
 *
 * This class collects three Vector6 in order to describe tangent kinematics of a body
 * in terms of its variation of pose, the tangent velocities and tangent accelerations
 */
struct TangentKinematics {

    /// The vector for \f$ \Delta \zeta \f$
    Vector6d m_Delta_zeta { Vector6d::Zero() };

    /// The vector for \f$ \Delta \eta \f$
    Vector6d m_Delta_twist { Vector6d::Zero() };

    /// The vector for \f$ \Delta \dot{\eta} \f$
    Vector6d m_Delta_acceleration { Vector6d::Zero() };
};



/*!
 * \brief skew computes the skew simmetric matrix associated to a vector in R^3
 * \param t_v the vector in R^3
 * \return the skew simmetric matrix associated to a vector in R^3
 */
Eigen::Matrix3d skew(const Eigen::Vector3d &t_v);


/*!
 * \brief antiSkew gives the inverse of the skew operator
 * \param t_skew_simmetric_matrix the skew simmetrix matrix to convert into a vector
 * \return the inverse of the skew operator
 *
 * This function takes a skew simmetric matrix and returns the corresponding vector.
 * The vector is extracted directly from the skew simmetric matrix
 */
Eigen::Vector3d antiSkew(const Eigen::Matrix3d &t_skew_simmetric_matrix);



/*!
 * \brief Rcal creates a block-diagonal matrix that projects a screw from one frame to another
 * \param t_aR_b the rotation matrix from the frame a to the frame b
 * \return a block-diagonal matrix that projects a screw from one frame to another
 */
Matrix6d Rcal(const Eigen::Matrix3d &t_aR_b);



/*!
 * \brief ad computes the adjoint transformation associated to a twist
 * \param t_twist the twist to transform
 * \return the adjoint transformation associated to a twist
 */
Matrix6d ad(const Vector6d &t_twist);


/*!
 * \brief ad computes the adjoint transformation associated to a twist
 * \param t_angular_component is the angular component of the twist
 * \param t_linear_component is the linear component of the twist
 * \return  the adjoint transformation associated to a twist
 */
Matrix6d ad(const Eigen::Vector3d t_angular_component,
            const Eigen::Vector3d t_linear_component);


/*!
 * \brief Ad computes the Adjoint transformation associated to a homogeneous transformation
 * \param t_aR_b the orientation of frame b with respect of frame a
 * \param t_r_ba the position of frame b with respect of frame a
 * \return the adjoint transformation associated to a homogeneous transformation
 *
 * This function computes the Adjoint transformation \f$Ad({^a\textbf{g}_b} )\f$.
 * With this function, aR_b expresses the orientation of a frame \f$\mathcal{F}_b\f$
 * with respect of a frame \f$\mathcal{F}_a\f$, while ar_b expresses the position of the
 * frame \f$\mathcal{F}_b\f$ with respect to the frame \f$\mathcal{F}_a\f$.
 */
Matrix6d Ad(const Eigen::Matrix3d &t_aR_b,
            const Eigen::Vector3d &t_r_ba);


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
 * \param t_r_ba the position of frame b with respect of frame a
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \return the derivative with respect to time of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the derivative with respect to time of the Adjoint transformation \f$Ad({^a\textbf{g}_b} )\f$.
 * We thus compute the map \f$\dot{Ad}({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} ) ad(\eta_{a/b})\f$.
 * With this function, aR_b expresses the orientation of a frame \f$\mathcal{F}_b\f$
 * with respect of a frame \f$\mathcal{F}_a\f$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_eta_ab expresses the relative twist between the frame \f$\mathcal{F}_b\f$ and the frame \f$\mathcal{F}_a\f$.
 * In other words: the twist of the frame \f$\mathcal{F}_b\f$ with respect to the frame \f$\mathcal{F}_a\f$.
 */
Matrix6d dotAd(const Eigen::Matrix3d &t_aR_b,
               const Eigen::Vector3d &t_r_ba,
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
 * \param t_r_ba the position of frame b with respect of frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \return the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the increment of the Adjoint transformation \f$Ad({^a\textbf{g}_b} )\f$.
 * We thus compute the map \f$\Delta Ad({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} ) ad(\Delta \zeta_{a/b})\f$.
 * With this function, aR_b expresses the orientation of a frame \f$\mathcal{F}_b\f$
 * with respect of a frame \f$\mathcal{F}_a\f$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_Delta_zeta_ab expresses the increment in the SE(3) pose associated to the homogeneous transformation
 * of frame b with respect to frame a.
 * In other words: it expresses the variation of the pose of frame \f$\mathcal{F}_b\f$ with respect to frame \f$\mathcal{F}_a\f$.
 */
Matrix6d DeltaAd(const Eigen::Matrix3d &t_aR_b,
                 const Eigen::Vector3d &t_r_ba,
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
 * \param t_r_ba the position of frame b with respect of frame a
 * \param t_Delta_zeta_ab the increment in the SE(3) pose associated to the homogeneous transformation
 * \param t_eta_ab the twist of the frame b with respect to the frame a
 * \param t_Delta_eta_ab the increment in se(3) for the twist of the frame b with respect to the frame a
 * \return the derivative, with respect to time, of the increment of the Adjoint transformation associated to the homogeneous transformation
 *
 * This function computes the derivative, with respect to time, of the increment of the Adjoint transformation \f$Ad({^a\textbf{g}_b} )\f$.
 * We thus compute the map \f$\Delta \dot{Ad}({^a\textbf{g}_b} ) = Ad({^a\textbf{g}_b} )
 * \left( ad(\Delta \zeta_{a/b})ad(\eta_{a/b}) + ad(\Delta \eta_{a/b}) \right)\f$.
 * With this function, aR_b expresses the orientation of a frame \f$\mathcal{F}_b\f$
 * with respect of a frame \f$\mathcal{F}_a\f$, while ar_b expresses the position of the
 * frame b with respect to the frame a.
 * The parameter t_Delta_zeta_ab expresses the increment in the SE(3) pose associated to the homogeneous transformation
 * of frame b with respect to frame a.
 * In other words: it expresses the variation of the pose of frame \f$\mathcal{F}_b\f$ with respect to frame \f$\mathcal{F}_a\f$.
 * The parameter t_eta_ab expresses the relative twist between the frame \f$\mathcal{F}_b\f$ and the frame \f$\mathcal{F}_a\f$.
 * In other words: the twist of the frame \f$\mathcal{F}_b\f$ with respect to the frame \f$\mathcal{F}_a\f$.
 * The parameter t_Delta_eta_ab expresses the increment in the relative twist between the frame \f$\mathcal{F}_b\f$ and the frame \f$\mathcal{F}_a\f$.
 * In other words: it expresses the variation of the twist of the frame \f$\mathcal{F}_b\f$ with respect to the frame \f$\mathcal{F}_a\f$.
 */
Matrix6d DeltaDotAd(const Eigen::Matrix3d &t_aR_b,
                    const Eigen::Vector3d &t_r_ba,
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
 * \param t_angle is the angle of rotation in radians [RAD]
 * \param t_axis is the normalised rotation axis
 * \return the quaternion rotated along the given axis
 */
Eigen::Quaterniond rotateAlongAxis(const double &t_angle,
                                   const Eigen::Vector3d &t_axis);



/*!
 * \brief getRx compute the rotation matrix associated to a rotation along the x axis
 * \param t_alpha the angle of rotation along the x axis in radians [RAD]
 * \return the rotation matrix associated to a rotation along the x axis
 */
Eigen::Matrix3d getRx(const double &t_alpha);


/*!
 * \brief getRy compute the rotation matrix associated to a rotation along the y axis
 * \param t_beta the angle of rotation along the y axis in radians [RAD]
 * \return the rotation matrix associated to a rotation along the y axis
 */
Eigen::Matrix3d getRy(const double &t_beta);


/*!
 * \brief getRz compute the rotation matrix associated to a rotation along the z axis
 * \param t_theta the angle of rotation along the z axis in radians [RAD]
 * \return the rotation matrix associated to a rotation along the z axis
 */
Eigen::Matrix3d getRz(const double &t_theta);



Eigen::Matrix3d expRodigues(const Eigen::Matrix3d &t_Omega);



Eigen::Matrix3d expRodigues(const Eigen::Vector3d &t_Omega);


/*!
 * \brief differenceInSO3 computes the difference in SO(3) between two orientation matrices
 * \param t_Ra the orientation matrix of frame a with respect to the reference frame
 * \param t_Rb the orientation matrix of frame b with respect to the reference frame
 * \return the vector in R^3 containing the difference in SO(3) between two orientation matrices
 *
 * This function takes two orientation, or rotation, matrices and computes the corresponding difference in SO(3).
 * These two rotation matrices are expressed with respect to the same reference frame.
 * The function computes the difference in the form of a vector, as solution of the following equation
 *
 * \f$\left[R_a^T R_b - R_a R_b^T \right]^\Vee\f$
 *
 */
Eigen::Vector3d differenceInSO3(const Eigen::Matrix3d &t_Ra,
                                const Eigen::Matrix3d &t_Rb);




/*!
 * \brief logSO3 compute the log in SO(3) between two rotation matrices
 * \param t_aR_b the rotation matrix expressing the orientation of frame b with respect to frame a
 * \return the log in SO(3) of the rotation matrix
 *
 * This function takes a relative rotations matrix as parameter and returns the error in SO(3) computed with
 * the log.
 *
 * It computes the error as :
 *
 * \f$ log_{SO(3)} \left( t\_{^aR_b} \right)  \f$
 */
Eigen::Vector3d logSO3(const Eigen::Matrix3d &t_aR_b);


/*!
 * \brief logSO3 compute the log in SO(3) between two rotation matrices
 * \param t_Ra the rotation matrix expressing the orientation of a frame a with respect to the reference frame
 * \param t_Rb the rotation matrix expressing the orientation of a frame b with respect to the reference frame
 * \return the log in SO(3) between two rotation matrices
 *
 * This function takes two rotations matrices as parameter and returns the error in SO(3) computed with
 * the log.
 *
 * Starting from t_Ra and t_Rb it computes the error as :
 *
 * \f$ log_{SO(3)} \left( t\_R_a^T t\_R_b \right) = log_{SO(3)} \left( t\_{^aR_b} \right) \f$
 *
 * This function serves as warper to the function logSO3(const Eigen::Matrix3d&)
 */
Eigen::Vector3d logSO3(const Eigen::Matrix3d &t_Ra,
                       const Eigen::Matrix3d &t_Rb);




std::pair<double, double> alpha_beta(const Eigen::Vector3d &t_Theta);

Eigen::Matrix3d expSO3(const Eigen::Vector3d &t_Theta);

Eigen::Matrix3d TSO3(const Eigen::Vector3d &t_Theta);

Eigen::Matrix3d TSO3_inverse(const Eigen::Vector3d &t_Theta);


Eigen::Matrix3d DeltaTSO3(const Eigen::Vector3d &t_Theta,
                          const Eigen::Vector3d &t_Delta_Theta);


}   //  namespace LieAlgebra


#endif // LIE_ALGEBRA_UTILITIES_H
