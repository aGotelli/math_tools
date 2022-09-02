#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"





namespace LieAlgebra {


SE3Pose::SE3Pose(const Eigen::Vector3d &t_position)
    : m_position(t_position)
{}


SE3Pose::SE3Pose(const Eigen::Quaterniond &t_quaternion,
                 const Eigen::Vector3d &t_position)
    : m_quaternion( t_quaternion ),
      m_position(t_position)
{}


SE3Pose::SE3Pose(const Eigen::Vector4d &t_quaternion,
                 const Eigen::Vector3d &t_position)
    : m_quaternion(Eigen::Quaterniond(t_quaternion(0),
                                      t_quaternion(1),
                                      t_quaternion(2),
                                      t_quaternion(3))),
      m_position(t_position)
{}


SE3Pose SE3Pose::operator*(const SE3Pose &t_other) const
{

    SE3Pose result_pose;
    result_pose.m_quaternion = m_quaternion * t_other.m_quaternion;

    result_pose.m_position = m_quaternion.toRotationMatrix()*t_other.m_position
                            + m_position;

    return result_pose;
}



SE3Pose SE3Pose::inverse() const
{
    return SE3Pose(m_quaternion.inverse(),
                   - m_quaternion.toRotationMatrix().transpose() * m_position);
}



std::string SE3Pose::toString()const
{
    std::stringstream message;

    message << "Pose :\n";
    message << "    - position [x, y, z] : [ " << m_position.x() << ", " << m_position.y() << ", " << m_position.z() << "]\n";
    message << "    - orientation : \n";
    message << "        - quaternion [w, x, y, z] : [ " << m_quaternion.w() << ", "  << m_quaternion.x() << ", "  << m_quaternion.y() << ", "  << m_quaternion.z() << "]\n";

    const auto R = getRotationMatrix();
    message << "                            | " << R.row(0) << " |\n";
    message << "        - rotation matrix : | " << R.row(1) << " |\n";
    message << "                            | " << R.row(2) << " |\n";
//    message << "            " << getRotationMatrix();
    message << "\n\n";


    return message.str();
}


Eigen::Matrix3d SE3Pose::getRotationMatrix() const
{
    return m_quaternion.toRotationMatrix();
}

Eigen::Matrix4d SE3Pose::getSE3PoseAsMatrix() const
{
    return (Eigen::Matrix4d() << getRotationMatrix(), m_position.transpose(),
                                 Eigen::RowVector3d::Zero(), 1).finished();
}



Screw::Screw(const Eigen::Vector3d &t_angular,
             const Eigen::Vector3d &t_linear) : m_angular(t_angular),
                                                m_linear(t_linear)
{}




Kinematics::Kinematics(const SE3Pose &t_pose)
    : m_pose(t_pose)
{}


std::string Kinematics::toString()const
{
    std::stringstream message;

    message << m_pose.toString();
    message << "Twist :\n";
    message << "    - angular [x, y, z] : [ " << m_twist(0) << ", " << m_twist(1) << ", " << m_twist(2) << "]\n";
    message << "    - linear  [x, y, z] : [ " << m_twist(3) << ", " << m_twist(4) << ", " << m_twist(5) << "]\n";
    message << "Acceleration :\n";
    message << "    - angular [x, y, z] : [ " << m_accelerations(0) << ", " << m_accelerations(1) << ", " << m_accelerations(2) << "]\n";
    message << "    - linear  [x, y, z] : [ " << m_accelerations(3) << ", " << m_accelerations(4) << ", " << m_accelerations(5) << "]\n";
    message << "\n\n";


    return message.str();
}


Eigen::Matrix3d skew(const Eigen::Vector3d &t_v)
{

    Eigen::Matrix3d v_hat;
    v_hat <<  0   ,  -t_v(2),   t_v(1),
            t_v(2),     0   ,  -t_v(0),
           -t_v(1),   t_v(0),     0   ;

    return v_hat;
}



Matrix6d ad(const Vector6d &t_twist)
{
    //  Decompose the strain
    const Eigen::Vector3d angular = t_twist.block<3,1>(0,0);
    const Eigen::Vector3d linear = t_twist.block<3,1>(3,0);

    Matrix6d ad(6,6);
    ad << skew(angular),     Eigen::Matrix3d::Zero(),
          skew(linear) ,         skew(angular) ;

    return ad;
}




Matrix6d Ad(const Eigen::Matrix3d &t_aR_b,
            const Eigen::Vector3d &t_ar_b)
{
    //  Decompose the strain

    Eigen::Matrix<double, 6, 6> Ad;
    Ad <<       t_aR_b       , Eigen::Matrix3d::Zero(),
          skew(t_ar_b)*t_aR_b,         t_aR_b         ;

    return Ad;
}


Matrix6d Ad(const SE3Pose &t_ag_b)
{
    return Ad(t_ag_b.getRotationMatrix(),
              t_ag_b.m_position);
}



Matrix6d dotAd(const Eigen::Matrix3d &t_aR_b,
               const Eigen::Vector3d &t_ar_b,
               const Vector6d &t_eta_ab)
{
    return Ad(t_aR_b, t_ar_b) * ad(t_eta_ab);
}



Matrix6d dotAd(const SE3Pose &t_ag_b,
               const Vector6d &t_eta_ab)
{
    return dotAd(t_ag_b.m_quaternion.toRotationMatrix(),
                 t_ag_b.m_position,
                 t_eta_ab);
}



Matrix6d dotAd(const Kinematics &t_relative_kinematics_ab)
{
    const auto [pose, twist, _] = t_relative_kinematics_ab;

    return dotAd(pose, twist);
}



Matrix6d DeltaAd(const Eigen::Matrix3d &t_aR_b,
                 const Eigen::Vector3d &t_ar_b,
                 const Vector6d &t_Delta_zeta_ab)
{
    return Ad(t_aR_b, t_ar_b)*ad(t_Delta_zeta_ab);
}



Matrix6d DeltaAd(const SE3Pose &t_ag_b,
                 const Vector6d &t_Delta_zeta_ab)
{
    return DeltaAd(t_ag_b.m_quaternion.toRotationMatrix(),
                   t_ag_b.m_position,
                   t_Delta_zeta_ab);
}



Matrix6d DeltaDotAd(const Eigen::Matrix3d &t_aR_b,
                    const Eigen::Vector3d &t_ar_b,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_eta_ab,
                    const Vector6d &t_Delta_eta_ab)
{
    return Ad(t_aR_b, t_ar_b) * ( ad(t_Delta_zeta_ab)*ad(t_eta_ab) + ad(t_Delta_eta_ab) );
}



Matrix6d DeltaDotAd(const SE3Pose &t_ag_b,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_eta_ab,
                    const Vector6d &t_Delta_eta_ab)
{
    return DeltaDotAd(t_ag_b.m_quaternion.toRotationMatrix(),
                      t_ag_b.m_position,
                      t_Delta_zeta_ab,
                      t_eta_ab,
                      t_Delta_eta_ab);
}



Matrix6d DeltaDotAd(const Kinematics &t_relative_kinematics_ab,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_Delta_eta_ab)
{
    const auto [pose, twist, _] = t_relative_kinematics_ab;

    return DeltaDotAd(pose,
                      t_Delta_zeta_ab,
                      twist,
                      t_Delta_eta_ab);
}



Matrix6d DeltaDotAd(const Kinematics &t_relative_kinematics_ab,
                    const TangentKinematics &t_relative_tangent_kinematics_ab)
{
    return DeltaDotAd(t_relative_kinematics_ab,
                      t_relative_tangent_kinematics_ab.m_Delta_zeta,
                      t_relative_tangent_kinematics_ab.m_Delta_twist);
}



Eigen::Quaterniond rotateAlongAxis(const double &t_angle, const Eigen::Vector3d &t_axis)
{
    const auto axis_projection = sin(t_angle/2)*t_axis;
    return Eigen::Quaterniond(cos(t_angle/2), axis_projection[0], axis_projection[1], axis_projection[2]);
}


}   //  namespace LieAlgebra
