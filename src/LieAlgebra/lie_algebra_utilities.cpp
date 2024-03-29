#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <eigen3/unsupported/Eigen/KroneckerProduct>

#include <iostream>

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


SE3Pose::SE3Pose(const Eigen::Matrix3d &t_R,
                 const Eigen::Vector3d &t_position)
    : m_quaternion( Eigen::Quaterniond(t_R) ),
      m_position(t_position)
{}



SE3Pose::SE3Pose(const Eigen::Vector3d &t_position,
                 const double &t_roll,
                 const double &t_pitch,
                 const double &t_yaw)
    : m_quaternion( [&](){
            Eigen::Quaterniond q;
            q = Eigen::AngleAxisd(t_roll, Eigen::Vector3d::UnitX())
                    * Eigen::AngleAxisd(t_pitch, Eigen::Vector3d::UnitY())
                    * Eigen::AngleAxisd(t_yaw, Eigen::Vector3d::UnitZ());
            return q;
        }() ),
      m_position(t_position)
{}

SE3Pose::SE3Pose(const Eigen::Vector3d &t_position,
                 const double &t_theta,
                 const Eigen::Vector3d &t_axis)
    : m_quaternion( rotateAlongAxis(t_theta, t_axis) ),
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

SE3Pose SE3Pose::Identity()
{
    const Eigen::Quaterniond q(1, 0, 0, 0);
    const auto r = Eigen::Vector3d::Zero();

    return SE3Pose(q, r);
}

SE3Pose SE3Pose::inverse() const
{
    return SE3Pose(m_quaternion.inverse(),
                   - m_quaternion.toRotationMatrix().transpose() * m_position);
}



std::string SE3Pose::toString(const std::string &t_indentation)const
{
    std::stringstream message;

    message << t_indentation << "Pose :\n";
    message << t_indentation << "    - position [x, y, z] :\n";
    message << t_indentation << "       [ " << m_position.x() << ", " << m_position.y() << ", " << m_position.z() << "]\n";
    message << t_indentation << "    - orientation : \n";
    message << t_indentation << "        - quaternion [w, x, y, z] :\n";
    message << t_indentation << "           [ " << m_quaternion.w() << ", "  << m_quaternion.x() << ", "  << m_quaternion.y() << ", "  << m_quaternion.z() << "]\n";

    const auto R = getRotationMatrix();
    message << t_indentation << "        - rotation matrix :\n";
    message << R;
    message << t_indentation << "\n\n";


    return message.str();
}


Eigen::Matrix3d SE3Pose::getRotationMatrix() const
{
    return m_quaternion.toRotationMatrix();
}

Eigen::Matrix4d SE3Pose::getSE3PoseAsMatrix() const
{  
    return (Eigen::Matrix4d() << getRotationMatrix(), m_position,
                                 Eigen::RowVector3d::Zero(), 1.0).finished();
}



//Screw::Screw(const Eigen::Vector3d &t_angular,
//             const Eigen::Vector3d &t_linear) : m_angular(t_angular),
//                                                m_linear(t_linear)
//{}




Kinematics::Kinematics(const SE3Pose &t_pose)
    : m_pose(t_pose)
{}


std::string Kinematics::toString(const std::string &t_indentation)const
{
    std::stringstream message;

    message << t_indentation << m_pose.toString(t_indentation);
    message << t_indentation << "Twist :\n";
    message << t_indentation << "    - angular [x, y, z] : [ " << m_twist(0) << ", " << m_twist(1) << ", " << m_twist(2) << "]\n";
    message << t_indentation << "    - linear  [x, y, z] : [ " << m_twist(3) << ", " << m_twist(4) << ", " << m_twist(5) << "]\n";
    message << t_indentation << "Acceleration :\n";
    message << t_indentation << "    - angular [x, y, z] : [ " << m_accelerations(0) << ", " << m_accelerations(1) << ", " << m_accelerations(2) << "]\n";
    message << t_indentation << "    - linear  [x, y, z] : [ " << m_accelerations(3) << ", " << m_accelerations(4) << ", " << m_accelerations(5) << "]\n";
    message << t_indentation << "\n\n";


    return message.str();
}



Kinematics Kinematics::operator*(const Kinematics &t_other)const
{
    Kinematics result;

    result.m_pose = this->m_pose*t_other.m_pose;

    result.m_twist = this->m_twist + t_other.m_twist;

    result.m_accelerations = this->m_accelerations + t_other.m_accelerations;

    return result;
}


Eigen::Matrix3d skew(const Eigen::Vector3d &t_v)
{

    Eigen::Matrix3d v_hat;
    v_hat <<  0   ,  -t_v(2),   t_v(1),
            t_v(2),     0   ,  -t_v(0),
           -t_v(1),   t_v(0),     0   ;

    return v_hat;
}


Eigen::Vector3d antiSkew(const Eigen::Matrix3d &t_skew_simmetric_matrix)
{


    return (Eigen::Vector3d() << -t_skew_simmetric_matrix(1, 2),
                                  t_skew_simmetric_matrix(0, 2),
                                 -t_skew_simmetric_matrix(0, 1)
            ).finished();
}



Matrix6d Rcal(const Eigen::Matrix3d &t_aR_b)
{
    return Eigen::KroneckerProduct(Eigen::Matrix2d::Identity(), t_aR_b);
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




Matrix6d ad(const Eigen::Vector3d t_angular_component,
            const Eigen::Vector3d t_linear_component)
{
    Matrix6d ad(6,6);
    ad << skew(t_angular_component),     Eigen::Matrix3d::Zero(),
          skew(t_linear_component) ,    skew(t_angular_component) ;

    return ad;
}




Matrix6d Ad(const Eigen::Matrix3d &t_aR_b,
            const Eigen::Vector3d &t_r_ba)
{
    //  Decompose the strain

    Eigen::Matrix<double, 6, 6> Ad;
    Ad <<       t_aR_b       , Eigen::Matrix3d::Zero(),
          skew(t_r_ba)*t_aR_b,         t_aR_b         ;

    return Ad;
}


Matrix6d Ad(const SE3Pose &t_ag_b)
{
    return Ad(t_ag_b.getRotationMatrix(),
              t_ag_b.m_position);
}



Matrix6d dotAd(const Eigen::Matrix3d &t_aR_b,
               const Eigen::Vector3d &t_r_ba,
               const Vector6d &t_eta_ab)
{
    return Ad(t_aR_b, t_r_ba) * ad(t_eta_ab);
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
                 const Eigen::Vector3d &t_r_ba,
                 const Vector6d &t_Delta_zeta_ab)
{
    return Ad(t_aR_b, t_r_ba)*ad(t_Delta_zeta_ab);
}



Matrix6d DeltaAd(const SE3Pose &t_ag_b,
                 const Vector6d &t_Delta_zeta_ab)
{
    return DeltaAd(t_ag_b.m_quaternion.toRotationMatrix(),
                   t_ag_b.m_position,
                   t_Delta_zeta_ab);
}



Matrix6d DeltaDotAd(const Eigen::Matrix3d &t_aR_b,
                    const Eigen::Vector3d &t_r_ba,
                    const Vector6d &t_Delta_zeta_ab,
                    const Vector6d &t_eta_ab,
                    const Vector6d &t_Delta_eta_ab)
{
    return Ad(t_aR_b, t_r_ba) * ( ad(t_Delta_zeta_ab)*ad(t_eta_ab) + ad(t_Delta_eta_ab) );
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



Eigen::Matrix3d getRx(const double &t_alpha)
{
    const auto c = cos(t_alpha);
    const auto s = sin(t_alpha);

    Eigen::Matrix3d R;
    R << 1,  0,  0,
         0,  c, -s,
         0,  s,  c;

    return R;
}



Eigen::Matrix3d getRy(const double &t_beta)
{
    const auto c = cos(t_beta);
    const auto s = sin(t_beta);

    Eigen::Matrix3d R;
    R << c,  0,  s,
         0,  1,  0,
        -s,  0,  c;

    return R;
}



Eigen::Matrix3d getRz(const double &t_theta)
{
    const auto c = cos(t_theta);
    const auto s = sin(t_theta);

    Eigen::Matrix3d R;
    R << c, -s,  0,
         s,  c,  0,
         0,  0,  1;

    return R;
}


Eigen::Matrix3d expRodigues(const Eigen::Matrix3d &t_Omega)
{
    const auto theta = antiSkew(t_Omega).norm();
    const auto omega = t_Omega / theta;

    const auto exp = Eigen::Matrix3d::Identity()
                    + omega * sin( theta )
                    + omega*omega * (1 - cos( theta ));

    if(exp.hasNaN())
        return Eigen::Matrix3d::Identity();

    return exp;
}



Eigen::Matrix3d expRodigues(const Eigen::Vector3d &t_Omega)
{
    const auto theta = t_Omega.norm();
    const auto omega = t_Omega / theta;


    if(omega.hasNaN())
        return Eigen::Matrix3d::Identity();


    return Eigen::AngleAxisd(theta, omega).toRotationMatrix();
//    const auto exp = Eigen::Matrix3d::Identity()
//                    + ( omega * sin( theta )
//                        +omega.cross(omega) * (1 - cos( theta )) ).cross(Eigen::Matrix3d::Identity());


//    return exp;
}





Eigen::Vector3d differenceInSO3(const Eigen::Matrix3d &t_Ra,
                                const Eigen::Matrix3d &t_Rb)
{
    return antiSkew(t_Ra.transpose()*t_Rb - t_Rb.transpose()*t_Ra);
}



Eigen::Vector3d logSO3(const Eigen::Matrix3d &t_aR_b)
{
    static_assert( std::numeric_limits<double>::has_quiet_NaN );

    //  Compute the angle
    const double theta = acos( 0.5*(t_aR_b.trace() - 1) );

    //  Compute the gain
    const double gain = [&](){
        const auto res = theta/(2*sin(theta));
        if(std::isnan(res))
            return 0.0;
        return res;
    }();

    return gain * antiSkew( t_aR_b - t_aR_b.transpose() );
}



Eigen::Vector3d logSO3(const Eigen::Matrix3d &t_Ra,
                       const Eigen::Matrix3d &t_Rb)
{
    //  Compute resulting rotation matrix
    const Eigen::Matrix3d aR_b = t_Ra.transpose() * t_Rb;

    //  Call the log funtion
    return logSO3(aR_b);
}




std::pair<double, double> alpha_beta(const Eigen::Vector3d &t_Theta)
{
    const double norm = t_Theta.norm();

    double alpha = sin(norm)/norm;
    double beta = (2 - 2*cos(norm))/(norm*norm);

    if(std::isnan(alpha) || std::isnan(beta) ||
       std::isinf(alpha) || std::isinf(beta)){
        alpha = 1 - ((norm*norm)/6.0);
        beta =  1 - ((norm*norm)/12.0);
    }

    return {alpha, beta};
}


Eigen::Matrix3d expSO3(const Eigen::Vector3d &t_Theta)
{
    const auto [alpha,beta] = alpha_beta( t_Theta );

    const Eigen::Matrix3d eye = Eigen::Matrix3d::Identity();


    const Eigen::Matrix3d hat_Theta = skew(t_Theta);

    Eigen::Matrix3d exp_SO3 = eye + alpha*hat_Theta + (beta/2.0)*hat_Theta*hat_Theta;

    return exp_SO3;

}



Eigen::Matrix3d TSO3(const Eigen::Vector3d &t_Theta)
{
    const auto [alpha,beta] = alpha_beta( t_Theta );

    const Eigen::Matrix3d eye = Eigen::Matrix3d::Identity();

    const double norm = t_Theta.norm();

    const Eigen::Matrix3d hat_Theta = skew(t_Theta);

    Eigen::Matrix3d T_SO3 = eye
                            - ( beta/2.0 )*hat_Theta
                            + ( (1.0-alpha)/(norm*norm) )*hat_Theta*hat_Theta;

    if(T_SO3.hasNaN())
        T_SO3 = eye - (1.0/2.0)*hat_Theta + (1.0/6.0)*hat_Theta*hat_Theta;

    return T_SO3;
}


Eigen::Matrix3d TSO3_inverse(const Eigen::Vector3d &t_Theta)
{
    const auto [alpha,beta] = alpha_beta( t_Theta );

    const Eigen::Matrix3d eye = Eigen::Matrix3d::Identity();

    const double norm = t_Theta.norm();

    const Eigen::Matrix3d hat_Theta = skew(t_Theta);


    Eigen::Matrix3d T_SO3_inverse = eye
                                    + ( 1.0/2.0 )*hat_Theta
                                    + ( 1.0/(norm*norm) ) * (1.0 - (alpha/beta))*hat_Theta*hat_Theta;


    if(T_SO3_inverse.hasNaN())
        T_SO3_inverse = eye;

    return T_SO3_inverse;
}



Eigen::Matrix3d DeltaTSO3(const Eigen::Vector3d &t_Theta,
                          const Eigen::Vector3d &t_Delta_Theta)
{

    const auto theta = t_Theta.norm();


    Eigen::Matrix3d Delta_TSO3;

    if(abs( theta ) <= 1e-5) [[unlikely]]
            Delta_TSO3 = -0.5*skew(t_Delta_Theta);
    else [[likely]] {
        const auto U = t_Theta/theta;
        const auto UT = U.transpose();


        const Eigen::Matrix3d eye = Eigen::Matrix3d::Identity();

        const Eigen::Vector3d f1 = t_Delta_Theta/theta;

        const double f2 = sin(0.5*theta)/(0.5*theta);

        const double f3 = f2*f2;

        const double cth = cos(theta);

        const double Uf1 = U.dot(f1);

        const auto [alpha, beta] = alpha_beta( t_Theta );

        Eigen::Matrix3d C1 = (cth-alpha)*Uf1*eye;
        Eigen::Matrix3d C2 = (1-alpha)*(t_Delta_Theta*UT-U*t_Delta_Theta.transpose());
        Eigen::Matrix3d C3 = (3*alpha-cth-2)*Uf1*(U*UT);
        Eigen::Matrix3d C4 = (f3-alpha)*Uf1*skew(t_Theta);
        Eigen::Matrix3d C5 = 0.5*f3*skew(t_Delta_Theta);

        Delta_TSO3 = C1 + C2 + C3 + C4 - C5;





    }



    return Delta_TSO3;


}





}   //  namespace LieAlgebra
