 #include "lie_algebra_utilities.hpp"










Eigen::Matrix3d skew(const Eigen::Vector3d &t_v)
{

    Eigen::Matrix3d v_hat;
    v_hat <<  0   ,  -t_v(2),   t_v(1),
            t_v(2),     0   ,  -t_v(0),
           -t_v(1),   t_v(0),     0   ;

    return v_hat;
}



Eigen::MatrixXd ad(const Eigen::VectorXd &t_strain)
{
    //  Decompose the strain
    const Eigen::Vector3d k = t_strain.block<3,1>(0,0);
    const Eigen::Vector3d gamma = t_strain.block<3,1>(3,0);

    Eigen::MatrixXd ad(6,6);
    ad << skew(k)    , Eigen::Matrix3d::Zero(),
          skew(gamma),          skew(k) ;

    return ad;
}




Eigen::MatrixXd Ad(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r)
{
    //  Decompose the strain

    Eigen::Matrix<double, 6, 6> Ad;
    Ad <<       t_R    , Eigen::Matrix3d::Zero(),
          skew(t_r)*t_R,          t_R ;

    return Ad;
}




Eigen::MatrixXd DeltaAd(const Eigen::Matrix3d &t_R,
                   const Eigen::Vector3d &t_r)
{
    //  Decompose the strain

    Eigen::Matrix<double, 6, 6> Ad;
    Ad <<       t_R    , Eigen::Matrix3d::Zero(),
          skew(t_r)*t_R,          t_R ;

    return Ad;
}




Eigen::Quaterniond rotateAlongAxis(const double &t_angle, const Eigen::Vector3d &t_axis)
{
    const auto axis_projection = sin(t_angle/2)*t_axis;
    return Eigen::Quaterniond(cos(t_angle/2), axis_projection[0], axis_projection[1], axis_projection[2]);
}
