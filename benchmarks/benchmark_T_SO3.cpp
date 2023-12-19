#include <benchmark/benchmark.h>


#include <Eigen/Dense>
#include <iostream>

#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"


int main(int argc, char *argv[])
{
    const Eigen::Quaterniond q = Eigen::Quaterniond::UnitRandom();
//    std::cout << "q : \n" << q << "\n\n" << std::endl;

    const Eigen::Matrix3d R = q.toRotationMatrix();
//    std::cout << "R : \n" << R << "\n\n" << std::endl;

    const Eigen::Vector3d Theta = Eigen::Vector3d(-0.85911, -2.44664, 0.51629); //::LieAlgebra::logSO3(R);
//    std::cout << "Theta : \n" << Theta << "\n\n" << std::endl;

    const Eigen::Vector3d Delta_Theta = Eigen::Vector3d(0.59688, 0.823295, -0.604897); //Eigen::Vector3d::Random();
//    std::cout << "Delta_Theta : \n" << Delta_Theta << "\n\n" << std::endl;

    const Eigen::Matrix3d Delta_T_SO3 = ::LieAlgebra::DeltaTSO3(Theta, Delta_Theta);



    ::benchmark::RegisterBenchmark("T SO3", [&](::benchmark::State &t_state){

        Eigen::Matrix3d T_SO3;
        while(t_state.KeepRunning()){
            T_SO3 = ::LieAlgebra::TSO3(Theta);
            ::benchmark::DoNotOptimize( T_SO3 );
        }

    });


    ::benchmark::RegisterBenchmark("Delta T SO3", [&](::benchmark::State &t_state){

        Eigen::Matrix3d Delta_T_SO3;
        while(t_state.KeepRunning()){
            Delta_T_SO3 = ::LieAlgebra::DeltaTSO3(Theta, Delta_Theta);
            ::benchmark::DoNotOptimize( Delta_T_SO3 );
        }

    });


    ::benchmark::Initialize(&argc, argv);


    ::benchmark::RunSpecifiedBenchmarks();




    return 0;
}
