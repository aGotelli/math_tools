
#include "math_tools/LieAlgebra/lie_algebra_utilities.hpp"

#include <benchmark/benchmark.h>

using namespace ::LieAlgebra;







int main(int argc, char *argv[])
{
    static constexpr unsigned int repetions = 50;

    const auto alpha_1 = M_PI/2;
    const auto alpha_2 = -M_PI/2;

    const Eigen::Matrix3d Rx_1 = ::LieAlgebra::getRx(alpha_1);
    const Eigen::Matrix3d Rx_2 = ::LieAlgebra::getRx(alpha_2);


    ::benchmark::RegisterBenchmark("Difference in SO3", [&](benchmark::State &t_state){

        while(t_state.KeepRunning())
            const auto res = differenceInSO3(Rx_1, Rx_2);

    })->Repetitions(repetions);

    ::benchmark::RegisterBenchmark("log in SO3", [&](benchmark::State &t_state){

        while(t_state.KeepRunning())
            const auto res = logSO3(Rx_1, Rx_2);

    })->Repetitions(repetions);



    ::benchmark::Initialize(&argc, argv);

    ::benchmark::RunSpecifiedBenchmarks();

    return 0;
}
