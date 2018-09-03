#include <iostream>
#include <random>
#include "math/vector3d.hpp"
#include "math/const.hpp"
#include "orbit/orbit.hpp"
#include "astro/lambert_problem.hpp"
#include "../tests/astro/lambert_test.hpp"
#include "../tests/astro/state_vectors_propagation_test.hpp"
#include "../tests/astro/orbit_test.hpp"
#include "../tests/math/acos_test.hpp"

int main()
{
    using namespace astro_cpp;
    std::ios::sync_with_stdio(false);
    std::cout.precision(16);

    //tests::kepler_solver_test(1000 * 1000);
    //tests::acos_test(10);
    //tests::test_case_size();
    //tests::lambert::lambert_problem(200 * 1000);

    tests::state_vectors_propagation(1000'000);
         // 4.54876098624392e-06

    return 0;
}
