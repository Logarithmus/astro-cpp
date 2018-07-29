#include <iostream>
#include "math/vector3d.hpp"
#include "math/const.hpp"
#include "astro/orbit.hpp"
#include "astro/lambert.hpp"
#include "../tests/astro/lambert_test.hpp"
#include "../tests/astro/orbit_test.hpp"
#include "../tests/math/acos_test.hpp"

int main()
{
    using namespace astro_cpp;

    /*
    //Eeloo orbit
    EllipticOrbit test_orbit(1.1723328e18, 9.011882e10, 0.26,
                             6.15 * math::DEG_TO_RAD,
                             50 * math::DEG_TO_RAD,
                             260 * math::DEG_TO_RAD,
                             3.14, 0.0);
    test_orbit.print();
    StateVectors state1 = test_orbit.state_vectors_from_time(500 * 6 * 3600); // 500 Kerbin days
    state1.print();
    StateVectors state2 = test_orbit.state_vectors_from_time(1500 * 6 * 3600); //1500 Kerbin days
    state2.print();

    LambertProblem test_lambert(state1.r, state2.r, 1000 * 6 * 3600, 1.1723328e18, false, 0);
    for (auto v: test_lambert.velocity)
    {
        v.first.print();
        v.second.print();
    }

    for (auto it: test_lambert.iter_count)
    {
        std::cout << it << '\n';
    }
    */

    std::ios::sync_with_stdio(false);
    //tests::acos_test(10 * 1000 * 1000);
    //tests::kepler_solver_test(10 * 1000 * 1000);
    //tests::test_case_size();
    //for (int i = 0; i < 5; ++i)
    //{
        tests::lambert_problem(1000 * 1000);

    //}
    return 0;
}
