#include "lambert_test.hpp"

#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include "math/const.hpp"
#include "astro/orbit.hpp"
#include "astro/lambert.hpp"
#include "astro/lambert_original.hpp"

namespace astro_cpp::tests
{

    struct TestCase
    {
        TestCase(EllipticOrbit& orbit, double t1, double t2, StateVectors& state1, StateVectors& state2);

        EllipticOrbit orbit;
        double t1, t2, delta_t;
        StateVectors state1, state2;
    };

    TestCase::TestCase(EllipticOrbit& orbit, double t1, double t2, StateVectors& state1, StateVectors& state2):
        orbit(std::move(orbit)), t1(t1), t2(t2), delta_t(t2 - t1), state1(std::move(state1)), state2(std::move(state2))
    {

    }

    void test_case_size()
    {
        std::cout << sizeof(EllipticOrbit) << ' ' << sizeof(TestCase) << ' ' << sizeof(LambertProblem) << std::endl;
    }

    void lambert_problem(const int N)
    {
        // Init random generator
        std::random_device rand_device;
        std::mt19937 rand_gen(rand_device());
        std::uniform_real_distribution<double> a(0.0, 1.0e16),
                                               e(0.0, 1.0),
                                               i(0.0, math::PI_OVER_2),
                                               W(0.0, math::TAU),
                                               w(0.0, math::TAU);

        // Generate random test cases
        const double mu = 1.1723328e18,
                     M0 = 0.0,
                     t0 = 0.0;
        std::vector<TestCase> test_cases;
        using namespace std::chrono;

        auto reserve_start = high_resolution_clock::now();
        test_cases.reserve(N);
        auto reserve_end = high_resolution_clock::now();
        long long reserve_time = duration_cast<nanoseconds>(reserve_end - reserve_start).count();
        std::cout << "Reserving took " << reserve_time << " ns" << std::endl
                  << static_cast<double>(reserve_time) / N << " ns per element\n" << std::endl;

        auto gen_start = high_resolution_clock::now();
        for (int j = 0; j < N; ++j)
        {
            EllipticOrbit orbit(mu, a(rand_gen), e(rand_gen), i(rand_gen), W(rand_gen), w(rand_gen), M0, t0);
            std::uniform_real_distribution<double> time1(0.0, 0.50001 * orbit.period()),
                                                   time2(0.50002 * orbit.period(), orbit.period());
            double t1 = time1(rand_gen), t2 = time2(rand_gen);
            StateVectors state1 = orbit.state_vectors_from_time(t1),
                         state2 = orbit.state_vectors_from_time(t2);
            test_cases.emplace_back(orbit, t1, t2, state1, state2);
        }
        auto gen_end = high_resolution_clock::now();
        long long gen_time = duration_cast<nanoseconds>(gen_end - gen_start).count();
        std::cout << "Generation took " << gen_time << " ns" << std::endl
                  << static_cast<double>(gen_time) / N << " ns per problem" << std::endl;

        // Solve them with new Lambert solver
        std::cout << "\n\nLAMBERT TEST START!" << std::endl;
        std::vector<LambertProblem> lambert_opt;
        lambert_opt.reserve(N);
        using namespace std::chrono;
        auto lambert_start = high_resolution_clock::now();
        for (auto test_case: test_cases)
        {
            lambert_opt.emplace_back(test_case.state1.r, test_case.state2.r, test_case.delta_t, mu, false, 0);
        }
        auto lambert_end = high_resolution_clock::now();
        long long elapsed = duration_cast<nanoseconds>(lambert_end - lambert_start).count();
        std::cout << "Lambert new: " << elapsed << " ns" << std::endl
                  << static_cast<double>(elapsed) / N << " ns per lambert" << std::endl;


        /*
        // Solve them with old Lambert solver
        std::vector<lambert_original> lambert_izzo;
        lambert_izzo.reserve(N);
        //int failures = 0;
        using namespace std::chrono;
        auto lambert_start = high_resolution_clock::now();
        for (auto test_case: test_cases)
        {
            lambert_izzo.emplace_back(test_case.state1.r, test_case.state2.r, test_case.delta_t, mu, false, 0);
        }
        auto lambert_end = high_resolution_clock::now();
        std::cout << "Lambert old: " << duration_cast<microseconds>(lambert_end - lambert_start).count() << std::endl;
        //using namespace std::chrono_literals;
        //std::this_thread::sleep_for(10s);
        */

        int failures = 0;
        std::cout.precision(15);

        auto check_start = high_resolution_clock::now();
        for (int j = 0; j < N; ++j)
        {
            double delta_v1_opt = (test_cases[j].state1.v - lambert_opt[j].velocity.back().first).mag(),
                   delta_v2_opt = (test_cases[j].state2.v - lambert_opt[j].velocity.back().second).mag(),
                   //delta_v1_izzo = (test_cases[j].state1.v - lambert_izzo[j].get_v1().back()).mag(),
                   //delta_v2_izzo = (test_cases[j].state2.v - lambert_izzo[j].get_v2().back()).mag(),
                   v1_kepler = test_cases[j].state1.v.mag(),
                   v2_kepler = test_cases[j].state2.v.mag(),
                   v1_opt = lambert_opt[j].velocity.back().first.mag(),
                   v2_opt = lambert_opt[j].velocity.back().second.mag();
                   //v1_izzo = lambert_izzo[j].get_v1().back().mag(),
                   //v2_izzo = lambert_izzo[j].get_v2().back().mag();

            if ( ( (delta_v1_opt / v1_kepler) > 1e-6) ||
                 ( (delta_v2_opt / v2_kepler) > 1e-6) )
                 //( (delta_v1_izzo / v1_kepler) > 1e-6) ||
                 //( (delta_v2_izzo / v2_kepler) > 1e-6)   )
            {
                ++failures;
                test_cases[j].orbit.print();
                std::cout << "t1 = " << test_cases[j].t1 << " (" << test_cases[j].t1 / test_cases[j].orbit.period() << " T)\n";
                std::cout << "t2 = " << test_cases[j].t2 << " (" << test_cases[j].t2 / test_cases[j].orbit.period() << " T)\n";

                std::cout << "\n    v1_kepler     = " << v1_kepler     << ";  v2_kepler      = " << v2_kepler
                          << "\n    v1_opt        = " << v1_opt        << ";  v2_opt         = " << v2_opt
                          //<< "\n    v1_izzo       = " << v1_izzo       << ";  v2_izzo        = " << v2_izzo
                          << "\n    delta_v1_opt  = " << delta_v1_opt  << ";  delta_v2_opt   = " << delta_v2_opt
                         // << "\n    delta_v1_izzo = " << delta_v1_izzo << ";  delta_v2_izzo  = " << delta_v2_izzo
                          << std::endl;
            }
        }
        auto check_end = high_resolution_clock::now();
        long long check_time = duration_cast<nanoseconds>(check_end - check_start).count();
        std::cout << "\nCheck: " << check_time << " ns" << std::endl
                  << static_cast<double>(check_time) / N << " ns per check" << std::endl;
        std::cout << "\nFailures: " << static_cast<double>(failures) / N * 100 << " %";
    }
} // namespace astro_cpp::tests
