#include "state_vectors_propagation_test.hpp"

#include <random>
#include <chrono>
#include <vector>
#include "math/const.hpp"
#include "orbit/elliptic_orbit.hpp"

namespace astro_cpp::tests
{
    struct ValidationData
    {
        ValidationData(const EllipticOrbit& orb, double t1, double t2, const StateVectors& state2);

        EllipticOrbit orbit;
        double t1, t2;
        StateVectors state2;
    };

    ValidationData::ValidationData(const EllipticOrbit& orb, double t1, double t2, const StateVectors& state2):
        orbit(orb), t1(t1), t2(t2), state2(state2)
    {

    }

    struct TestCase
    {
        TestCase(const StateVectors& state1, double t1, double t2);

        StateVectors state1;
        double delta_t;
    };

    TestCase::TestCase(const StateVectors& state1, double t1, double t2):
        state1(state1), delta_t(t2 - t1)
    {

    }

    void state_vectors_propagation(const int N)
    {
        // Init random generator
        std::random_device rand_device;
        using namespace std::chrono;
        auto seed = high_resolution_clock::now().time_since_epoch().count();
        std::mt19937 rand_gen(seed);
        std::uniform_real_distribution<double> a(0.0, 1.0e13),
                                               e(0.1, 0.3),
                                               i(0.0, math::PI),
                                               W(0.0, math::TAU),
                                               w(0.0, math::TAU);

        // Generate random test cases
       // const double mu = 1.1723328e18,
        const double mu = 1.32712440018e20,
                     M0 = 0.0,
                     t0 = 0.0;
        std::vector<TestCase> test_cases;
        std::vector<ValidationData> validation_data;

        auto reserve_start = high_resolution_clock::now();
        test_cases.reserve(N);
        validation_data.reserve(N);
        auto reserve_end = high_resolution_clock::now();
        long long reserve_time = duration_cast<nanoseconds>(reserve_end - reserve_start).count();
        std::cout << "Reserving took " << reserve_time << " ns" << std::endl
                  << static_cast<double>(reserve_time) / N << " ns per element\n" << std::endl;

        auto gen_start = high_resolution_clock::now();
        for (int j = 0; j < N; ++j)
        {
            EllipticOrbit orbit(mu, a(rand_gen), e(rand_gen), i(rand_gen), W(rand_gen), w(rand_gen), M0, t0);
            std::uniform_real_distribution<double> time1(0.0, 1e-8 * orbit.period()),
                                                   time2(0.0, orbit.period());
            double t1 = time1(rand_gen), t2 = time2(rand_gen);
            StateVectors state1 = orbit.state_vectors_from_time(t1),
                         state2 = orbit.state_vectors_from_time(t2);
            validation_data.emplace_back(orbit, t1, t2, state2);
            test_cases.emplace_back(state1, t1, t2);
        }
        auto gen_end = high_resolution_clock::now();
        long long gen_time = duration_cast<nanoseconds>(gen_end - gen_start).count();
        std::cout << "Generation took " << gen_time << " ns" << std::endl
                  << static_cast<double>(gen_time) / N << " ns per problem" << std::endl;

        // State 1 propagation
        std::cout << "\n\nKEPLER PROPAGATION START !!!" << std::endl;
        std::vector<StateVectors> state2_propagated;
        state2_propagated.reserve(N);
        using namespace std::chrono;
        auto kepler_start = high_resolution_clock::now();
        for (auto test_case: test_cases)
        {
            state2_propagated.emplace_back(test_case.state1.propagated(mu, test_case.delta_t));
        }
        auto kepler_end = high_resolution_clock::now();
        long long elapsed = duration_cast<nanoseconds>(kepler_end - kepler_start).count();
        std::cout << "Propagation: " << elapsed << " ns" << std::endl
                  << static_cast<double>(elapsed) / N << " ns per problem" << std::endl;


        // Check results
        int failures = 0;
        std::cout.precision(15);

        double max_delta_r2 = 0.0;
        unsigned int j_max = 0;
        auto check_start = high_resolution_clock::now();
        for (int j = 0; j < N; ++j)
        {
            const double delta_r2 = (validation_data[j].state2.r - state2_propagated[j].r).mag();
               //          delta_v2 = (validation_data[j].state2.v - state2_propagated[j].v).mag(),
             //            r2_orbit = validation_data[j].state2.r.mag(),
             //            v2_orbit = validation_data[j].state2.v.mag(),
            //             r2_prop = state2_propagated[j].r.mag(),
            //             v2_prop = state2_propagated[j].v.mag();

            /*
            if ( ( (delta_r2 / r2_orbit) > 1e-10) ||
                 ( (delta_v2 / v2_orbit) > 1e-10) )
            {
                ++failures;
                validation_data[j].orbit.print();
                std::cout << "t1 = " << validation_data[j].t1 << " (" << validation_data[j].t1 / validation_data[j].orbit.period() << " T)\n";
                std::cout << "t2 = " << validation_data[j].t2 << " (" << validation_data[j].t2 / validation_data[j].orbit.period() << " T)\n";

                std::cout << "\n    r2_orbit = " << r2_orbit << ";  v2_orbit = " << v2_orbit
                          << "\n    r2_prop  = " << r2_prop  << ";  v2_prop  = " << v2_prop
                          << "\n    delta_r2 = " << delta_r2 << ";  delta_v2 = " << delta_v2
                          << std::endl;

            }
            */

            if (delta_r2 > max_delta_r2)
            {
                max_delta_r2 = delta_r2;
                j_max = j;
            }
        }
        auto check_end = high_resolution_clock::now();
        long long check_time = duration_cast<nanoseconds>(check_end - check_start).count();
        std::cout << "\nCheck: " << check_time << " ns" << std::endl
                  << static_cast<double>(check_time) / N << " ns per check" << std::endl;
        std::cout << "\nMax delta_r2: " << max_delta_r2 << std::endl;

        const double delta_r2 = (validation_data[j_max].state2.r - state2_propagated[j_max].r).mag(),
                     delta_v2 = (validation_data[j_max].state2.v - state2_propagated[j_max].v).mag(),
                     r2_orbit = validation_data[j_max].state2.r.mag(),
                     v2_orbit = validation_data[j_max].state2.v.mag(),
                     r2_prop = state2_propagated[j_max].r.mag(),
                     v2_prop = state2_propagated[j_max].v.mag();

        validation_data[j_max].orbit.print();
        std::cout << "t1 = " << validation_data[j_max].t1 << " (" << validation_data[j_max].t1 / validation_data[j_max].orbit.period() << " T)\n";
        std::cout << "t2 = " << validation_data[j_max].t2 << " (" << validation_data[j_max].t2 / validation_data[j_max].orbit.period() << " T)\n";

        std::cout << "\n    r2_orbit = " << r2_orbit << ";  v2_orbit = " << v2_orbit
                  << "\n    r2_prop  = " << r2_prop  << ";  v2_prop  = " << v2_prop
                  << "\n    delta_r2 = " << delta_r2 << ";  delta_v2 = " << delta_v2
                  << std::endl;

        std::cout << "\nFailures: " << static_cast<double>(failures) / N * 100 << " %";
    }
}
