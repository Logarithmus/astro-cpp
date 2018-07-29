#include "acos_test.hpp"

#include <random>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

namespace astro_cpp::tests
{
    void acos_test(const int N)
    {
        std::random_device rand_device;
        std::mt19937 rand_gen(rand_device());
        std::uniform_real_distribution<double> cosine(-1.0, 1.0);

        std::vector<double> cosines;
        cosines.reserve(N);
        for (int i = 0; i < N; ++i)
        {
            cosines.emplace_back(cosine(rand_gen));
        }


        //std::acos()
        int i = 0;
        using namespace std::chrono;
        auto t1 = high_resolution_clock::now();
        //while (i < N)
        //{
        //    sum += std::acos(cosines[i++]);
       // }
        auto t2 = high_resolution_clock::now();
       // long long time_std = duration_cast<nanoseconds>(t2 - t1).count();
       // std::cout << sum << std::endl
        //          << "std::acos(): " << time_std << std::endl << std::endl
       //          << static_cast<double>(time_std) / N << " ns / acos" << std::endl << std::endl;


        //fast_acos()
        i = 0;
        double sum = 0.0;
//        double max_abs_error = 0.0;
        t1 = high_resolution_clock::now();
        while (i < N)
        {
            //max_abs_error = std::max(max_abs_error, std::abs(fast_acos_arora_russel(cosines[i]) - std::acos(cosines[i])));
            //++i;
            sum += fast_acos_nvidia(cosines[i++]);
        }
        t2 = high_resolution_clock::now();
        long long time_fast = duration_cast<nanoseconds>(t2 - t1).count();
        std::cout << "Max abs error: " << sum << std::endl
                  << "fast_acos(): " << time_fast << std::endl << std::endl
                  << static_cast<double>(time_fast) / N << " ns / acos" << std::endl << std::endl;
    }
} // namespace astro_cpp::tests
