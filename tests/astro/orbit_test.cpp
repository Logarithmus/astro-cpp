#include "orbit_test.hpp"

#include <iostream>
#include <cstdio>
#include <random>
#include <chrono>
#include <vector>
#include "astro/orbit.hpp"
#include "math/const.hpp"

namespace astro_cpp::tests
{
    double kepler_solver_landis(double M, double e)
    {
        // Solving cubic equation based on Pade approximant for sin(E) (more info in the paper)
        const double alpha = (3.0 * math::PI_SQR + 1.6 * math::PI * (math::PI - std::abs(M)) / (1.0 + e)) / (math::PI_SQR - 6),
                     d = 3.0 * (1.0 - e) + alpha * e,
                     M_sqr = M * M,
                     alpha_d = alpha * d,
                     q = 2.0 * alpha_d * (1.0 - e) - M_sqr,
                     q_sqr = q * q,
                     r = 3.0 * alpha_d * (d - 1.0 + e) * M + M_sqr * M,
                     pre_w = std::abs(r) + std::sqrt(q_sqr*q + r*r),
                     w = std::cbrt(pre_w * pre_w);

        // Finally get our initial guess
        double E = ((2.0 * r * w) / (w*w + w*q + q_sqr) + M) / d;

        // 3rd order Householder's method
        double sin_E = std::sin(E),
               cos_E = std::cos(E),
               f = E - e * sin_E - M,
               Df = 1.0 - e * cos_E,
               DDf = e * sin_E,
               DDDf = e * cos_E,
               Df_Df = Df * Df,
               f_f = f * f;
        E = E - (6.0 * f * Df_Df - 3.0 * f_f * DDf) /
                (6.0 * Df_Df * Df - 6.0 * f * Df * DDf + f_f * DDDf);

        return E;
    }

    double kepler_solver_newton(double M, double e)
    {
        double E = ((M < 0) ? M - e : M + e),
               E0;
        int iters = 0;
        do
        {
            ++iters;
            E0 = E;
            E = E0 - (E0 - e * std::sin(E0) - M) / (1 - e * std::cos(E0));
        } while ((std::abs(E - E0) > 1e-14) && (iters < 15));

        return E;
    }

    void kepler_solver_test(const int N)
    {
        //double M[N], e[N];
        std::vector<double> M, e;
        M.reserve(N);
        e.reserve(N);
        std::random_device rand_device;
        std::mt19937 rand_gen(rand_device());
        std::uniform_real_distribution M_rand(-10 * math::TAU, 10 * math::TAU),
                                       e_rand(0.0, 1.0);

        for (int i = 0; i < N; ++i)
        {
            M.emplace_back(M_rand(rand_gen));
            e.emplace_back(e_rand(rand_gen));
        }

        double sum = 0.0;
        using namespace std::chrono;
        std::cout.precision(15);
        auto t1 = high_resolution_clock::now();
        for (int i = 0; i < N; ++i)
        {
            M[i] = std::fmod(M[i], math::TAU);
            if (M[i] > math::PI)
            {
                M[i] -= math::TAU;
            }
            else if (M[i] < -math::PI)
            {
                M[i] += math::TAU;
            }

            double E = kepler_solver_landis(M[i], e[i]),
                   delta = M[i] - (E - e[i] * std::sin(E));
           // std::cout << M[i] << std::endl;
            sum += M[i];
            if (delta > 1e-14)
            {
                std::cout << "M = " << M[i]
                      << "\ne = " << e[i]
                      << "\nE = " << E
                      << "\ndelta = " << delta << std::endl;
            }
        }
        auto t2 = high_resolution_clock::now();
        long long delta_t = duration_cast<nanoseconds>(t2 - t1).count();
        std::cout << "Mean = " << sum / N
                  << "\nTotal time: " << delta_t << " ns\n"
                  << static_cast<double>(delta_t) / N << " ns per operation\n";
    }
} // namespace astro_cpp::tests
